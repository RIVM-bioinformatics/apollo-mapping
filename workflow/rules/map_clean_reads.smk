
rule bwa_mem:
    input:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        ref=MultiReferenceProvider.get_ref_path,
        # "fake" input needed for DAG construction (and bwa index files should exist)
        #idx = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".sa",
        idx_bwa = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".sa",
        idx_fal = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".fal",
        idx_fai = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".fai",
        idx_bed = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".bed",
    output:
        #sam=temp(OUT + "/mapped_reads/raw/{ref_type}/{sample}.sam"),
        sam=OUT + "/mapped_reads/raw/{ref_type}/{sample}.sam",
    log:
        OUT + "/log/bwa_mem/{sample}__{ref_type}.log",
    message:
        "Mapping reads for {wildcards.sample} (on '{wildcards.ref_type}')"
    params:
        bases_per_batch="100000000",
        verbosity="3",
        softclip_supp_aln="-Y",
        rgid="'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'",
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/bwa:0.7.17"
    threads:
        config["threads"]["bwa"]
    resources:
        mem_gb=config["mem_gb"]["bwa"],
    shell:
        """
bwa mem \
-K {params.bases_per_batch} \
-v {params.verbosity} \
-t {threads} \
{params.softclip_supp_aln} \
-R {params.rgid} \
{input.ref} \
{input.r1} {input.r2} 2>{log} 1>{output}
        """

rule sam_to_sorted_bam:
    input:
        sam=OUT + "/mapped_reads/raw/{ref_type}/{sample}.sam",
    output:
        bam=OUT + "/mapped_reads/sorted/{ref_type}/{sample}.bam",
    message:
        "Convert sam to sorted bam for {wildcards.sample} (mapped on '{wildcards.ref_type}')"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/sam_to_sorted_bam/{sample}__{ref_type}.log",
    threads:
        config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
samtools view -b -@ {threads} {input.sam} 2>{log} | \
samtools sort -@ {threads} - 1> {output.bam} 2>>{log}
        """


rule MarkDuplicates:
    input:
        OUT + "/mapped_reads/sorted/{ref_type}/{sample}.bam",
    output:
        bam=OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam",
        metrics=OUT + "/mapped_reads/duprem/{ref_type}/{sample}.metrics",
    message:
        "Marking and removing optical duplicates for {wildcards.sample} (mapped on '{wildcards.ref_type}')"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    params:
        use_singularity=config["use_singularity"],
    threads:
        config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    log:
        OUT + "/log/MarkDuplicates/{sample}__{ref_type}.log",
    shell:
        """
if [ {params.use_singularity} == True ]
then
    EXEC=\"java -jar /usr/picard/picard.jar\"
else
    EXEC=picard
fi

$EXEC \
MarkDuplicates \
INPUT={input} \
OUTPUT={output.bam} \
METRICS_FILE={output.metrics} \
VALIDATION_STRINGENCY=SILENT \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
ASSUME_SORT_ORDER="queryname" \
CLEAR_DT="false" \
ADD_PG_TAG_TO_READS=false 2>&1>{log}
        """


rule base_index_bam:
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
samtools index {input.bam} 2>&1>{log}
        """

use rule base_index_bam as index_bam with:
    log:
        OUT + "/log/index_bam/{sample}__{ref_type}.log",
    input:
        bam=rules.MarkDuplicates.output.bam,
    output:
        bai=str(rules.MarkDuplicates.output.bam)+".bai",
    message:
        "Indexing bam file of {wildcards.sample} (mapped on '{wildcards.ref_type}')"


rule bam_bifurcate_accessions:
    input:
        bam=rules.MarkDuplicates.output.bam,
        # required DAG order trigger
        _bai = rules.index_bam.output.bai,
        # Mind the input is actually <fa>".accessions.yaml";
        # doing MultiReferenceProvider.get_ref_path+".accessions.yaml" gives
        # TypeError unsupported operand type(s) for +: 'method' and 'str'
        fa=MultiReferenceProvider.get_ref_path
    output:
        bam_nuclear=OUT + "/mapped_reads/final/{ref_type}-{sample}-nuclear.bam",
        bam_mitochondrial=OUT + "/mapped_reads/final/{ref_type}-{sample}-mitochondrial.bam",
    message:
        "Filtering out mitochondrial accessions in BAM of {wildcards.sample} (mapped on '{wildcards.ref_type}')"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    params:
        yaml = lambda wildcards, input: f"{input.fa}.accessions.yaml"
    threads:
        1
    log:
        OUT + "/log/bam_mito_filtered/{sample}__{ref_type}.log",
    shell:
        """
accs_nuclear=$(cat {params.yaml}  | yq ".nuclear" | sed 's/^- //' | tr "\n" " " | sed 's/ $//');
accs_mitochondrial=$(cat {params.yaml}  | yq ".mitochondrial" | sed 's/^- //' | tr "\n" " " | sed 's/ $//');
# usefull to place in the logs ...
printf "accs_nuclear      : '$accs_nuclear'\n"
printf "accs_mitochondrial: '$accs_mitochondrial'\n"
if [ -f {output.bam_nuclear} ]; then
    # already copied;
    ok=1
elif [ "$accs_mitochondrial" == "" ]; then
    cp -p {input.bam} {output.bam_nuclear};
    # empty mitochondrial bam file; read the comment in the else statement
    samtools view -@ 1 -bH {input.bam} > {output.bam_mitochondrial};
else
    # !important! realize "^@<accessions>" are unchanged!
    # Possibele solution, if desired is to use "samtools reheader <in.header.sam> <in.bam>"
    # In this/our use-case, removing the reads withou header modification will do.
    samtools view -@ 1 -bh {input.bam} $accs_nuclear > {output.bam_nuclear};
    samtools view -@ 1 -bh {input.bam} $accs_mitochondrial > {output.bam_mitochondrial};
fi
        """


use rule index_bam as index_bam_nuclear with:
    log:
        OUT + "/log/index_bam/{sample}__{ref_type}-genomic.log",
    input:
        bam = rules.bam_bifurcate_accessions.output.bam_nuclear
    output:
        bai = str(rules.bam_bifurcate_accessions.output.bam_nuclear) + ".bai"

use rule index_bam as index_bam_mitochondrial with:
    log:
        OUT + "/log/index_bam/{sample}__{ref_type}-mitochondrial.log",
    input:
        bam = rules.bam_bifurcate_accessions.output.bam_mitochondrial
    output:
        bai = str(rules.bam_bifurcate_accessions.output.bam_mitochondrial) + ".bai"
