
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


rule index_bam:
    input:
        bam=OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam",
    output:
        bai=OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam.bai",
    message:
        "Indexing bam file of {wildcards.sample} (mapped on '{wildcards.ref_type}')"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    log:
        OUT + "/log/index_bam/{sample}__{ref_type}.log",
    shell:
        """
samtools index {input.bam} 2>&1>{log}
        """
