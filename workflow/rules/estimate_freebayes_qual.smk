# DRY and concise: reuse existing rules to guarantee identical mapping + variant calling settings
include: "map_clean_reads.smk"
include: "generate_reference_indices.smk"
include: "call_variants.smk"

# Variables concerning the simulated genome;
# In case you adjust these, please make sure all divisions yield integers, not floats!
SIMULATED_GENOME_NT_SIZE=100000
SIMULATED_GENOME_NUM_SCAFS=2
SIMULATED_GENOME_SCAF_SIZE=int(SIMULATED_GENOME_NT_SIZE/SIMULATED_GENOME_NUM_SCAFS)
SIMULATED_GENOME_SCAF_SIZE_KB=int(SIMULATED_GENOME_SCAF_SIZE/1000)

# DRY: reused prefix for all output files (genome-2x50kb)
_sfx = str(SIMULATED_GENOME_NUM_SCAFS)+"x"+str(SIMULATED_GENOME_SCAF_SIZE_KB)+"kb"
SIMULATED_GENOME_PREFIX = OUT + "/simulated/genome-"+ _sfx
SIMULATED_GENOME_FNAME = SIMULATED_GENOME_PREFIX + ".fa"

# ----------------------------------------------------------------------------------------
# SeQan / mason: generate a simulated minified genome and simulated PE reads from it
# ----------------------------------------------------------------------------------------

rule_name="mason_genome"
rule mason_genome:
    output:
        SIMULATED_GENOME_FNAME
    params:
        length_genome=SIMULATED_GENOME_NT_SIZE,
        length_scaf=SIMULATED_GENOME_SCAF_SIZE,
    conda:
        "../envs/mason.yaml"
    log:
        OUT + "/log/simulated/" + rule_name + ".log"
    shell:
        """
        # create a 100kb genome, consisting of 2 scaffolds
        mason_genome -l {params.length_scaf} -l {params.length_scaf} -o {output};
        sed 's/^>/>scaf/' {output} > {output}.patched;
        mv {output}.patched {output};
        """

rule_name="mason_variator"
rule mason_variator:
    input:
        SIMULATED_GENOME_FNAME
    output:
        fa=SIMULATED_GENOME_PREFIX + "-mutated.fa",
        vcf=SIMULATED_GENOME_PREFIX + "-mutated.vcf"
    params:
        snp_rate=0.001
    conda:
        "../envs/mason.yaml"
    message:
        "introduce random variation in 100kb simulated genome"
    log:
        OUT + "/log/simulated/" + rule_name + ".log"
    shell:
        """
        mason_variator \
          -ir {input} \
          -ov {output.vcf} \
          -of {output.fa} \
          --snp-rate {params.snp_rate}
        """


# Helper functions to parse required, single-value parameters from larger input file

from pathlib import Path
def _parse_metrics_block(filepath:Path) -> dict:
    """ parse the ## METRICS CLASS block from a Picard result file """
    with filepath.open(mode="r",encoding="utf-8") as f:
        for line in f.readlines():
            if line.startswith("## METRICS CLASS"):
                keys = f.readline().strip().split("\t")
                vals = f.readline().strip().split("\t")
                break
    data = dict(zip(keys,vals))
    return data


def get_dynamic_param_read_length(wildcards, input):
    """ """
    data = _parse_metrics_block(input._data_for_read_length)
    return int(data['MEAN_ALIGNED_READ_LENGTH'])

    with open(input._data_for_read_length, "r") as f:
        # parse METRICS CLASS block
        for line in f.readlines():
            if line.startswith("## METRICS CLASS"):
                keys = f.readline().strip().split("\t")
                vals = f.readline().strip().split("\t")
                break
        data = dict(zip(keys,vals))
        return int(data['MEAN_ALIGNED_READ_LENGTH'])

def get_dynamic_param_read_depth(wildcards, input):
    """ """
    data = _parse_metrics_block(input._data_for_read_length)
    return int(data['MEDIAN_COVERAGE'])


rule_name="mason_simulator"
rule mason_simulator:
    input:
        fa=SIMULATED_GENOME_FNAME,
        vcf=SIMULATED_GENOME_PREFIX + "-mutated.vcf",
        _data_for_read_depth = OUT + "/qc_mapping/CollectWgsMetrics/{sample}__{ref_type}.txt",
        _data_for_read_length = OUT + "/qc_mapping/CollectAlignmentSummaryMetrics/{sample}__{ref_type}.txt"
    output:
        R1 = OUT + "/simulated/data/{sample}_on_{ref_type}_reads_L1.fq",
        R2 = OUT + "/simulated/data/{sample}_on_{ref_type}_reads_R2.fq"
    params:
        length_genome=SIMULATED_GENOME_NT_SIZE,
        read_length=150, #get_dynamic_param_read_length,
        depth=50, #get_dynamic_param_read_depth,
    conda:
        "../envs/mason.yaml"
    message:
        "simulate PE fastq data for {wildcards.sample} on '{wildcards.ref_type}' from {input.fa}"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"
    shell:
        """
        N=$(echo {params.length_genome} {params.read_length} {params.depth} | awk '{{ print int($1/(2*$2)*$3) }}');
        mason_simulator \
            -n $N \
            --seq-technology illumina \
            --seq-mate-orientation FR \
            --illumina-read-length {params.read_length} \
            -o  {output.R1} \
            -or {output.R2} \
            -ir {input.fa} \
            -iv {input.vcf};
        """

# ----------------------------------------------------------------------------------------
# required indices (inherited rules only)
# ----------------------------------------------------------------------------------------

rule_name="est_qual_bwa_index"
use rule bwa_index_ref as est_qual_bwa_index_ref with:
    # !important! only overrule required directives; keep all others
    input:
        SIMULATED_GENOME_FNAME
    output:
        SIMULATED_GENOME_FNAME+".sa"
    message:
        "Indexing reference genome [bwa index]: "+ SIMULATED_GENOME_FNAME
    log:
        OUT + "/log/simulated/" + rule_name + ".log"

rule_name="est_qual_samtools_faidx_ref"
use rule samtools_faidx_ref as est_qual_samtools_faidx_ref with:
    # !important! only overrule required directives; keep all others
    input:
        SIMULATED_GENOME_FNAME
    output:
        SIMULATED_GENOME_FNAME+".fai"
    log:
        OUT + "/log/simulated/" + rule_name + ".log"
    message:
        "Indexing reference genome [samtools faidx]: "+ SIMULATED_GENOME_FNAME

# ----------------------------------------------------------------------------------------
# per-sample analyses
# ----------------------------------------------------------------------------------------

rule_name="est_qual_bwa_mem"
use rule bwa_mem as est_qual_bwa_mem with:
    # !important! only overrule required directives; keep all others
    input:
        r1=OUT + "/simulated/data/{sample}_on_{ref_type}_reads_L1.fq",
        r2=OUT + "/simulated/data/{sample}_on_{ref_type}_reads_R2.fq",
        ref=SIMULATED_GENOME_FNAME,
        _trigger_dag_bwa_index=SIMULATED_GENOME_FNAME + ".sa",
        _trigger_dag_samtools_faidx=SIMULATED_GENOME_FNAME + ".fai"
    output:
        sam=OUT + "/simulated/data/{sample}_on_{ref_type}.sam"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"

rule_name="est_qual_sam_to_sorted_bam"
use rule sam_to_sorted_bam as est_qual_sam_to_sorted_bam with:
    # !important! only overrule required directives; keep all others
    input:
        sam=OUT + "/simulated/data/{sample}_on_{ref_type}.sam"
    output:
        bam=OUT + "/simulated/data/{sample}_on_{ref_type}.bam"
    message:
        "Convert sam to sorted bam for {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"


rule_name="est_qual_index_bam"
use rule index_bam as est_qual_index_bam with:
    # !important! only overrule required directives; keep all others
    input:
        bam=OUT + "/simulated/data/{sample}_on_{ref_type}.bam"
    output:
        bai=OUT + "/simulated/data/{sample}_on_{ref_type}.bam.bai"
    message:
        "generate samtools index for {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"

rule_name="est_qual_freebayes"
use rule freebayes as est_qual_freebayes with:
    # !important! only overrule required directives; keep all others
    input:
        bam=OUT + "/simulated/data/{sample}_on_{ref_type}.bam",
        ref=SIMULATED_GENOME_FNAME,
        _trigger_dag_samtools_index=OUT + "/simulated/data/{sample}_on_{ref_type}.bam.bai",
    output:
        vcf=OUT + "/simulated/data/{sample}_on_{ref_type}.vcf"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"
    message:
        "Call variants using freebayes {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"

rule_name="est_qual_filter_vcf"
rule est_qual_filter_vcf:
    input:
        OUT + "/simulated/data/{sample}_on_{ref_type}.vcf"
    output:
        OUT + "/simulated/data/{sample}_on_{ref_type}_filtered.vcf"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"
    message:
        "Filter variants using bcftools {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools filter  -i 'INFO/AF >= 0.10' {input} | bcftools view -v snps  -o  {output}
        """


rule_name="est_qual_median_and_stdv"
rule est_qual_median_and_stdv:
    input:
        OUT + "/simulated/data/{sample}_on_{ref_type}_filtered.vcf"
    output:
        OUT + "/simulated/data/{sample}_on_{ref_type}_priors.yaml"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"
    message:
        "Write median and stdv of VCF quality score to yaml from {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/python_pandas_env.yaml"
    shell:
        """
        # don't use bcftools view -H {input}; saves the dependancy of bcftools + python
        grep -v "^#" {input}  | cut -f 6 | python workflow/scripts/array2nsmmmmsd.py >  {output}
        """
