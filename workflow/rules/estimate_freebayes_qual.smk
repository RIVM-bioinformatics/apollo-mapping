# DRY and concise: reuse existing rules to guarantee identical mapping + variant calling settings
include: "map_clean_reads.smk"
include: "generate_reference_indices.smk"
include: "call_variants.smk"

# Python Imports
from pathlib import Path

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

# ------------------------------------------------------------------------------------------------ #
# base rules "applied by function" (can't define 'conda:' in base rule is "use rule X as Y with:")
# ------------------------------------------------------------------------------------------------ #

from snakemake.rules import Rule

CONDA_ENV_MASON = "../envs/mason.yaml"
CONDA_ENV_MASON = "../envs/mason.yaml"

def apply_mason_rule_defaults(rule_obj:Rule) -> None:
    """ DRY: extend a rule that relies on mason with its defaults; can't define 'conda:' in base rule """
    rule_obj.threads = config["threads"]["other"]
    rule_obj.resources = {"mem_gb": config["mem_gb"]["other"]}
    rule_obj.conda = CONDA_ENV_MASON
    rule_obj.container = "" #"docker://PATH/TO/MASON/CONTAINER:<PINNEDVERSION>"

def apply_bcftools_rule_defaults(rule_obj:Rule) -> None:
    """ DRY: extend a rule that relies on mason with its defaults; can't define 'conda:' in base rule """
    rule_obj.threads = config["threads"]["other"]
    rule_obj.resources = {"mem_gb": config["mem_gb"]["other"]}
    rule_obj.conda = "../envs/bcftools.yaml"
    rule_obj.container = "" #"docker://PATH/TO/MASON/CONTAINER:<PINNEDVERSION>"

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
    log:
        OUT + "/log/simulated/" + rule_name + ".log"
    conda:
        # !important! snakemake inspects rules upon initialization for which conda envs to build
        #             so, in the first "mason" rule, specify it explicitly.
        CONDA_ENV_MASON
    shell:
        """
        # create a 100kb genome, consisting of 2 scaffolds
        mason_genome -l {params.length_scaf} -l {params.length_scaf} -o {output};
        sed 's/^>/>scaf/' {output} > {output}.patched;
        mv {output}.patched {output};
        """

apply_mason_rule_defaults(rules.mason_genome)

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

apply_mason_rule_defaults(rules.mason_variator)


# Helper functions to parse required, single-value parameters from larger input file

def _parse_picard_metrics_block(filepath:Path,at_max_num_lines:int=100) -> dict:
    """ parse the ## METRICS CLASS block from a Picard result file """
    with filepath.open(mode="r",encoding="utf-8") as f:
        lineId: int = 0
        while True:
            line = f.readline().strip()
            lineId+=1
            if lineId >= at_max_num_lines:
                break
            elif line.startswith("## METRICS CLASS"):
                keys = f.readline().strip().split("\t")
                vals = f.readline().strip().split("\t")
                break
    data = dict(zip(keys,vals))
    return data

def get_dynamic_param_read_length(wildcards, input):
    """ """
    data = _parse_picard_metrics_block(Path(input._data_for_read_length))
    return int(float(data['MEAN_ALIGNED_READ_LENGTH']))

def get_dynamic_param_read_depth(wildcards, input):
    """ """
    data = _parse_picard_metrics_block(Path(input._data_for_read_depth))
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
        read_length=get_dynamic_param_read_length,
        depth=get_dynamic_param_read_depth
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

apply_mason_rule_defaults(rules.mason_simulator)

# ----------------------------------------------------------------------------------------
# required indices (mostly inherited rules)
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

rule_name="bcftools_bgzip_index_vcf"
rule bcftools_bgzip_index_vcf:
    input:
        SIMULATED_GENOME_PREFIX + "-mutated.vcf"
    output:
        vcf_gz = SIMULATED_GENOME_PREFIX + "-mutated.vcf.gz",
        vcf_csi = SIMULATED_GENOME_PREFIX+ "-mutated.vcf.gz.csi"
    log:
        OUT + "/log/simulated/" + rule_name + ".log"
    message:
        "Bgzip and index VCF file: "+ SIMULATED_GENOME_FNAME
    shell:
        """
        bcftools view  {input} -O z -o {output.vcf_gz}; bcftools index {output.vcf_gz};
        """

apply_bcftools_rule_defaults(rules.bcftools_bgzip_index_vcf)

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

rule_name="est_qual_freebayes_bgzip_index_vcf"
use rule bcftools_bgzip_index_vcf as est_qual_freebayes_bgzip_index_vcf with:
    input:
        OUT + "/simulated/data/{sample}_on_{ref_type}.vcf"
    output:
        vcf_gz = OUT + "/simulated/data/{sample}_on_{ref_type}.vcf.gz",
        vcf_csi = OUT + "/simulated/data/{sample}_on_{ref_type}.vcf.gz.csi"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"
    message:
        "Bgzip and index VCF file from {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"

rule_name="est_qual_filter_vcf"
rule est_qual_filter_vcf:
    input:
        vcf_freebayes = OUT + "/simulated/data/{sample}_on_{ref_type}.vcf.gz",
        vcf_mutations = SIMULATED_GENOME_PREFIX + "-mutated.vcf.gz",
    output:
        OUT + "/simulated/data/{sample}_on_{ref_type}_filtered.vcf"
    log:
        OUT + "/log/simulated/{sample}-on-{ref_type}-" + rule_name + ".log"
    message:
        "Filter variants using bcftools {wildcards.sample} on '{wildcards.ref_type}' ["+rule_name+"]"
    shell:
        ## realize this filter that "discards highly unlikely SNPs" is not stringent enough
        #bcftools filter  -i 'INFO/AF >= 0.10' {input} | bcftools view -v snps  -o  {output}
        ## realize that there is occurence of False Negative SNPs:
        ## - occurring in  input.vcf_mutations
        ## - not called in input.vcf_freebayes
        ## - in the few tested cases so far, FN% can be 1-2%, mostly in less deep sequenced dataset
        """
        bcftools view -T {input.vcf_mutations} {input.vcf_freebayes} | bcftools view -v snps -o {output}
        """

apply_bcftools_rule_defaults(rules.est_qual_filter_vcf)


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
    conda:
        "../envs/python_pandas_env.yaml"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
        ## don't use bcftools view -H {input}; saves the dependancy of bcftools + python
        grep -v "^#" {input}  | cut -f 6 | python workflow/scripts/array2nsmmmmsd.py >  {output}
        """
