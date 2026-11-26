# Imports
import yaml as YAML
from typing import Dict
from snakemake.rules import Rule

# DRY: env name variables
CONDA_ENV_BCF = "../envs/bcftools.yaml"
CONDA_ENV_DEEPTOOLS = "../envs/deeptools.yaml"
CONDA_ENV_CURVE_FITTING = "../envs/curve_fitting.yaml"
#CONTAINER_BCF = "docker://staphb/bcftools:1.16" #"docker://PATH/TO/BCF/CONTAINER:<PINNEDVERSION>"

def apply_bcftools_rule_defaults(rule_obj:Rule) -> None:
    """ DRY: extend a rule that relies on mason with its defaults; can't define 'conda:' in base rule """
    # TODO: realize identical code exists in estimate_freebayes_qual.smk
    #       Refactor into generic "workflow/helpers.py"
    rule_obj.threads = config["threads"]["other"]
    rule_obj.resources = {"mem_gb": config["mem_gb"]["other"]}

def apply_deeptool_rule_defaults(rule_obj:Rule) -> None:
    """ DRY: extend a rule that relies on mason with its defaults; can't define 'conda:' in base rule """
    # TODO: realize identical code exists in estimate_freebayes_qual.smk
    #       Refactor into generic "workflow/helpers.py"
    rule_obj.threads = config["threads"]["other"]
    rule_obj.resources = {"mem_gb": config["mem_gb"]["other"]}

def get_filter_thresholds(wildcards, input:list) -> Dict[str,float]:
    """ Read """
    # sadly "input" is a list, so we need to explicitly rely on the order
    # TODO for Sofie:   pitch in the "read" fitted data yamls,
    #                   and check if data was fitted succesfully or not.
    #                   If not, fallback to priors
    vcf, fname_yaml_cov, fname_yaml_qual = input[0:3]

    yaml_cov = YAML.safe_load(open(fname_yaml_cov))
    yaml_qual = YAML.safe_load(open(fname_yaml_qual))
    mapped_on = list(set(fname_yaml_cov.split("_")).intersection(['species','strain']))[0]
    ppf_th_cov = config['filter_variants_'+mapped_on]

    d = {
        "min_QUAL": yaml_qual['ppf_table'][float(ppf_th_cov['min-ppf-QUAL'])],
        "max_QUAL": yaml_qual['ppf_table'][float(ppf_th_cov['max-ppf-QUAL'])],
        "min_DP": yaml_cov['ppf_table'][float(ppf_th_cov['min-ppf-DP'])],
        "max_DP": yaml_cov['ppf_table'][float(ppf_th_cov['max-ppf-DP'])],
    }
    return d


def get_blacklist_bed(wildcards) -> str:
    """ get the full path to the blacklist for variant filtering """
    checkpoint_path = checkpoints.assign_reference.get(sample=wildcards.sample).output.reference

    # TODO: make sure blacklist filename is DYNAMICALLY obtained
    #       At this moment there is only a single multispecies, so current hard-coded return will work just fine.
    #       In rules.filter_variants.run it is checked (config["trigger_multiclade_masking_workflow"])
    #       if the bedfile is being used or is silently omitted from usage.
    return f"{OUT}/reference/species/cauris-GCA_002759435.3-blacklist.bed"

apply_deeptools_rule_defaults(rules.make_coverage_track)

rule_name="make_coverage_track"
rule make_coverage_track:
    input:
        bam= rules.bam_bifurcate_accessions.output.bam_nuclear
    output:
        bigwig=OUT + "/variants/coverage_tracks/{sample}.bw",
    message:
        "Generating coverage track for {wildcards.sample}."
    conda:
        CONDA_ENV_DEEPTOOLS
    log:
        OUT + "/log/filter_variants/make_coverage_track_{sample}.log",
    shell:
        """
        bamCoverage -b {input} -o {output} 2>&1>{log}
        """

rule_name = "curve_fitting"
rule curve_fitting:
    input:
        vcf=rules.vcflib_breakmulti_wave.output.vcf,
        bigwig=rules.make_coverage_track.output.bigwig,
    output:
        yaml_cov_qual_fitted=OUT + "/curve_fitting/fitted_ppf_table_{sample}.yaml",
        plot_cov=OUT + "/curve_fitting/{sample}_model_fittings_coverage.html",
        plot_qual=OUT + "/curve_fitting/{sample}_model_fittings_score.html",
    message:
        "Performing curve fitting on score and coverage values for {wildcards.sample}"
    conda:
        CONDA_ENV_CURVE_FITTING
    log:
        OUT + "/curve_fitting/curve_fitting_{sample}.log",
    shell:
        """
        python workflow/scripts/curve_fitting.py \
            -v {input.vcf} \
            -b {input.bigwig} \
            -o {output.yaml_cov_qual_fitted} \
            --plot_cov {output.plot_cov} \
            --plot_qual {output.plot_qual} \
            > {log} 2>&1
        """

rule_name="filter_variants"
rule filter_variants:
    input:
        vcf=rules.vcflib_breakmulti_wave.output.vcf,
        # TODO: need to pitch in the real fitted data
        ##yaml_cov_fitted=rules.curvefitting_on_coverage.output.yaml,
        ##yaml_qual_fitted=rules.curvefitting_on_quality_score.output.yaml,
        yaml_cov_prior=rules.est_cov_distribution_specs.output,
        yaml_qual_prior=rules.est_qual_distribution_specs.output,
        bed=get_blacklist_bed,
    output:
        vcf=OUT + "/variants/filtered/{sample}_on_{ref_type}-labeled.vcf",
    conda:
        CONDA_ENV_BCF
    message:
        "Filtering variants for {wildcards.sample} on {wildcards.ref_type} based on modelled thresholds (and blacklisting)",
    log:
        OUT + "/log/filter_variants/" + rule_name + "{sample}-on-{ref_type}.log",
    params:
        # TODO: add skipping/stearing into apollo_mapping.py
        omit_qual_thresholding = config.get("omit_qual_thresholding",False),
        omit_cov_thresholding = config.get("omit_cov_thresholding",False),
        ths = get_filter_thresholds,
        # overtake parameter from config (which isn't known to script: directive
        trigger_multiclade_masking_workflow = str(config.get("trigger_multiclade_masking_workflow",False)),
    script:
        "../scripts/filter_variants.py"

rule_name="bcftools_simplify_fails"
rule bcftools_simplify_fails:
    input:
        vcf = OUT + "/variants/filtered/{sample}_on_{ref_type}-labeled.vcf",
    output:
        vcf = OUT + "/variants/filtered/{sample}_on_{ref_type}-simplified.vcf",
    conda:
        CONDA_ENV_BCF
    log:
        OUT + "/log/filter_variants/" + rule_name + "{sample}-on-{ref_type}.log",
    shell:
        """
        bcftools filter -e 'FILTER != "PASS" && FILTER != "."' -s FAIL {input.vcf} -o {output.vcf}
        """

apply_bcftools_rule_defaults(rules.bcftools_simplify_fails)

rule_name="bcftools_only_pass"
rule bcftools_only_pass:
    input:
        vcf = OUT + "/variants/filtered/{sample}_on_{ref_type}-simplified.vcf",
    output:
        vcf = OUT + "/variants/filtered/{sample}_on_{ref_type}-passed.vcf",
    conda:
        CONDA_ENV_BCF
    log:
        OUT + "/log/filter_variants/" + rule_name + "{sample}-on-{ref_type}.log",
    shell:
        """
        bcftools view -f PASS -v snps {input.vcf} -o {output.vcf}
        """

apply_bcftools_rule_defaults(rules.bcftools_only_pass)


# Legacy rule from apollo snp
# rule extract_snps_only:
#     input:
#         vcf=OUT + "/variants/{sample}.vcf",
#         ref=OUT + "/reference/reference.fasta",
#     output:
#         vcf=OUT + "/variants/snps/{sample}.snps.vcf",
#     message:
#         "Keeping only SNPs for {wildcards.sample}"
#     conda:
#         "../envs/gatk_picard.yaml"
#     container:
#         "docker://broadinstitute/gatk:4.4.0.0"
#     log:
#         OUT + "/log/extract_snps_only/{sample}.log",
#     threads: config["threads"]["filter_variants"]
#     resources:
#         mem_gb=config["mem_gb"]["filter_variants"],
#     shell:
#         """
# gatk SelectVariants \
# -V {input.vcf} \
# -O {output.vcf} \
# -R {input.ref} \
# --select-type-to-include SNP 2>&1>{log}
#         """
