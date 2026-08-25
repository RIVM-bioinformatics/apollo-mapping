# Imports
import yaml as YAML
from typing import Dict

from snakemake.rules import Rule

def apply_bcftools_rule_defaults(rule_obj:Rule) -> None:
    """ DRY: extend a rule that relies on mason with its defaults; can't define 'conda:' in base rule """
    # TODO: realize identical code exists in estimate_freebayes_qual.smk
    #       Refactor into generic "workflow/helpers.py"
    rule_obj.threads = config["threads"]["other"]
    rule_obj.resources = {"mem_gb": config["mem_gb"]["other"]}
    rule_obj.conda = "../envs/bcftools.yaml"
    rule_obj.container = "" #"docker://PATH/TO/MASON/CONTAINER:<PINNEDVERSION>"

def get_filter_thresholds(wildcards, input:list) -> Dict[str,float]:
    """ Read """
    # sadly "input" is a list, so we need to explicitly rely on the order
    # TODO: pitch in the "read" fitted data yamls,
    #       and check if data was fitted succesfully or not.
    #       If not, fallback to priors
    vcf, fname_yaml_cov, fname_yaml_qual = input[0:3]

    class AttrDict(dict):
        """ mixin between dictionary and object, with both options for key lookup """
        # TODO: refactor / move elsewhere!
        # https://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute
        def __init__(self, *args, **kwargs):
            super(AttrDict,self).__init__(*args,**kwargs)
            self.__dict__ = self

    input_dict = {key: getattr(input,key) for key in input.keys()}
    yaml_cov = YAML.safe_load(open(fname_yaml_cov))
    yaml_qual = YAML.safe_load(open(fname_yaml_qual))
    mapped_on = list(set(fname_yaml_cov.split("_")).intersection(['species','strain']))[0]
    ppf_th_cov = config['filter_variants_'+mapped_on]
    d = AttrDict({
        "min_QUAL": yaml_qual['ppf_table'][float(ppf_th_cov['min-ppf-QUAL'])],
        "max_QUAL": yaml_qual['ppf_table'][float(ppf_th_cov['max-ppf-QUAL'])],
        "min_DP": yaml_cov['ppf_table'][float(ppf_th_cov['min-ppf-DP'])],
        "max_DP": yaml_cov['ppf_table'][float(ppf_th_cov['max-ppf-DP'])],
    })
    return d


def get_blacklist_bed(wildcards) -> str:
    """ get the full path to the blacklist for variant filtering """
    checkpoint_path = checkpoints.assign_reference.get(sample=wildcards.sample).output.reference

    # TODO: make sure blacklist filename is DYNAMICALLY obtained
    #       At this moment there is only a single multispecies, so current hard-coded return will work just fine.
    #       In rules.filter_variants.run it is checked (config["trigger_multiclade_masking_workflow"])
    #       if the bedfile is being used or is silently omitted from usage.
    return f"{OUT}/reference/species/cauris-GCA_002759435.3-blacklist.bed"


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
    message:
        "Filtering variants for {wildcards.sample} on {wildcards.ref_type} based on modelled thresholds (and blacklisting)",
    log:
        OUT + "/log/filter_variants/" + rule_name + "{sample}-on-{ref_type}.log",
    params:
        # TODO: add skipping/stearing into apollo_mapping.py
        omit_qual_thresholding = config.get("omit_qual_thresholding",False),
        omit_cov_thresholding = config.get("omit_cov_thresholding",False),
        ths = get_filter_thresholds
    run:
        # gather extra header lines
        extra_header_lines = []

        # Always filter indels and complex.
        steps = [
            "bcftools filter -e 'TYPE=\"indel\"' -s Indel -m+ {input.vcf}",
            "bcftools filter -e 'INFO/TYPE=\"complex\"' -s Complex -m+"
        ]

        # In 'species" mode AND multiclade, filter on the blacklist
        # TODO: !!important!! at this point the pipeline will crash in some conditions
        #       e.g. non C.auris ;-) e.g. not-multiclade.
        #       Need to get fixed ASAP ... but can't be done without extra input
        #       Probably it's needed to add to input: ref=MultiReferenceProvider.get_ref_path
        #       Based on this, we can check if it's linked to a blacklist
        if config["trigger_multiclade_masking_workflow"] == "True" and wildcards.ref_type == "species":
            header_line = '##FILTER=<ID=Blacklisted,Description="Variant overlaps with softclipped and/or segmental variable regions">'
            steps.append(
                """bcftools annotate -a <(cat {input.bed} | awk '{{ OFS=\"\\t\"; print $1,$2+1,$3,\"Blacklisted\" }}') """
                f"-c CHROM,FROM,TO,FILTER "
                f'--header-lines <(echo \'{header_line}\')'
            )

        # Conditionally (default=yes) filter on quality
        if not params.omit_qual_thresholding:
            header_qual = '##FILTER=<ID=QualFit,Description="Variant quality (QUAL) is outside the fitted model thresholds">'
            extra_header_lines.append(header_qual)
            steps.append("bcftools filter -e 'QUAL < {params.ths.min_QUAL} || QUAL > {params.ths.max_QUAL}' -s QualFit -m+ "
                )

        # Conditionally (default=yes) filter on coverage
        if not params.omit_cov_thresholding:
            header_cov = '##FILTER=<ID=CovFit,Description="Variant coverage (INFO/DP) is outside the fitted model thresholds">'
            extra_header_lines.append(header_cov)
            steps.append("bcftools filter -e 'INFO/DP < {params.ths.min_DP} || INFO/DP > {params.ths.max_DP}' -s CovFit -m+ "
                )

        if extra_header_lines:
            # add extra header lines; better done using (blanco) annotate, not bcftools reheader
            header_lines = "\\n".join(extra_header_lines)
            steps.append(
                f"bcftools annotate --header-lines <(echo -e '{header_lines}')"
            )

        ## Add write-to-output-file to the final command, merge CLI commands and add stderr logfile writing
        steps[-1] += " -o {output.vcf}"
        full_command = " | ".join(steps)
        full_command += " > {log} 2>&1"

        # Execute the full command!
        shell(full_command)

rule_name="bcftools_simplify_fails"
rule bcftools_simplify_fails:
    input:
        vcf = OUT + "/variants/filtered/{sample}_on_{ref_type}-labeled.vcf",
    output:
        vcf = OUT + "/variants/filtered/{sample}_on_{ref_type}-simplified.vcf",
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