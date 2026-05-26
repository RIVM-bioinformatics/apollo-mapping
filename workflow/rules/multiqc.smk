
class MultiReferenceProviderMultiQcMappingInput(MultiReferenceProvider):
    """ added class method that dynamically generates multiqc targets, supporting main species & optional strain reference """
    @classmethod
    def get_multiqc_targets(cls, wildcards) -> List[str]:
        """ return a list of 'rule all' mapping specs target results for species (and additional strain-specific mappings) """

        def fill_template_outfiles(sample:str,ref_type:str) -> List[str]:
            outfile_templates = [
                f"{OUT}/qc_mapping/CollectAlignmentSummaryMetrics/{sample}__{ref_type}.txt",
                f"{OUT}/qc_mapping/CollectWgsMetrics/{sample}__{ref_type}.txt"
                ]
            return outfile_templates

        return cls._get_rule_all_targets_given_templates(fill_template_outfiles)


class MultiReferenceProviderMultiQcVariantCallingInput(MultiReferenceProvider):
    """ added class method that dynamically generates multiqc targets, supporting main species & optional strain reference """
    @classmethod
    def get_multiqc_targets(cls, wildcards) -> List[str]:
        """ return a list of 'rule all' variant calling target results for species (and additional strain-specific mappings) """

        def fill_template_outfiles(sample:str,ref_type:str) -> List[str]:
            outfile_templates = [
                #f"{OUT}/qc_mapping/CollectAlignmentSummaryMetrics/{sample}__{ref_type}.txt",
                #f"{OUT}/qc_mapping/CollectWgsMetrics/{sample}__{ref_type}.txt"
                # TODO: debuild to freebayes ouput specs
                #expand(
                #    OUT + "/qc_variant_calling/bcftools_stats/{sample}.txt",
                #    sample=SAMPLES,
                #),
                ## OUT + "/qc_variant_calling/report_filter_status_mqc.tsv",
                #OUT + "/qc_variant_calling/report_allelefreq_mqc.tsv",
                #expand(
                #    OUT + "/qc_mapping/insertsize/{sample}_metrics.txt",
                #    sample=SAMPLES,
                #),
                #expand(
                #    OUT + "/qc_variant_calling/VariantEval/{sample}.txt",
                #    sample=SAMPLES,
                #),
                ]
            return outfile_templates

        return cls._get_rule_all_targets_given_templates(fill_template_outfiles)



rule multiqc_base:
    # basic config for multiqc rules
    # please specify input: output: params: log: for inherited rules,
    # optionally overrule message:
    message:
        "Generating multiqc report"
    threads:
        config["threads"]["multiqc"]
    conda:
        "../envs/multiqc.yaml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
    resources:
        mem_gb=config["mem_gb"]["multiqc"],
    shell:
        """
        multiqc --interactive --force --config {params.config_file} \
        -o {params.output_dir} \
        -n multiqc.html {input} &> {log}
        """

use rule multiqc_base as multiqc_fastq with:
    log:
        OUT + "/log/multiqc/multiqc_fastq.log",
    params:
        config_file="config/multiqc_config.yaml",
        output_dir=OUT + "/multiqc/fastq",
    output:
        OUT + "/multiqc/fastq/multiqc.html",
    input:
        # simple 1:1 relation of sample to report
        expand(
            OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.zip",
            sample=SAMPLES,
            read="R1 R2".split(),
        ),
        expand(
            OUT + "/clean_fastq/{sample}_fastp.json",
            sample=SAMPLES,
        )


use rule multiqc_base as multiqc_mapping with:
    log:
        OUT + "/log/multiqc/mapping/multiqc.log",
    params:
        config_file="config/multiqc_config.yaml",
        output_dir=OUT + "/multiqc/mapping",
    input:
        MultiReferenceProviderMultiQcMappingInput.get_multiqc_targets,
    output:
        OUT + "/multiqc/mapping/multiqc.html"

use rule multiqc_base as multiqc_variantcalling with:
    log:
        OUT + "/log/multiqc/variantcalling/multiqc.log",
    params:
        config_file="config/multiqc_config.yaml",
        output_dir=OUT + "/multiqc/variantcalling",
    input:
        MultiReferenceProviderMultiQcVariantCallingInput.get_multiqc_targets,
    output:
        OUT + "/multiqc/variantcalling/multiqc.html"


use rule multiqc_base as multiqc_faked_all with:
    log:
        OUT + "/log/multiqc/multiqc.log",
    params:
        config_file="config/multiqc_config.yaml",
        output_dir=OUT + "/multiqc",
    input:
        rules.multiqc_fastq.input,
        rules.multiqc_mapping.input,
        rules.multiqc_variantcalling.input,
        # trigger execution of the sub-workflow multiqc reports
        OUT + "/multiqc/fastq/multiqc.html",
        OUT + "/multiqc/mapping/multiqc.html",

        ####OUT + "/multiqc/variantcalling/multiqc.html",

        # TODO: add various other information to multiqc
        #expand(
        #    OUT + "/identify_impurity/{sample}/{sample}_bracken_species.kreport2",
        #    sample=SAMPLES,
        #),
    output:
        OUT + "/multiqc/multiqc.html"
