
# define & imports external module qc_mapping.smk
import builtins
builtins.OUT = OUT
builtins.config = config
builtins.SAMPLES = SAMPLES
module juno_qc_mapping:
    snakefile: "../../juno_mapping/workflow/rules/qc_mapping.smk"

# TODO: move elsewhere!
class AttrDict(dict):
    """ mixin between dictionary and object, with both options for key lookup """
    # https://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

def apolloified_qc_mapping(rule_name:str) -> AttrDict:
    """ DRY helper function to manipulate rule IO from juno to apollo

    Reason is that input in rules in juno_mapping/workflow/rules/qc_mapping.smk is of this format:

        input:
            bam=OUT + "/mapped_reads/duprem/{sample}.bam",
            ref=OUT + "/reference/reference.fasta",

    In multi-species, multiclade aware apollo-mapping, this has to be extended on.

        input:
            bam=OUT + "/mapped_reads/duprem/{sample}__{ref_type}.bam",
            ref = MultiReferenceProvider.get_ref_path,

    As a consequence, log: should co-adopt alongside as well. If not done (correctly) errors like these appear:

        SyntaxError:
        Not all output, log and benchmark files of rule apollo_CollectGcBiasMetrics contain the same wildcards.
        This is crucial though, in order to avoid that two or more jobs write to the same file.

    """
    # TODO: there's this code in ViroConstrictor:
    # https://github.com.mcas.ms/RIVM-bioinformatics/ViroConstrictor/blob/release/version_1.7.0/
    #       ViroConstrictor/workflow/helpers/generic_workflow_methods.py#L252-L301
    # I tested it (20260421) and it doesn't work here.
    # Reason is 99.9% certain snakemake's version (requiring 9.x, apollo now on 7.x)
    # Once this is settled, rule_name variables can be dynamically obtained,
    # mostly likely thereby deprecating the need of rule_name = "rule_name" variable
    # just above each transformed rule (juno->apollo).
    # Dream scenario is that it would work from within **THIS** function.
    # Tricky though, since inspect might get messed up by being called from within this helper function,
    # and not from within the actual rule block.

    return AttrDict( {
        "input":  OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam",
        "log":    OUT + f"/log/{rule_name}/{{sample}}__{{ref_type}}.log",
        "output_dir": str(OUT + f"/qc_mapping/{rule_name}"),
        "output_picard_txt": OUT + f"/qc_mapping/{rule_name}/{{sample}}__{{ref_type}}.txt",
    } )

# shorter variable name
AQC = apolloified_qc_mapping

# define the (re-used/overwritten) rules

rule_name = "insertsize"
use rule get_insert_size from juno_qc_mapping as apollo_get_insert_size with:
    log: AQC(rule_name).log
    input:
        # fake input to trigger the DAG ---> eventual multiqc target in rule all
        # TODO: once rule all multiqc works, we can remove these fake targets here
        _fake1 = OUT + "/qc_mapping/CollectQualityYieldMetrics/{sample}__{ref_type}.txt",
        _fake2 = OUT + "/qc_mapping/CollectWgsMetrics/{sample}__{ref_type}.txt",
        _fake3 = OUT + "/qc_mapping/CollectGcBiasMetrics/{sample}__{ref_type}.txt",
        _fake4 = OUT + "/qc_mapping/CollectAlignmentSummaryMetrics/{sample}__{ref_type}.txt",
        _fake5 = OUT + "/qc_mapping/bbtools/per_sample/{sample}__{ref_type}_MinLenFiltSummary.tsv",
        #_fake6 = OUT + "/qc_mapping/bbtools/bbtools_summary_report.tsv",
        #_fake7 = OUT + "/qc_mapping/bbtools/bbtools_scaffolds.tsv",
        bam = AQC(rule_name).input
    message:
        #"Calculating insert size for {{get_rule_name()}} {wildcards.sample} (mapped on '{wildcards.ref_type}')"
        "Calculating insert size for {wildcards.sample} (mapped on '{wildcards.ref_type}')"
    output:
        txt = str(AQC(rule_name).output_dir) + "/{sample}__{ref_type}_metrics.txt",
        pdf = str(AQC(rule_name).output_dir) + "/{sample}__{ref_type}_report.pdf",

rule_name = "pileup_contig_metrics"
use rule pileup_contig_metrics from juno_qc_mapping as apollo_pileup_contig_metrics with:
    log: AQC(rule_name).log
    input: bam = AQC(rule_name).input
    message: "Making pileup and calculating contig metrics for {wildcards.sample} (mapped on '{wildcards.ref_type}')"
    output:
        summary = OUT + "/qc_mapping/bbtools/per_sample/{sample}__{ref_type}_MinLenFiltSummary.tsv",
        perScaffold = OUT + "/qc_mapping/bbtools/per_sample/{sample}__{ref_type}_perMinLenFiltScaffold.tsv",

rule_name = "parse_bbtools"
use rule parse_bbtools from juno_qc_mapping as apollo_parse_bbtools with:
    input:  OUT + "/qc_mapping/bbtools/per_sample/{sample}__{ref_type}__perMinLenFiltScaffold.tsv",

rule_name = "parse_bbtools_summary"
use rule parse_bbtools_summary from juno_qc_mapping as apollo_parse_bbtools_summary with:
    input:  OUT + "/qc_mapping/bbtools/per_sample/{sample}__{ref_type}__MinLenFiltScaffold.tsv",

# DRY helper snippet; some picard mapping stats subtools require the reference sequence,
# which can differ per sample. 'wilcards' is needed here,
# (which differs from usage when within actual rules)
PICARD_INPUT_INCLUDING_REFERENCE = lambda wildcards: {
    "bam": AQC(rule_name).input,
    "ref": MultiReferenceProvider.get_ref_path(wildcards)
}

rule_name = "CollectQualityYieldMetrics"
use rule CollectQualityYieldMetrics from juno_qc_mapping as apollo_CollectQualityYieldMetrics with:
    log: AQC(rule_name).log
    output: txt = AQC(rule_name).output_picard_txt
    input: bam = AQC(rule_name).input

rule_name = "CollectWgsMetrics"
use rule CollectWgsMetrics from juno_qc_mapping as apollo_CollectWgsMetrics with:
    log: AQC(rule_name).log
    output: txt = AQC(rule_name).output_picard_txt
    input:
        # !important! don't input mapped_reads/duprep; sorted bam required here
        bam = OUT + "/mapped_reads/sorted/{ref_type}/{sample}.bam",
        ref = MultiReferenceProvider.get_ref_path,

rule_name = "CollectAlignmentSummaryMetrics"
use rule CollectAlignmentSummaryMetrics from juno_qc_mapping as apollo_CollectAlignmentSummaryMetrics with:
    log: AQC(rule_name).log
    output: txt = AQC(rule_name).output_picard_txt
    #input: unpack(PICARD_INPUT_INCLUDING_REFERENCE)
    input:
        bam = AQC(rule_name).input,
        ref = MultiReferenceProvider.get_ref_path #(wildcards)


rule_name = "CollectGcBiasMetrics"
use rule CollectGcBiasMetrics from juno_qc_mapping as apollo_CollectGcBiasMetrics with:
    log: AQC(rule_name).log
    #input: unpack(PICARD_INPUT_INCLUDING_REFERENCE)
    input:
        bam = AQC(rule_name).input,
        ref = MultiReferenceProvider.get_ref_path #(wildcards)
    output:
        txt = str(AQC(rule_name).output_dir) + "/{sample}__{ref_type}.txt",
        pdf = str(AQC(rule_name).output_dir) + "/{sample}__{ref_type}.pdf",
        summary = str(AQC(rule_name).output_dir) + "/{sample}__{ref_type}.summary.txt",


if False:

        rule pileup_contig_metrics:
            input:
                bam=OUT + "/mapped_reads/duprem/{sample}.bam",
            output:
                summary=OUT + "/qc_mapping/bbtools/per_sample/{sample}_MinLenFiltSummary.tsv",
                perScaffold=OUT
                + "/qc_mapping/bbtools/per_sample/{sample}_perMinLenFiltScaffold.tsv",
            message:
                "Making pileup and calculating contig metrics for {wildcards.sample}."
            conda:
                "../envs/bbtools.yaml"
            container:
                "docker://staphb/bbtools:39.01"
            log:
                OUT + "/log/qc_mapping/bbtools_{sample}.log",
            threads: config["threads"]["bbtools"]
            resources:
                mem_gb=config["mem_gb"]["bbtools"],
            shell:
                """
                pileup.sh in={input.bam} \
                    out={output.perScaffold} \
                    secondary=f \
                    samstreamer=t 2> {output.summary}
                cp {output.summary} {log}
                """


        rule parse_bbtools:
            input:
                expand(
                    OUT + "/qc_mapping/bbtools/per_sample/{sample}_perMinLenFiltScaffold.tsv",
                    sample=SAMPLES,
                ),
            output:
                OUT + "/qc_mapping/bbtools/bbtools_scaffolds.tsv",
            message:
                "Parsing the results of bbtools (pileup contig metrics)."
            threads: config["threads"]["bbtools"]
            resources:
                mem_gb=config["mem_gb"]["bbtools"],
            log:
                OUT + "/log/qc_mapping/pileup_contig_metrics_combined.log",
            script:
                "../scripts/parse_bbtools.py"


        rule parse_bbtools_summary:
            input:
                expand(
                    OUT + "/qc_mapping/bbtools/per_sample/{sample}_MinLenFiltSummary.tsv",
                    sample=SAMPLES,
                ),
            output:
                OUT + "/qc_mapping/bbtools/bbtools_summary_report.tsv",
            message:
                "Parsing the results of bbtools (pileup contig metrics) and making a multireport."
            threads: config["threads"]["bbtools"]
            resources:
                mem_gb=config["mem_gb"]["bbtools"],
            log:
                OUT + "/log/qc_mapping/pileup_contig_metrics_combined.log",
            shell:
                "python workflow/scripts/parse_bbtools_summary.py -i {input} -o {output} > {log}"


        rule get_insert_size:
            input:
                bam=OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam",
            output:
                txt=OUT + "/qc_mapping/insertsize/{sample}__{ref_type}_metrics.txt",
                pdf=OUT + "/qc_mapping/insertsize/{sample}__{ref_type}_report.pdf",
            log:
                OUT + "/log/get_insert_size/{sample}__{ref_type}.log",
            message:
                "Calculating insert size for {wildcards.sample}"
            conda:
                "../envs/gatk_picard.yaml"
            container:
                "docker://broadinstitute/picard:2.27.5"
            threads:
                config["threads"]["picard"]
            resources:
                mem_gb=config["mem_gb"]["picard"],
            shell:
                """
        java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
        I={input.bam} \
        O={output.txt} \
        H={output.pdf} 2>&1>{log}
                """


        rule CollectAlignmentSummaryMetrics:
            input:
                bam=OUT + "/mapped_reads/duprem/{sample}.bam",
                ref=OUT + "/reference/reference.fasta",
            output:
                txt=OUT + "/qc_mapping/CollectAlignmentSummaryMetrics/{sample}.txt",
            container:
                "docker://broadinstitute/picard:2.27.5"
            log:
                OUT + "/log/CollectAlignmentSummaryMetrics/{sample}.log",
            threads: config["threads"]["picard"]
            resources:
                mem_gb=config["mem_gb"]["picard"],
            shell:
                """
        java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics -I {input.bam} -R {input.ref} -O {output}
                """


        rule CollectWgsMetrics:
            input:
                bam=OUT + "/mapped_reads/duprem/{sample}.bam",
                ref=OUT + "/reference/reference.fasta",
            output:
                txt=OUT + "/qc_mapping/CollectWgsMetrics/{sample}.txt",
            container:
                "docker://broadinstitute/picard:2.27.5"
            log:
                OUT + "/log/CollectWgsMetrics/{sample}.log",
            threads: config["threads"]["picard"]
            resources:
                mem_gb=config["mem_gb"]["picard"],
            shell:
                """
        java -jar /usr/picard/picard.jar CollectWgsMetrics -I {input.bam} -R {input.ref} -O {output}
                """
