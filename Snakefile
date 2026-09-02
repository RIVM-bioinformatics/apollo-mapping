import yaml
import os
from typing import List, Dict, Any, Union, Optional
from collections.abc import Callable
from workflow.helpers.generic_workflow_methods import check_rule_environments, check_rule_directives, assign_mem_gb_default
##from workflow.helpers.generic_workflow_methods import get_rule_name

# define paths to workflow & scripts directories (helpers & tools not in conda environment!)
workflow_path = os.path.join(workflow.basedir, "workflow")
scripts_path = os.path.join(workflow_path, "scripts")

# Sample wildcard is constrained to all characters except "/"
# Otherwise files with the same extension in subdirs match as well
wildcard_constraints:
    sample="[^\/]+",
    ref_type="(species|strain)"


sample_sheet = config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

for param in ["threads", "mem_gb"]:
    for k in config[param]:
        config[param][k] = int(config[param][k])

print(SAMPLES)
print("config:",config)

OUT = config["output_dir"]
PATH_TO_REFERENCES = config["reference_genomes_dir"]

# some rules are conditional, or have conditional output files. Capture these for the final rule all
CONDITIONAL_TARGETS = []

localrules:
    all,
    assign_reference,
    copy_reference_species_genomes,
    copy_reference_strain_genomes,
    copy_blacklist_bed,
    # Realize these two below will often trigger this message in Snakemake:
    # "localrules directive specifies rules that are not present in the Snakefile:"
    copy_external_species_genome,
    copy_external_strain_genome


# TODO: refactor/move to a more generic place (other file or even other repo)
class MultiReferenceProvider:
    """ class that based on YSON produced by the checkpoint match_ref orchestrates further scheduling

    Core purpose is to:
        - keep each (PE) fastq connected to its designated reference assembly
        - keep each (PE) fastq connected to its designated reference and strain assembly (if applicable)
        - keep track, given the above, what
    """

    # some config/variables that are best placed here (too)
    SKIP_KRAKEN: bool = config['kraken_db_dir'] == "None"
    SKIP_REFERENCE_SELECTION: bool =  config['skip_reference_selection'] == "True"
    SKIP_MULTICLADE_STRAIN_MAPPING: bool = config['skip_multiclade_strain_mapping'] == "True"
    TRIGGER_MULTICLADE_MASKING_WORKFLOW: bool = config['trigger_multiclade_masking_workflow'] == "True"
    EXTERIOR_FASTA: bool = config['exterior_fasta'] == "True"

    @staticmethod
    def get_yaml_content(sample_id: str):
        """ use the explicit dictionary access to ensure Snakemake tracks the dependency """
        checkpoint_path = checkpoints.assign_reference.get(sample=sample_id).output.reference

        if os.path.exists(checkpoint_path):
            with open(checkpoint_path, "r") as f:
                return yaml.safe_load(f)
        return None

    @classmethod
    def get_species_fasta_path(cls, wildcards):
        """ get local path to the central species reference fitted for this sample (as assigned by match_ref) """
        data = cls.get_yaml_content(wildcards.sample)
        # !important! Explicitly name the wildcard variable (for fasta copying rule)
        reference_species_basename = data["species"]["fasta"]
        return os.path.join(OUT, "reference", "species", data["species"]["fasta"])

    @classmethod
    def get_strain_fasta_path(cls, wildcards):
        """ get local path to clade-specific reference (from a multi-clade group) fitted for this sample (as assigned by match_ref) """
        # We MUST use wildcards directly here as Snakemake passes the object
        data = cls.get_yaml_content(wildcards.sample)
        if data and data.get("strain"):
            # !important! Explicitly name the wildcard variable (for fasta copying rule)
            reference_strain_basename = data["strain"]["fasta"]
            # Ensure the path matches the copy rule (e.g., 'strains')
            return os.path.join(OUT, "reference", "strains", data["strain"]["fasta"])
        return []


    @classmethod
    def get_elegiable_strain_targets(cls, wildcards) -> List[str]:
        """ return a list of rule all target results for samples with an additional strain-specific reference

        Accepts 'wildcards' so Snakemake can call it as a deferred input function.
        """
        output_targets = []
        # We still use the global SAMPLES list here
        for s in SAMPLES:
            try:
                # Use the internal helper that handles the checkpoint check
                data = cls.get_yaml_content(s)
                if data and data.get("strain"):
                    output_targets.append(
                        os.path.join(OUT, f"reference/{s}__strain__just_checking_if_it_works.txt")
                    )
            except (snakemake.exceptions.IncompleteCheckpointException, Exception):
                continue
        return output_targets

    @classmethod
    def get_ref_path(cls, wildcards) -> str:
        """ return the absolute path to 'species' or 'strain' reference fasta """
        data = cls.get_yaml_content(wildcards.sample)
        if wildcards.ref_type == "species":
            return os.path.join(OUT, "reference", "species", data["species"]["fasta"])
        elif wildcards.ref_type == "strain" and data.get("strain"):
            return os.path.join(OUT, "reference", "strains", data["strain"]["fasta"])
        else:
            # !important!
            return ""


    @classmethod
    def _get_rule_all_targets_given_templates(cls,fill_template_outfiles_function:Callable) -> List[str]:
        """ return a list of 'rule all' target results for species and additional strain-specific mappings """
        targets = []
        for sample in SAMPLES:
            # all samples are mapped to their corresponding "species" reference
            # !important! this is the crucial ref_type wildcard variable to discriminate species/strains
            ref_type = "species"
            targets.extend(fill_template_outfiles_function(sample,ref_type))
            try:
                # Take into account if checkpoint has been passed already
                data = cls.get_yaml_content(sample)
                if data.get("strain"):
                    ref_type = "strain"
                    targets.extend(fill_template_outfiles_function(sample,ref_type))

            except (snakemake.exceptions.IncompleteCheckpointException, KeyError):
                # Snakemake will re-evaluate after checkpoint has been passed
                continue
        return targets

    @classmethod
    def get_rule_all_targets(cls, wildcards) -> List[str]:
        """ return a list of 'rule all' target results for species and additional strain-specific mappings """
        def fill_template_outfiles(sample:str,ref_type:str):
            # TODO: once pipeline development has finished,
            #       only the truely final exterior leaves of the DAG needs to get specified
            #outfile_template = f"{OUT}/qc_mapping/insertsize/{sample}__{ref_type}_metrics.txt"
            #return [ outfile_template ]
            outfile_templates = [
                f"{OUT}/qc_mapping/insertsize/{sample}__{ref_type}_metrics.txt",
                f"{OUT}/simulated/data/{sample}_on_{ref_type}_priors_QUAL.yaml",
                f"{OUT}/simulated/data/{sample}_on_{ref_type}_priors_COV.yaml",
                f"{OUT}/mapped_reads/final/{ref_type}-{sample}-nuclear.bam.bai",
                f"{OUT}/mapped_reads/final/{ref_type}-{sample}-mitochondrial.bam.bai",
                f"{OUT}/variants/filtered/{sample}_on_{ref_type}-passed.vcf",
            ]
            return outfile_templates

        return cls._get_rule_all_targets_given_templates(fill_template_outfiles)

include: "workflow/rules/identify_species.smk"
include: "workflow/rules/generate_reference_indices.smk"
# !important! realize rule itself is configured to act on --skip-kraken
include: "workflow/rules/identify_impurity_using_kraken.smk"

if True:
    # fastq cleaning & summarizing stats on it
    include: "workflow/rules/clean_fastq.smk"
else:
    # TODO: dissect rule sort_paired_fastq in juno-mapping, we don't need neither want it in apollo-mapping
    # TODO: once done, apollo-mapping can inherit these rules from juno-mapping
    # fastq cleaning & summarizing stats on it
    include: "juno_mapping/workflow/rules/fastqc_raw_data.smk"
    include: "juno_mapping/workflow/rules/fastqc_clean_data.smk"
    include: "juno_mapping/workflow/rules/clean_fastq.smk"

# ------------------------------------------------------------------------------------------------------ #
# !important!   Realize that from this point onwards, the DAG will branch (in a multi-clade scenario).
#               In will map on "species" and it will map on "clade", which is handled by
#               the helper class MultiReferenceProvider.
# ------------------------------------------------------------------------------------------------------ #

# mapping part of the pipeline
include: "workflow/rules/map_clean_reads.smk"

if MultiReferenceProvider.TRIGGER_MULTICLADE_MASKING_WORKFLOW:
    """ trigger the multiclade masking workflow; adds some targets to rule_all """
    class MultiReferenceProviderConcatSoftclippedInput(MultiReferenceProvider):
        """ added class method that dynamically generates input for per-clade softclip concatenation """
        @classmethod
        def get_rule_all_targets(cls, wildcards) -> List[str]:
            """ extend this classmethod take these additional rule all targets """
            targets = super().get_rule_all_targets(wildcards)
            def fill_template_outfiles(sample:str,ref_type:str) -> List[str]:
                outfile_templates = [
                    f"{OUT}/softclipped/multiclade.{ref_type}-softclipped.bg",
                    f"{OUT}/softclipped-{ref_type}/{sample}.{ref_type}.access-softclipped.bw",
                    f"{OUT}/softclipped-{ref_type}/{sample}.{ref_type}.softclipped.bw",
                    f"{OUT}/coverage/{ref_type}/{sample}.genome.bw",
                    ]
                return outfile_templates
            return targets + cls._get_rule_all_targets_given_templates(fill_template_outfiles)

    # overrule MultiReferenceProvider to become the augmented one
    MultiReferenceProvider = MultiReferenceProviderConcatSoftclippedInput
    include: "workflow/rules/identify_softclipped_regions.smk"

# !imporant! rule all is needed BEFORE:
# - include workflow/rules/qc_mapping.smk, which imports juno-mapping/workflow/rules/qc_mapping.smk
# - definition of PREFIXES

rule all:
    input:
        #expand(OUT + "/reference/{sample}-match-ref-reference.tsv", sample=SAMPLES),
        #expand(OUT + "/reference/{sample}-match-ref-taxid.tsv", sample=SAMPLES),
        #expand(OUT + "/reference/{sample}-samtools-stats.tsv", sample=SAMPLES),
        #expand(OUT + "/reference/{sample}-references.yml", sample=SAMPLES),
        MultiReferenceProvider.get_rule_all_targets,

        # temporarily set fastqc pipeline to output targets
        expand(OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.html", sample=SAMPLES.keys(), read=["R1", "R2"]),
        expand(OUT + "/clean_fastq/{sample}_fastp.json", sample=SAMPLES.keys()),
        OUT + "/multiqc/multiqc.html",

# EOF mapping part of the pipeline with (multi)QC
include: "workflow/rules/qc_mapping.smk"

# variant calling part of the pipeline
# !important! realize MultiReferenceProvider in case of multi-clade analyses
include: "workflow/rules/estimate_freebayes_qual.smk"
include: "workflow/rules/call_variants.smk"
include: "workflow/rules/filter_variants.smk"

# completion part of the pipeline
include: "workflow/rules/multiqc.smk"


# EOF workflow construction: QC
check_rule_environments(workflow=workflow,env_type="conda")     # preferably strict=True
check_rule_environments(workflow=workflow,env_type="container") # preferably strict=True
check_rule_directives(
    workflow=workflow,
    directives=["message", "threads", "resources.mem_gb", "log"],
    strict=False
)
assign_mem_gb_default(workflow=workflow,default_mem_gb=config['mem_gb']['other'])

if False:
    # TODO: IDEA/CONCEPT: pre-define "prefixes" to be used for log files, output (sub)directories etc.
    #       There are other ways of doing thus using get_rule_name(),
    #       but the POC code below does a similar job

    def clean_snakefile_origin(smk_path:str,joiner:str="__") -> str:
        """ Convert the *.smk file path to its de-suffixed name; snakefile itself returns empty string """
        if workflow.main_snakefile == smk_path:
            return ""
        elif smk_path.lower() == "snakefile":
            return ""
        return os.path.basename(smk_path).replace('.smk', '')+joiner

    # 2. Gebruik een dictionary-benadering voor je paden
    # Plaats dit ONDERAAN je Snakefile (na alle includes!)
    # Herschrijf naar AttrDict: PREFIX.rule_name.outdir, PREFIX.rule_name.logdir, PREFIX.rule_name.module_name
    PREFIXES = {
        r.name: f"results/{clean_snakefile_origin(r.snakefile)}{r.name}"
        for r in workflow.rules
    }

    print("PREFIXES:",PREFIXES)

