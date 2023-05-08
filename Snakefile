import yaml

# Sample wildcard is constrained to all characters except "/"
# Otherwise files with the same extension in subdirs match as well
wildcard_constraints:
   sample = '[^\/]+'

sample_sheet=config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

for param in ["threads", "mem_gb"]:
    for k in config[param]:
        config[param][k] = int(config[param][k])

# print(SAMPLES)
# print(config)

OUT = config["output_dir"]

localrules:
    all, copy_ref

include: "workflow/rules/identify_species.smk"
include: "workflow/rules/clean_fastq.smk"
include: "workflow/rules/map_clean_reads.smk"
include: "workflow/rules/qc_mapping.smk"
include: "workflow/rules/call_variants.smk"
include: "workflow/rules/multiqc.smk"

rule all:
    input:
        expand(OUT + "/variants/{sample}.vcf", sample = SAMPLES),
        expand(OUT + "/variants/snps/{sample}.snps.vcf", sample = SAMPLES),
        OUT + "/qc_mapping/bbtools/bbtools_scaffolds.tsv",
        OUT + "/qc_mapping/bbtools/bbtools_summary_report.tsv",
        OUT + "/multiqc/multiqc.html",
        OUT + "/qc_mapping/report_allelefreq_multiallelic.tsv",
        OUT + "/qc_mapping/report_filter_status.tsv"
