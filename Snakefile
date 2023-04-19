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
        expand(OUT + "/identify_species/{sample}/{sample}_species_content.txt", sample = SAMPLES),
        expand(OUT + "/qc_mapping/allele_frequency/{sample}.tsv", sample = SAMPLES),
        expand(OUT + "/qc_mapping/samtools_stats/{sample}_metrics.txt", sample = SAMPLES),
        expand(OUT + "/qc_mapping/insertsize/{sample}_metrics.txt", sample = SAMPLES),
        expand(OUT + "/qc_clean_fastq/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["pR1", "pR2"],),
        expand(OUT + "/multiqc/multiqc.html"),
