rule get_filter_status:
    input:
        vcf = OUT + "/variants/marked/{sample}.vcf",
    output:
        tsv = OUT + "/qc_variant_calling/get_filter_status/{sample}.tsv"
    message: "Writing filter status of variants to table for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/logs/get_filter_status/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk VariantsToTable -V {input.vcf} \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F DP \
-F FILTER \
--show-filtered \
-O {output.tsv} 2>&1>{log}
        """

rule combine_filter_status:
    input:
        tsv=expand(
            OUT + "/qc_variant_calling/get_filter_status/{sample}.tsv", sample=SAMPLES
        ),
        header="files/report_filter_status.mqc",
    output:
        OUT + "/qc_variant_calling/report_filter_status_mqc.tsv",
    message: "Combining variant QC reports"
    log:
        OUT + "/logs/combine_filter_status.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb = config["mem_gb"]["other"]
    shell:
        """
python workflow/scripts/combine_variant_tables.py --input {input.tsv} --output {output} --fields FILTER --mqc {input.header}
        """

rule count_allelefreq_multiallelic:
    input:
        vcf = OUT + "/variants/{sample}.vcf",
    output:
        tsv = OUT + "/qc_variant_calling/count_allelefreq_multiallelic/{sample}.tsv"
    message: "Writing allele frequency and multi-allelic status of variants to table for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/logs/count_allelefreq_multiallelic/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk VariantsToTable -V {input.vcf} \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F DP \
-F AF \
-F MULTI-ALLELIC \
-O {output.tsv} 2>&1>{log}
        """

rule combine_allelefreq_multiallelic:
    input:
        expand(OUT + "/qc_variant_calling/count_allelefreq_multiallelic/{sample}.tsv", sample = SAMPLES)
    output:
        OUT + "/qc_variant_calling/report_allelefreq_multiallelic.tsv"
    message: "Combining variant QC reports"
    log:
        OUT + "/logs/combine_allelefreq_multiallelic.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb = config["mem_gb"]["other"]
    shell:
        """
python workflow/scripts/combine_variant_tables.py --input {input} --output {output} --fields AF MULTI-ALLELIC
        """

rule bcftools_stats:
    input:
        ref=OUT + "/reference/reference.fasta",
        vcf=OUT + "/variants/{sample}.vcf",
    output:
        txt=OUT + "/qc_variant_calling/bcftools_stats/{sample}.txt",
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    log:
        OUT + "/log/bcftools_stats/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
bcftools stats \
--fasta-ref {input.ref} \
{input.vcf} \
1>{output.txt} \
2>{log}
        """


rule variant_evaluation:
    input:
        vcf = OUT + "/variants/marked/{sample}.vcf",
        ref = OUT + "/reference/reference.fasta",
    output:
        eval = OUT + "/qc_variant_calling/VariantEval/{sample}.txt",
    message: "Variant evaluation for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/log/variant_evaluation/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk VariantEval \
-R {input.ref} \
-O {output.eval} \
--eval {input.vcf}\
        """