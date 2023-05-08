rule samtools_stats:
    input:
        bam = OUT + "/mapped_reads/duprem/{sample}.bam",
    output:
        txt = OUT + "/qc_mapping/samtools_stats/{sample}_metrics.txt",
    message: "Calculating samtools stats for {wildcards.sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/logs/samtools_stats/{sample}.log"
    threads:
        config["threads"]["samtools_stats"]
    resources:
        mem_gb = config["mem_gb"]["samtools_stats"]
    shell:
        """
samtools stats {input.bam} > {output.txt} 2>{log}
        """


rule get_insert_size:
    input:
        bam = OUT + "/mapped_reads/duprem/{sample}.bam",
    output:
        txt = OUT + "/qc_mapping/insertsize/{sample}_metrics.txt",
        pdf = OUT + "/qc_mapping/insertsize/{sample}_report.pdf",
    message: "Calculating insert size for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    log:
        OUT + "/logs/get_insert_size/{sample}.log"
    threads:
        config["threads"]["picard"]
    resources:
        mem_gb = config["mem_gb"]["picard"]
    shell:
        """
java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
I={input.bam} \
O={output.txt} \
H={output.pdf} 2>&1>{log}
        """

rule get_filter_status:
    input:
        vcf = OUT + "/variants/marked/{sample}.vcf",
    output:
        tsv = OUT + "/qc_mapping/get_filter_status/{sample}.tsv"
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
        expand(OUT + "/qc_mapping/get_filter_status/{sample}.tsv", sample = SAMPLES)
    output:
        OUT + "/qc_mapping/report_filter_status.tsv"
    message: "Combining variant QC reports"
    log:
        OUT + "/logs/combine_filter_status.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb = config["mem_gb"]["other"]
    shell:
        """
python workflow/scripts/combine_variant_tables.py --input {input} --output {output} --fields FILTER
        """

rule count_allelefreq_multiallelic:
    input:
        vcf = OUT + "/variants/{sample}.vcf",
    output:
        tsv = OUT + "/qc_mapping/count_allelefreq_multiallelic/{sample}.tsv"
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
        expand(OUT + "/qc_mapping/count_allelefreq_multiallelic/{sample}.tsv", sample = SAMPLES)
    output:
        OUT + "/qc_mapping/report_allelefreq_multiallelic.tsv"
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