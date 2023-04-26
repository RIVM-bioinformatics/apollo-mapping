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

rule pileup_contig_metrics:
    input:
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
    output:
        summary=OUT
        + "/qc_mapping/bbtools/per_sample/{sample}_MinLenFiltSummary.tsv",
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
            OUT
            + "/qc_mapping/bbtools/per_sample/{sample}_MinLenFiltSummary.tsv",
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

rule extract_allele_frequencies:
    input:
        vcf = OUT + "/variants/{sample}.vcf",
    output:
        tsv = OUT + "/qc_mapping/allele_frequency/{sample}.tsv"
    message: "Writing allele frequency of variants to table for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/logs/extract_allele_frequencies/{sample}.log"
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
-O {output.tsv} 2>&1>{log}
        """

# rule get_filter_reasons_variants:
#     input:
