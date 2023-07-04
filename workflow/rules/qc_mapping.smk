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
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
    output:
        txt=OUT + "/qc_mapping/insertsize/{sample}_metrics.txt",
        pdf=OUT + "/qc_mapping/insertsize/{sample}_report.pdf",
    message:
        "Calculating insert size for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    log:
        OUT + "/log/get_insert_size/{sample}.log",
    threads: config["threads"]["picard"]
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
        bam=OUT + "/mapped_reads/sorted/{sample}.bam",
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
