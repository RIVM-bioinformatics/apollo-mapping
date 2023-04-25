rule multiqc:
    input:
        expand(
            OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.zip",
            sample=SAMPLES,
            read="R1 R2".split(),
        ),
        expand(
            OUT + "/clean_fastq/{sample}_fastp.json",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/qc_mapping/samtools_stats/{sample}_metrics.txt",
            sample=SAMPLES,
            read="R1 R2".split(),
        ),
        expand(
            OUT + "/qc_mapping/insertsize/{sample}_metrics.txt",
            sample=SAMPLES,
            read="R1 R2".split(),
        ),
        # expand(
        #     OUT + "/qc_mapping/allele_frequency/{sample}.tsv",
        #     sample=SAMPLES,
        #     read="R1 R2".split(),
        # ),
        expand(
            OUT + "/identify_species/{sample}/{sample}_bracken_species.kreport2",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/variants/evaluation/{sample}.txt",
            sample=SAMPLES,
        ),
        # OUT + "/qc_de_novo_assembly/quast/report.tsv",
        # OUT + "/qc_de_novo_assembly/checkm/checkm_report.tsv",
        # expand(OUT + "/log/clean_fastq/clean_fastq_{sample}.log", sample=SAMPLES),

    output:
        OUT + "/multiqc/multiqc.html",
        phred=OUT + "/multiqc/multiqc_data/multiqc_data.json",
        seq_len=OUT + "/multiqc/multiqc_data/multiqc_fastqc.txt",
    message:
        "Generating multiqc report"
    log:
        OUT + "/log/multiqc/multiqc.log",
    threads: config["threads"]["multiqc"]
    conda:
        "../envs/multiqc.yaml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
    params:
        config_file= "config/multiqc_config.yaml",
        output_dir= OUT + "/multiqc",
    resources:
        mem_gb = config["mem_gb"]["multiqc"],
    shell:
        """
        multiqc --interactive --force --config {params.config_file} \
        -o {params.output_dir} \
        -n multiqc.html {input} &> {log}
        """