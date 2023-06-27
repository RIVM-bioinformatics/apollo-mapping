rule clean_fastq:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
    output:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        unpaired=OUT + "/clean_fastq/{sample}_unpaired_joined.fastq.gz",
        html=OUT + "/clean_fastq/{sample}_fastp.html",
        json=OUT + "/clean_fastq/{sample}_fastp.json",
    message:
        "Filtering reads for {wildcards.sample}."
    conda:
        "../envs/qc_and_clean.yaml"
    container:
        "docker://biocontainers/fastp:v0.20.1_cv1"
    threads: config["threads"]["fastp"]
    resources:
        mem_gb=config["mem_gb"]["fastp"],
    log:
        OUT + "/log/clean_fastq/clean_fastq_{sample}.log",
    params:
        mean_quality=config["mean_quality_threshold"],
        window_size=config["window_size"],
        min_length=config["min_read_length"],
    shell:
        """
fastp --in1 {input.r1} \
--in2 {input.r2} \
--out1 {output.r1} \
--out2 {output.r2} \
--unpaired1 {output.unpaired} \
--unpaired2 {output.unpaired} \
--html {output.html} \
--json {output.json} \
--report_title "FastP report for sample {wildcards.sample}" \
--detect_adapter_for_pe \
--thread {threads} \
--cut_right \
--cut_window_size {params.window_size} \
--cut_mean_quality {params.mean_quality} \
--correction \
--length_required {params.min_length} > {log} 2>&1
        """


rule qc_raw_fastq:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][wildcards.read],
    output:
        html=OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.html",
        zip=OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip",
    message:
        "Running FastQC on pre-trimmed reads for {wildcards.sample}."
    conda:
        "../envs/qc_and_clean.yaml"
    container:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: config["threads"]["fastqc"]
    resources:
        mem_gb=config["mem_gb"]["fastqc"],
    params:
        output_dir=OUT + "/qc_raw_fastq/",
    log:
        OUT + "/log/qc_raw_fastq/qc_raw_fastq_{sample}_{read}.log",
    shell:
        """
bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log} > {log} 2>&1
        """


rule qc_clean_fastq:
    input:
        OUT + "/clean_fastq/{sample}_p{read}.fastq.gz",
    output:
        html=OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.html",
        zip=OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.zip",
    message:
        "Running FastQC after filtering/trimming {wildcards.sample}."
    conda:
        "../envs/qc_and_clean.yaml"
    container:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: config["threads"]["fastqc"]
    resources:
        mem_gb=config["mem_gb"]["fastqc"],
    log:
        OUT + "/log/qc_clean_fastq/qc_clean_fastq_{sample}_{read}.log",
    params:
        output_dir=OUT + "/qc_clean_fastq/",
    shell:
        """
        if [ -s {input} ]
        then
            fastqc --quiet --outdir {params.output_dir} {input} >> {log} 
        else  
            touch {output}
        fi
        """
