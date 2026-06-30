rule make_coverage_track:
    input:
        OUT + "/mapped_reads/duprem/{sample}.bam",
    output:
        bigwig=OUT + "/variants/coverage_tracks/{sample}.bw",
    message:
        "Generating coverage track for {wildcards.sample}."
    conda:
        "../envs/deeptools.yaml"
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0"
    threads: config["threads"]["deeptools"]
    resources:
        mem_gb=config["mem_gb"]["deeptools"],
    log:
        OUT + "/log/filter_variants/make_coverage_track_{sample}.log",
    shell:
        """
        bamCoverage -b {input} -o {output} 2>&1>{log}
        """

# rule curvefitting_to_determine_filter_values:
#     input:
#         bigwig=OUT + "/variants/coverage_tracks/{sample}.bw
#         vcf=OUT + "/variants/raw/{sample}.vcf"
#     output:
#         filter_values=OUT + "/variants/coverage_tracks/{sample}_filter_values.txt",
#     message:
#         "Performing curve fitting for {wildcards.sample} to determine coverage filter values.",
#     conda:
#         "../envs/curve_fitting.yaml"
#     threads: config["threads"]["curve_fitting"]
#     resources:
#         mem_gb=config["mem_gb"]["curve_fitting"],
#     log:
#         OUT + "/log/filter_variants/curve_fitting_{sample}.log",
#     shell:
#         """python workflow/scripts/curve_fitting.py -b {input.bigwig} -v {input.vcf} -o {output.filter_values} 2>&1>{log}"""
        


# rule filter_variants:
#     input:
#         vcf=OUT + "/variants/raw/{sample}.vcf",
#         filter_values=OUT + "/variants/coverage_tracks/{sample}_filter_values.txt",
#     output:
#         vcf=OUT + "/variants/filtered/{sample}.vcf",
#     message:
#         "Filtering variants for {wildcards.sample} based on coverage.",
#     conda:
#         "../envs/bcftools.yaml"
#     threads: config["threads"]["filter_variants"]
#     resources:
#         mem_gb=config["mem_gb"]["filter_variants"],
#     log:
#         OUT + "/log/filter_variants/filter_variants_{sample}.log",
#     shell:
#         """python workflow/scripts/filter_variants.py -i {input.vcf} -f {input.filter_values} -o {output.vcf} 2>&1>{log}"""

# Legacy rule from apollo snp
# rule extract_snps_only:
#     input:
#         vcf=OUT + "/variants/{sample}.vcf",
#         ref=OUT + "/reference/reference.fasta",
#     output:
#         vcf=OUT + "/variants/snps/{sample}.snps.vcf",
#     message:
#         "Keeping only SNPs for {wildcards.sample}"
#     conda:
#         "../envs/gatk_picard.yaml"
#     container:
#         "docker://broadinstitute/gatk:4.4.0.0"
#     log:
#         OUT + "/log/extract_snps_only/{sample}.log",
#     threads: config["threads"]["filter_variants"]
#     resources:
#         mem_gb=config["mem_gb"]["filter_variants"],
#     shell:
#         """
# gatk SelectVariants \
# -V {input.vcf} \
# -O {output.vcf} \
# -R {input.ref} \
# --select-type-to-include SNP 2>&1>{log}
#         """