
# This now calls variants with ploidy 2
# Otherwise, it is not able to detect heterozygosity
# Run this standard with ploidy 2? Or do another check prior to this to check heteroyzgysogity?
rule call_variants:
    input:
        bam = OUT + "/mapped_reads/duprem/{sample}.bam",
        ref = OUT + "/reference/reference.fasta",
        bai = OUT + "/mapped_reads/duprem/{sample}.bam.bai",
        gatk_dict = OUT + "/reference/reference.dict",
        samtools_index = OUT + "/reference/reference.fasta.fai"
    output:
        vcf = OUT + "/variants/raw/{sample}.vcf"
    message: "Calling variants for {wildcards.sample}"
    params:
        ploidy = 2
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/log/haplotypecaller/{sample}.log"
    threads:
        config["threads"]["call_variants"]
    resources:
        mem_gb = config["mem_gb"]["call_variants"]
    shell:
        """
gatk HaplotypeCaller \
--input {input.bam} \
--output {output.vcf} \
--reference {input.ref} \
--sample-ploidy {params.ploidy} 2>&1>{log}
        """

rule mark_low_confidence_variants:
    input:
        vcf = OUT + "/variants/raw/{sample}.vcf",
        ref = OUT + "/reference/reference.fasta",
    output:
        vcf = OUT + "/variants/marked/{sample}.vcf"
    message: "Marking low confidence variants for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/log/mark_low_confidence_variants/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk VariantFiltration \
-R {input.ref} \
-V {input.vcf} \
-O {output.vcf} \
--filter-name "depth_min_10" \
--filter-expression "DP < 10" \
--filter-name "fs_max_60" \
--filter-expression "FS > 60.0" \
--filter-name "varconf_min_2" \
--filter-expression "QD < 2.0" \
--filter-name "mapqual_min_40" \
--filter-expression "MQ < 40.0" 2>&1>{log}
        """

rule remove_low_confidence_variants:
    input:
        vcf = OUT + "/variants/marked/{sample}.vcf",
        ref = OUT + "/reference/reference.fasta",
    output:
        vcf = OUT + "/variants/{sample}.vcf"
    message: "Remove low confidence variants for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/log/remove_low_confidence_variants/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk SelectVariants \
-V {input.vcf} \
-O {output.vcf} \
-R {input.ref} \
--exclude-filtered 2>&1>{log}
        """

rule extract_snps_only:
    input:
        vcf = OUT + "/variants/{sample}.vcf",
        ref = OUT + "/reference/reference.fasta",
    output:
        vcf = OUT + "/variants/snps/{sample}.snps.vcf",
    message: "Keeping only SNPs for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.4.0.0"
    log:
        OUT + "/log/extract_snps_only/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk SelectVariants \
-V {input.vcf} \
-O {output.vcf} \
-R {input.ref} \
--select-type-to-include SNP 2>&1>{log}
        """

rule variant_evaluation:
    input:
        vcf = OUT + "/variants/marked/{sample}.vcf",
        ref = OUT + "/reference/reference.fasta",
    output:
        eval = OUT + "/variants/evaluation/{sample}.txt",
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