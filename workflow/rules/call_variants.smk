# This now calls variants with ploidy 2
# Otherwise, it is not able to detect heterozygosity
# Run this standard with ploidy 2? Or do another check prior to this to check heteroyzgysogity?
rule call_variants:
    input:
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
        ref=OUT + "/reference/reference.fasta",
        bai=OUT + "/mapped_reads/duprem/{sample}.bam.bai",
        #gatk_dict=OUT + "/reference/reference.dict",
        #samtools_index=OUT + "/reference/reference.fasta.fai",
    output:
        vcf=OUT + "/variants/raw/{sample}.vcf",
    message:
        "Calling variants for {wildcards.sample}",
    params:
        ploidy=config["Freebayes"]["ploidy"],
        haplotype_length=config["Freebayes"]["haplotype_length"],
        min_repeat_size=config["Freebayes"]["min_repeat_size"], 
        min_repeat_entropy=config["Freebayes"]["min_repeat_entropy"], 
    conda:
        "../envs/freebayes.yaml"
    container:
        "docker://quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0"
    log:
        OUT + "/log/freebayes/{sample}.log",
    threads: config["threads"]["call_variants"]
    resources:
        mem_gb=config["mem_gb"]["call_variants"],
    shell:
        """
freebayes \
-f {input.ref} \
-p {params.ploidy} \
--haplotype-length {params.haplotype_length} \
--min-repeat-size {params.min_repeat_size} \
--min-repeat-entropy {params.min_repeat_entropy} \
{input.bam} >{output.vcf} 2>{log}
        """

