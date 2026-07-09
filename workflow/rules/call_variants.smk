# This now calls variants with ploidy 2
# Otherwise, it is not able to detect heterozygosity
# Run this standard with ploidy 2? Or do another check prior to this to check heteroyzgysogity?
#from numpy.f2py import rules

rule freebayes:
    input:
        #bam=OUT + "/mapped_reads/duprem/{sample}.bam",
        #bai=OUT + "/mapped_reads/duprem/{sample}.bam.bai",
        bam= rules.bam_bifurcate_accessions.output.bam_nuclear,
        _bai= rules.index_bam_nuclear.output.bai,
        ref=MultiReferenceProvider.get_ref_path,
        #gatk_dict=OUT + "/reference/reference.dict",
        #samtools_index=OUT + "/reference/reference.fasta.fai",
    output:
        vcf=OUT + "/variants/raw/{ref_type}-{sample}.vcf",
    message:
        "Calling variants for {wildcards.sample} [mapped on ref_type'{wildcards.ref_type}']",
    params:
        ploidy=config["freebayes"]["ploidy"],
        haplotype_length=config["freebayes"]["haplotype_length"],
        min_repeat_size=config["freebayes"]["min_repeat_size"],
        min_repeat_entropy=config["freebayes"]["min_repeat_entropy"],
    conda:
        "../envs/freebayes.yaml"
    container:
        "docker://quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0"
    log:
        OUT + "/log/freebayes/{ref_type}-{sample}.log",
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

# Discredited bcfools norm. There are other tools that do a better job!
rule bcftools_norm:
    input:
        vcf=rules.freebayes.output.vcf,
        ref=MultiReferenceProvider.get_ref_path,
    output:
        vcf=OUT + "/variants/norm/bcftools-norm-{ref_type}-{sample}.vcf",
    message:
        "Normalizing called variants for {wildcards.sample} [mapped on ref_type'{wildcards.ref_type}']",
    conda:
        "../envs/bcftools.yaml"
    log:
        OUT + "/log/bcftools-norm/{ref_type}-{sample}.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
        bcftools norm -f {input.ref} -a -m -any {input.vcf} -o {output.vcf}
        """

# Now "normalize" variation landscape.
# This is an example of variant types AFTER vcfbreakmulti from raw freebayes results:
# Especially the "mnp" category we want to get rid of, which vcflib does excellently!
#  22034 TYPE=complex
#   3455 TYPE=del
#   3333 TYPE=ins
#   1878 TYPE=mnp
# 131590 TYPE=snp
rule vcflib_breakmulti_wave:
    input:
        vcf=rules.freebayes.output.vcf,
    output:
        vcf=OUT + "/variants/norm/{ref_type}-{sample}.vcf",
    message:
        "Normalizing called variants for {wildcards.sample} [mapped on ref_type'{wildcards.ref_type}']",
    conda:
        "../envs/vcflib.yaml"
    log:
        OUT + "/log/vcflib_breakmulti_allelicprimitives/{ref_type}-{sample}.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        ## DEPRECATED! vcfwave handles fixing types natively
        # vcfbreakmulti {input.vcf} | vcfallelicprimitives -k -g | workflow/scripts/fix_types.sh > {output.vcf} 2> {log}
        r"""
        vcfbreakmulti {input.vcf} | vcfwave | vcfallelicprimitives \
            | awk '{{gsub(/TYPE=snp,snp/, "TYPE=snp"); gsub(/TYPE=ins,ins/, "TYPE=ins"); gsub(/TYPE=del,del/, "TYPE=del"); print}}' \
            > {output.vcf} 2> {log}
        """
