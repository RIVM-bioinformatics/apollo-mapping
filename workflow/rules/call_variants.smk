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

rule bcftools_norm:
    input:
        vcf=rules.freebayes.output.vcf,
    output:
        vcf=OUT + "/variants/norm/{ref_type}-{sample}.vcf",
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
        bcftools norm -a -m -any {input.vcf} -o {output.vcf}
        """



if False:

    rule model_vcf_filtering_params:
        input:
            vcf = rules.bcftools_norm.output.vcf,
            bigwig = rules.make_coverage_track.output.bigwig,
            yaml = rules.est_qual_distribution_specs.output,
        output:
            yaml = str(rules.freebayes.output.vcf) + ".modeled-fit-thresholds.yaml"
        message:
            "Performing curve fitting for {wildcards.sample} [on {wildcards.ref_type}] to determine coverage and QUAL filter values."
        params:
            qual_prior_mean = "TODO-OBTAIN-FROM-{input.yaml}",
            qual_prior_stdv = "TODO-OBTAIN-FROM-{input.yaml}"
        shell:
            """python workflow/scripts/curve_fitting.py -b {input.bigwig} -v {input.vcf} -o {output.yaml} 2>&1>{log} """

    rule filter_variants_species:
        wildcard_constraints:
            ref_type="species"
        input:
            vcf = rules.freebayes.output.vcf,
            yaml = rules.model_vcf_filtering_params.output.yaml,
            fa=MultiReferenceProvider.get_ref_path
        output:
            vcf=OUT + "/variants/filtered/{ref_type}/{sample}.vcf",
        params:
            bed = "{input.fa}.variantcalling-blacklist.bed",
        shell:
            """
            custom_script_to_filter_variants.py {input.vcf} --settings {input.yaml} --blacklist {params.bed}
            """

    rule filter_variants_strain:
        wildcard_constraints:
            ref_type="strain"
        input:
            vcf=rules.freebayes.output.vcf,
            yaml=rules.model_vcf_filtering_params.output.yaml,
        output:
            vcf=OUT + "/variants/filtered/{ref_type}/{sample}.vcf",
        shell:
            """
            custom_script_to_filter_variants.py {input.vcf} --settings {input.yaml}
            """

