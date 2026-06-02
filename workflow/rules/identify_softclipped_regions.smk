""" rules that analyse QC'ed mapped BAM files for regions having exhibiting access softclipping on the (central) reference genome fasta """

# A general note on wildcard_constraints:
# - apollo-mapping in multiclade mode maps fastq to ref_type "species" and "strain"
# - 'softclipping' information only relevant for clade-to-central-reference aka "species"

# ------------------------------------------------------------------------------------------------ #
# base rules
# ------------------------------------------------------------------------------------------------ #

from snakemake.rules import Rule

def apply_bedtools_rule_defaults(rule_obj:Rule) -> None:
    """ extend a rule that relies on bedtools with its defaults; can't define 'conda:' in base rule """
    rule_obj.threads = config["threads"]["other"]
    rule_obj.resources = {"mem_gb": config["mem_gb"]["other"]}
    rule_obj.conda = "../envs/bedtools.yaml"
    rule_obj.container = "" #"docker://quay.io/biocontainers/bedtools:XXXXXXXXX"
    rule_obj.wildcard_constraints = {"ref_type": "species"}

rule bedgraph_to_bigwig_base:
    # basic config for bedgraph to bigwig conversion rules
    # please specify for inherited rules the directives input: output: log: (optionally params: wildcard_constraints:)
    # optionally overrule message:
    wildcard_constraints:
        ref_type="species"
    message:
        "Converting BedGraph to BigWig format"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"]
    conda:
        "../envs/ucsc_bedgraphtobigwig.yaml"
    #container:
    #    TODO: container needed here!!
    #    "docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
    shell:
        """ bedGraphToBigWig {input.bedgraph} {input.fa}.fal {output} """

# ------------------------------------------------------------------------------------------------ #
# actual rules
# ------------------------------------------------------------------------------------------------ #

rule bam_to_cov_bedgraph:
    input:
        bam=OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam",
        fa=MultiReferenceProvider.get_species_fasta_path
    output:
        OUT + "/coverage/{ref_type}/{sample}.genome.bg"
    message:
        "Generate coverage BedGraph from BAM for {wildcards.sample} on {input.fa}"
    log:
        OUT + "/log/coverage/{sample}_{ref_type}_bam_to_cov_bedgraph.log"
    shell:
        """ bedtools bamtobed -i {input.bam}  | bedtools genomecov -i stdin -g {input.fa}.fal -bg > {output} """

apply_bedtools_rule_defaults(rules.bam_to_cov_bedgraph)


rule bam_to_softclipped_bed:
    input:
        OUT + "/mapped_reads/duprem/{ref_type}/{sample}.bam"
    output:
        OUT + "/softclipped-{ref_type}/{sample}.softclipped.bed"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    conda:
        "../envs/bedtools_samtools.yaml"
    message:
        "Generate {wildcards.sample}.softclipped.bed from BAM"
    log:
        OUT + "/log/softclipped/{sample}_{ref_type}=bam_to_softclipped_bed.log"
    shell:
        """ ./workflow/scripts/convert-sam-to-softclipped-bed.py {input} > {output} """


rule bedgraph_lowcovered_regions:
    input:
        bg=OUT + "/coverage/{ref_type}/{sample}.genome.bg",
        fa=MultiReferenceProvider.get_species_fasta_path,
    output:
        OUT + "/coverage/{ref_type}/{sample}.lowcovered.bg"
    log:
        OUT + "/log/coverage/{sample}_{ref_type}__softclipped_bed_to_bedgraph.log"
    params:
        max_coverage=2,
        reference=MultiReferenceProvider.get_species_fasta_path,
    message:
        "Converting bedgraph: {wildcards.sample} on {input.fa} to lowcovered.bg"
        ###"Converting bedgraph [{wildcards.sample} on {wildcards.ref_type}/{params.reference}] to lowcovered.bg"
    shell:
        """ ( ( awk '{{ if ($4<={params.max_coverage}) {{ print $0 }} }}' {input.bg} ) && \
              ( bedtools subtract -b {input.bg} -a {input.fa}.bed | awk '{{ print $0"\t"0 }}' ) ) \
              | bedtools sort """

apply_bedtools_rule_defaults(rules.bedgraph_lowcovered_regions)


rule bed_to_bedgraph_softclipped:
    input:
        bed=OUT + "/softclipped-{ref_type}/{sample}.softclipped.bed",
        fa=MultiReferenceProvider.get_species_fasta_path,
    output:
        OUT + "/softclipped-{ref_type}/{sample}.{ref_type}.softclipped.bg"
    log:
        OUT + "/log/coverage/{sample}_{ref_type}_softclipped_bed_to_bedgraph.log"
    message:
        "Converting softclipped bed to bedgraph: {wildcards.sample} on {input.fa}"
    shell:
        """ bedtools genomecov -i {input.bed} -g {input.fa}.fal -bg > {output} """

apply_bedtools_rule_defaults(rules.bed_to_bedgraph_softclipped)


rule bedgraphs_to_access_softclipped:
    input:
        bg1=OUT + "/coverage/{ref_type}/{sample}.genome.bg",
        bg2=OUT + "/softclipped-{ref_type}/{sample}.{ref_type}.softclipped.bg",
        # factually unused
        fa=MultiReferenceProvider.get_species_fasta_path,
    output:
        OUT+ "/softclipped-{ref_type}/{sample}.{ref_type}.access-softclipped.bg"
    log:
        OUT + "/log/softclipped/{sample}_{ref_type}_access_softclipped_bedgraph.log"
    params:
        access_ratio=0.75,
    message:
        "Converting coverage & softclipped bedgraph to access-softclipped bedgraph: {wildcards.sample} on {input.fa}"
    shell:
        ## better to use bedtools unionbedg & awk fiddling;
        ## --binSize 1 == slooooooow
        ## bigwigCompare -b1 cladeV-softclipped.bw -b2 cladeV-coverage.bw --operation subtract --binSize 1
        ##               --outFileFormat bedgraph -o verschil.bg
        """ bedtools unionbedg -i {input.bg1} {input.bg2} | \
            awk '{{ if ($5>=($4*{params.access_ratio})) {{ print $0 }} }}' | cut -f 1-3,5 > {output} """

apply_bedtools_rule_defaults(rules.bedgraphs_to_access_softclipped)


use rule bedgraph_to_bigwig_base as bedgraph_to_bigwig_genome with:
    input:
        bedgraph=OUT + "/coverage/{ref_type}/{sample}.genome.bg",
        fa=MultiReferenceProvider.get_species_fasta_path,
    output:
        OUT + "/coverage/{ref_type}/{sample}.genome.bw"
    log:
        OUT + "/log/coverage/{sample}_{ref_type}_bedgrahp_to_bigwig_genome.log"

use rule bedgraph_to_bigwig_base as bedgraph_to_bigwig_softclipped with:
    input:
        bedgraph=OUT + "/softclipped-{ref_type}/{sample}.{ref_type}.softclipped.bg",
        fa=MultiReferenceProvider.get_species_fasta_path,
    output:
        OUT + "/softclipped-{ref_type}/{sample}.{ref_type}.softclipped.bw",
    log:
        OUT + "/log/coverage/{sample}_{ref_type}_bedgrahp_to_bigwig_softclipped.log"

use rule bedgraph_to_bigwig_base as bedgraph_to_bigwig_access_softclipped with:
    input:
        bedgraph=OUT + "/softclipped-{ref_type}/{sample}.{ref_type}.access-softclipped.bg",
        fa=MultiReferenceProvider.get_species_fasta_path,
    output:
        OUT + "/softclipped-{ref_type}/{sample}.{ref_type}.access-softclipped.bw"
    log:
        OUT + "/log/coverage/{sample}_{ref_type}_bedgrahp_to_bigwig_access_softclipped.log"


rule concatenate_softclipped_to_bed:
    input:
        MultiReferenceProviderConcatSoftclippedInput.get_concat_softclip_targets
    output:
        OUT + "/softclipped/multiclade.{ref_type}-softclipped.bed"
    log:
        OUT + "/log/softclipped/softclipped-concatenated-bed-{ref_type}.log"
    message:
        "concatenate all per-clade softclipped regions into united BED format"
    shell:
        # concatenate all per-clade access-softclipped regions and add {sample} as extra 5th column
        # using bedtools merge, calculate mean/min/max of each concatenated region,
        # and report in the final column which sample(s) contributed to it
        """ grep -P "\\t" {input} | sed 's|^.\+softclipped-species/\([^\.]\+\)\.species\.access-softclipped.bg:\(.\+$\)|\\2\\t\\1|' \
            | bedtools sort | bedtools merge -prec 2 -c 4,4,4,5 -o mean,min,max,distinct > {output} """

apply_bedtools_rule_defaults(rules.concatenate_softclipped_to_bed)


rule concatenate_softclipped_to_bedgraph:
    input:
        OUT + "/softclipped/multiclade.{ref_type}-softclipped.bed"
    output:
        OUT + "/softclipped/multiclade.{ref_type}-softclipped.bg"
    log:
        OUT + "/log/softclipped/softclipped-concatenated-bedgraph-{ref_type}.log"
    params:
        min_count=3,
        min_length=10
    message:
        "concatenate all per-clade softclipped regions into simplified BedGraph format"
    shell:
        # drop all softclipped regions having <=min_count and/or being <=min_length, and write to BedGraph (4 columns)
        """ cut -f 1-4 {input} \
            | awk '{{ if ($4>={params.min_count} && $3-$2>={params.min_length}) {{ print $0 }} }}' \
            > {output} """

apply_bedtools_rule_defaults(rules.concatenate_softclipped_to_bedgraph)


