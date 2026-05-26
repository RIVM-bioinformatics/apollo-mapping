""" rules that generate indices from (uncompressed) reference genome fasta files """

# A general comment on wildcard_constraints:
# - constrain input to reference folders;
# - in theory "other" fasta files might get picked up by the snakemake scheduler

rule base_generate_fasta_indices:
    # base rule directives for fasta index generation
    #log:
    #    OUT + "/log/indices/{ref_type}_{prefix}_{FUTURE_INJECT_RULE_NAME_HERE}.log"
    #message:
    #    "Indexing reference genome [{FUTURE_INJECT_RULE_NAME_HERE}]: reference/{wildcards.ref_type}/{wildcards.prefix}.fa"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],

rule bwa_index_ref:
    input:
        OUT + "/reference/{ref_type}/{prefix}.fa"
    output:
        OUT + "/reference/{ref_type}/{prefix}.fa.sa"
    wildcard_constraints:
        ref_type="species|strains"
    message:
        "Indexing reference genome [samtools faidx]: reference/{wildcards.ref_type}/{wildcards.prefix}.fa"
    #conda:
    #    "../envs/bwa_samtools.yaml"
    #container:
    #    "docker://staphb/bwa:0.7.17"
    log:
        OUT + "/log/indices/{ref_type}_{prefix}_bwa_index.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """ bwa index {input} 2> {log} """

rule samtools_faidx_ref:
    input:
        OUT + "/reference/{ref_type}/{prefix}.fa"
    output:
        OUT + "/reference/{ref_type}/{prefix}.fa.fai"
    wildcard_constraints:
        ref_type="species|strains"
    message:
        "Indexing reference genome [samtools faidx]: reference/{wildcards.ref_type}/{wildcards.prefix}.fa"
    #conda:
    #    "../envs/bwa_samtools.yaml"
    #container:
    #    "docker://staphb/bwa:0.7.17"
    log:
        OUT + "/log/indices/{ref_type}_{prefix}_samtools_faidx.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
        samtools faidx {input}
        #samtools faidx {input} 2> {log}
        """

rule fai_to_fal:
    input:
        OUT + "/reference/{ref_type}/{prefix}.fa.fai"
    output:
        OUT + "/reference/{ref_type}/{prefix}.fa.fal"
    wildcard_constraints:
        ref_type="species|strains"
    message:
        "Indexing reference genome [fai to fal]: reference/{wildcards.ref_type}/{wildcards.prefix}.fa"
    log:
         OUT + "/log/indices/{ref_type}_{prefix}_fai_to_fal.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """ cut -f 1-2 {input} 2> {log} > {output} """


rule fai_to_bed:
    input:
        OUT + "/reference/{ref_type}/{prefix}.fa.fai"
    output:
        OUT + "/reference/{ref_type}/{prefix}.fa.bed"
    wildcard_constraints:
        ref_type="species|strains"
    message:
        "Indexing reference genome [fai to bed]: reference/{wildcards.ref_type}/{wildcards.prefix}.fa"
    log:
         OUT + "/log/indices/{ref_type}_{prefix}_fai_to_bed.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """ cut -f 1-2 {input} | awk '{ OFS="\t"; print $1,0,$2 }' 2> {log} > {output} """