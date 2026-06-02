import os

if config.get("skip_reference_selection",None) == "False":
    # vanilla: minimap2 --> match_ref --> assign_reference
    rule minimap2:
        input:
            # Realize match-ref is not critically depending on cleaned fastq data.
            # Avdb prefers to do match-ref on "raw" fastq data (although that might be a personal choice ...)
            # In case this is decided on differently, change to clean_fastq/ data as input
            #r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
            #r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
            r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
            index = config["identify_species_index"]
        output:
            sam=OUT + "/reference/{sample}__minimap2.sam"
        message:
            "map reads to identify-species-index {input.index} for sample {wildcards.sample}"
        #conda:
        #    "../envs/bwa_samtools.yaml"
        log:
            OUT + "/log/identify_species/minimap2_{sample}.log",
        threads: 10 #config["threads"]["other"]
        resources:
            mem_gb= 4 #config["mem_gb"]["other"],
        shell:
            """
    minimap2 -ax sr {input.index} {input.r1} {input.r2} > {output.sam}
            """

    rule match_ref:
        input:
            sam=OUT + "/reference/{sample}__minimap2.sam"
        output:
            out1=OUT + "/reference/{sample}-match-ref-reference.tsv",
            out2=OUT + "/reference/{sample}-match-ref-taxid.tsv",
            out3=OUT + "/reference/{sample}-samtools-stats.tsv",
            # !important! this attribute is used for snakemake scheduling
            reference=OUT + "/reference/{sample}-references.yml"
        params:
            prefix=OUT + "/reference/{sample}",
        message:
            "identify species using apollo-match-reference for sample {wildcards.sample}"
        #conda:
        #    "../envs/bwa_samtools.yaml"
        #container:
        #    "docker://staphb/bwa:0.7.17"
        log:
            OUT + "/log/identify_species/match-ref_{sample}.log",
        threads:
            1  #config["threads"]["other"]
        resources:
            mem_gb=8    #config["mem_gb"]["other"],
        shell:
            """
    apollo-match-reference {input.sam} {params.prefix} 2>&1>{log}
            """

    # !important!
    CONDITIONAL_TARGETS += rules.minimap2.output
    CONDITIONAL_TARGETS += rules.match_ref.output


if config.get("species_reference",None) not in (None,"None"):
    # regardless if match_ref performed: forced_ref --> assign_reference
    rule forced_ref:
        input:
            r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        params:
            outdir=OUT + "/reference",
            prefix=OUT + "/reference/{sample}",
            # overtake those user-specified to CLI entry point apollo_mapping.py
            species_reference=config["species_reference"],
            clade_reference=config["clade_reference"],
            exterior_fasta=config["exterior_fasta"],
        message:
            "set enforced reference species using apollo-reference/match-ref.py for sample {wildcards.sample}"
        output:
            # !important! this attribute is used for snakemake scheduling
            reference=OUT + "/reference/{sample}-forced_references.yml"
        log:
            OUT + "/log/identify_species/forced-ref_{sample}.log",
        run:
            cmd_create_dir="mkdir -p {params.outdir}; "
            cmd="""apollo-match-reference {params.prefix} --species-reference {params.species_reference}"""
            if params.clade_reference not in (None,"None"):
                cmd+= " --clade-reference {params.clade_reference}"
            if params.exterior_fasta == "True":
                cmd+= " --custom"
            cmd+= " 2>&1>{log}"
            shell(cmd_create_dir+cmd)

    # !important!
    CONDITIONAL_TARGETS += rules.forced_ref.output

if config.get("species_reference",None) not in (None,"None") and config.get("skip_reference_selection",None) == "False":
    # --species <name>, but match_ref has been executed ( no --skip-reference-selection )
    checkpoint assign_reference:
        input:
            overruled=OUT + "/reference/{sample}-references.yml",
            reference=OUT + "/reference/{sample}-forced_references.yml"
        output:
            reference=OUT + "/reference/{sample}-assigned_references.yml"
        message:
            "assign reference genome for sample {wildcards.sample} [scenario: forced]"
        shell:
            """ cp {input.reference} {output.reference} """

elif config.get("species_reference",None) not in (None,"None"):
    # ( --species <name> or --reference_fasta <fasta> ) AND --skip-reference-selection
    checkpoint assign_reference:
        input:
            reference=OUT + "/reference/{sample}-forced_references.yml"
        output:
            reference=OUT + "/reference/{sample}-assigned_references.yml"
        message:
            "assign reference genome for sample {wildcards.sample} [scenario: skipped]"
        shell:
            """ cp {input.reference} {output.reference} """

else:
    # vanilla: switch to the reference as indicated by rule match_ref
    checkpoint assign_reference:
        input:
            reference=OUT + "/reference/{sample}-references.yml"
        output:
            reference=OUT + "/reference/{sample}-assigned_references.yml"
        message:
            "assign reference genome for sample {wildcards.sample} [scenario: vanilla]"
        shell:
            """ cp {input.reference} {output.reference} """

# !important!
CONDITIONAL_TARGETS += rules.assign_reference.output

################################################################################################################

rule copy_reference_species_genomes:
    # only copy "unique" reference species genome
    input:
        lambda wildcards: PATH_TO_REFERENCES + "/" + wildcards.reference_species_basename
    output:
        OUT + "/reference/species/{reference_species_basename, [^/]+\.(?:fa|fasta|fna)}"
    shell:
        """if [ ! -f {output} ]; then sleep 2; cp {input} {output}; fi"""

rule copy_reference_strain_genomes:
    # only copy "unique" reference strain genome
    input:
        lambda wildcards: PATH_TO_REFERENCES + "/" + wildcards.reference_strain_basename
    output:
        OUT + "/reference/strains/{reference_strain_basename, [^/]+\.(?:fa|fasta|fna)}"
    shell:
        """if [ ! -f {output} ]; then sleep 2; cp {input} {output}; fi"""


if MultiReferenceProvider.EXTERIOR_FASTA:

    rule copy_external_species_genome:
        input:
            config["species_reference"]
        output:
            os.path.join(OUT, "reference/species", os.path.basename(config["species_reference"]))
        message:
            "copy externally provided reference genome"
        shell:
            """cp {input} {output}"""

if MultiReferenceProvider.EXTERIOR_FASTA and config["clade_reference"] != "None":

    rule copy_external_strain_genome:
        input:
            config["clade_reference"]
        input:
            os.path.join(OUT, "reference/strains", os.path.basename(config["clade_reference"]))
        message:
            "copy externally provided multiclade strain genome"
        shell:
            """cp {input} {output}"""

################################################################################################################

if False:

    ruleorder: use_forced_reference_yaml > use_matched_ref_yaml

    rule use_forced_reference_yaml:
        input: OUT + "/reference/{sample}-FORCED-REFERENCES.yaml"
        output: reference = OUT + "/reference/{sample}-references.yaml"
        shell: "cp {input} {output}"

    rule use_matched_ref_yaml:
        input: OUT + "reference/{sample}-match-ref-references.yaml"
        output: reference = OUT + "/reference/{sample}-references.yaml"
        shell: "cp {input} {output}"

    #checkpoint assigned_reference:
    #    input: []
    #    shell: []
    #    output: reference = OUT + "/reference/{sample}-references.yaml"




################################################################################################################


if False:

    def get_reference_species_genome(wildcards):
        """ Force snakemake to wait untill checkpoint is done in order to get assigned reference

        """
        target_yaml = checkpoints.match_ref.get(**wildcards).output.reference

        # can only return results when checkpoint file exists
        if os.path.exists(target_yaml):
            with open(target_yaml) as f:
                data = yaml.safe_load(f)
                # !important! will be added to wildcards
                reference_species_basename = data["species"]["fasta"]
                # !important! path should exists
                return os.path.join( OUT, "reference", "species", data["species"]["fasta"] )

        # Snakemake will call function again once checkpoint has been passed.
        return "WAITING_FOR_CHECKPOINT ... "

    def get_optional_reference_strain_genome(wildcards):
        """ Force snakemake to wait untill checkpoint is done in order to get assigned reference

        """
        target_yaml = checkpoints.match_ref.get(**wildcards).output.reference

        # can only return results when checkpoint file exists
        if os.path.exists(target_yaml):
            with open(target_yaml) as f:
                data = yaml.safe_load(f)
                if data.get("strain"):
                    # !important! will be added to wildcards
                    reference_strain_basename = data["strain"]["fasta"]
                    # !important! path should exists
                    return os.path.join( OUT, "reference", "strains", data["strain"]["fasta"] )
                else:
                    # !important! Not elegiable
                    return []

        # Snakemake will call function again once checkpoint has been passed.
        return "WAITING_FOR_CHECKPOINT ... "

if False:

    rule that_has_reference_genome_dynamically:
        input:
            r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
            reference=MultiReferenceProvider.get_species_fasta_path
            #reference=get_reference_species_genome
        output:
            #sam=OUT + "/reference/{sample}__minimap2.sam"
            txt=OUT + "/reference/{sample}__just_checking_if_it_works.txt"
        message:
            "let's say this is the rule alike map_clean_reads.smk --> bwa_mem"
        threads: 10 #config["threads"]["other"]
        resources:
            mem_gb= 4 #config["mem_gb"]["other"],
        shell:
            """
            echo "{input.r1} {input.r2} {input.reference}" > {output.txt}
            """

    rule that_has_reference_strain_dynamically:
        input:
            r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
            reference=MultiReferenceProvider.get_strain_fasta_path
            #reference=get_optional_reference_strain_genome
        output:
            #sam=OUT + "/reference/{sample}__minimap2.sam"
            txt=OUT + "/reference/{sample}__strain__just_checking_if_it_works.txt"
        message:
            "let's say this is the rule alike map_clean_reads.smk --> bwa_mem"
        threads: 10 #config["threads"]["other"]
        resources:
            mem_gb= 4 #config["mem_gb"]["other"],
        shell:
            """
            echo "{input.r1} {input.r2} {input.reference}" > {output.txt}
            """

    rule mapping_pipeline_start:
        input:
            r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
            r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
            reference=MultiReferenceProvider.get_ref_path
        output:
            txt=OUT + "/mapped_reads/{sample}/{ref_type}/done.txt"
        shell:
            "echo 'Mapping {wildcards.sample} {input.r1} {input.r2} ({wildcards.ref_type}) using {input.reference}' > {output.txt}"


if False:

        rule bwa_index_ref:
            input:
                "{prefix}.fa"
            output:
                "{prefix}.fa.sa"
            message:
                "Indexing reference genome: {wildcards.prefix}"
            #conda:
            #    "../envs/bwa_samtools.yaml"
            #container:
            #    "docker://staphb/bwa:0.7.17"
            log:
                "{prefix}_bwa_index.log"
            threads: config["threads"]["other"]
            resources:
                mem_gb=config["mem_gb"]["other"],
            shell:
                """
                #bwa index {input} 2> {log}
                echo "YAHOOOOOOOOOOOO {input}" > {output}
                """

        rule samtools_faidx_ref:
            input:
                OUT + "/reference/{folder}/{prefix}.fa"
            output:
                OUT + "/reference/{folder}/{prefix}.fa.fai"
            wildcard_constraints:
                # constrain input to reference folders!
                folder="species|strains"
            message:
                "Indexing reference genome [samtools faidx]: reference/{wildcards.folder}/{wildcards.prefix}.fa"
            #conda:
            #    "../envs/bwa_samtools.yaml"
            #container:
            #    "docker://staphb/bwa:0.7.17"
            #log:
            #    "{prefix}_bwa_index.log"
            #threads: config["threads"]["other"]
            #resources:
            #    mem_gb=config["mem_gb"]["other"],
            shell:
                """
                samtools faidx {input}
                #samtools faidx {input} 2> {log}
                """

if False:

        rule mapping_pipeline_start:
            input:
                r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
                r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
                reference=MultiReferenceProvider.get_ref_path,
                # "fake" input needed for DAG construction (and bwa index files should exist)
                idx = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".sa",
                _idx2 = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".fai",
                _idx3 = lambda wildcards: MultiReferenceProvider.get_ref_path(wildcards) + ".fal"
            output:
                txt=OUT + "/mapped_reads/{sample}/{ref_type}/done.txt"
            shell:
                "echo 'Mapping {wildcards.sample} {input.r1} {input.r2} ({wildcards.ref_type}) using {input.reference}' > {output.txt}"


        rule fake_sam_to_sorted_bam:
            input:
                #sam=OUT + "/mapped_reads/raw/{sample}.sam",
                sam=OUT + "/mapped_reads/{sample}/{ref_type}/done.txt",
            output:
                #bam=OUT + "/mapped_reads/sorted/{sample}.bam",
                bam=OUT + "/mapped_reads/{sample}/{ref_type}/done.txt.sorted.bam",
            message:
                "Convert sam to sorted bam for {wildcards.sample}"
            log:
                OUT + "/log/sam_to_sorted_bam/{sample}_{ref_type}.log",
            shell:
                """
                #samtools view -b -@ {threads} {input.sam} 2>{log} | \
                #samtools sort -@ {threads} - 1> {output.bam} 2>>{log}
                cat {input.sam} > {output.bam}
                """
