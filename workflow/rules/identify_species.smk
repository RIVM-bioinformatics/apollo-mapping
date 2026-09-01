import os
import shutil

# get path in the main conda environment to apollo-match-reference executable
APOLLO_MATCH_REFERENCE_BINARY = shutil.which("apollo-match-reference")
if not APOLLO_MATCH_REFERENCE_BINARY:
    raise RuntimeError("can't find apollo-match-reference in main environment")

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
        conda:
            "../envs/minimap2.yaml"
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
        conda:
            # TODO: strictly spoken, only samtools is needed in the env
            "../envs/bwa_samtools.yaml"
        container:
            "docker://staphb/bwa:0.7.17"
        log:
            OUT + "/log/identify_species/match-ref_{sample}.log",
        threads:
            1  #config["threads"]["other"]
        resources:
            mem_gb=8    #config["mem_gb"]["other"],
        shell:
            """
    {APOLLO_MATCH_REFERENCE_BINARY} {input.sam} {params.prefix} >{log} 2>&1
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
            #outdir=OUT + "/reference",
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
        conda:
            # TODO: strictly spoken, only samtools is needed in the env
            "../envs/bwa_samtools.yaml"
        #run:
        #    cmd_create_dir="mkdir -p {params.outdir}; "
        #    cmd="""{APOLLO_MATCH_REFERENCE_BINARY} {params.prefix} --species-reference {params.species_reference}"""
        #    if params.clade_reference not in (None,"None"):
        #        cmd+= " --clade-reference {params.clade_reference}"
        #    if params.exterior_fasta == "True":
        #        cmd+= " --custom"
        #    cmd+= " 2>&1>{log}"
        #    shell(cmd_create_dir+cmd)
        shell:
            """
            # base command
            CMD="{APOLLO_MATCH_REFERENCE_BINARY} {params.prefix} --species-reference {params.species_reference}"
            # optionally, expand with clade_reference
            if [ "{params.clade_reference}" != "None" ] && [ -n "{params.clade_reference}" ]; then
                CMD="$CMD --clade-reference {params.clade_reference}"
            fi
            # optionally, expand with exterior_fasta
            if [ "{params.exterior_fasta}" = "True" ]; then
                CMD="$CMD --custom"
            fi
            # execute
            eval $CMD > {log} 2>&1
            """

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
        params:
            # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
            conda_not_needed=True,
            container_not_needed=True
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
        params:
            # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
            conda_not_needed=True,
            container_not_needed=True
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
        params:
            # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
            conda_not_needed=True,
            container_not_needed=True
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
        fa = OUT + "/reference/species/{reference_species_basename, [^/]+\.(?:fa|fasta|fna)}",
        yaml = OUT + "/reference/species/{reference_species_basename, [^/]+\.(?:fa|fasta|fna)}.accessions.yaml",
    params:
        # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
        conda_not_needed=True,
        container_not_needed=True
    shell:
        # realize typically many:1 relationship of sample:reference, so file can exist already
        """
        if [ ! -f {output.fa} ]; then
            sleep 2;
            cp {input} {output.fa};
            cp {input}.accessions.yaml {output.yaml};
            # make sure files are (re)writable by snakemake
            chmod ug+w {output.fa};
            chmod ug+w {output.yaml};
        fi
        """

rule copy_blacklist_bed:
    # TODO: Once there are multiple species with a multi-strain definition,
    #       it should depend on {reference_species_basename}
    #       Please read the general to does that mention the multiclade blacklist approach
    input:
        config["apollo_reference_dir"] + "/multiclade" + "/cauris-GCA_002759435.3-blacklist.bed",
        config["apollo_reference_dir"] + "/multiclade" + "/cauris-GCA_002759435.3-vs-fastq-blacklist.bed",
        config["apollo_reference_dir"] + "/multiclade" + "/cauris-GCA_002759435.3-vs-WGA-blacklist.bed",
    output:
        bed = OUT + "/reference/species/cauris-GCA_002759435.3-blacklist.bed"
    params:
        blackistdir = config["apollo_reference_dir"] + "/multiclade",
        refoutdir = OUT + "/reference/species",
        # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
        conda_not_needed=True,
        container_not_needed=True
    shell:
        """
        # TODO: refactor to not copying ALL blacklist files in the future
        cp {params.blackistdir}/*.bed {params.refoutdir};
        # make sure files are (re)writable by snakemake
        chmod ug+w {params.refoutdir}/*.bed;
        """

rule copy_reference_strain_genomes:
    # only copy "unique" reference strain genome
    input:
        lambda wildcards: PATH_TO_REFERENCES + "/" + wildcards.reference_strain_basename
    output:
        fa = OUT + "/reference/strains/{reference_strain_basename, [^/]+\.(?:fa|fasta|fna)}",
        yaml = OUT + "/reference/strains/{reference_strain_basename, [^/]+\.(?:fa|fasta|fna)}.accessions.yaml"
    params:
        # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
        conda_not_needed=True,
        container_not_needed=True
    shell:
        # realize typically many:1 relationship of sample:reference, so file can exist already
        """
        if [ ! -f {output.fa} ]; then
            sleep 2;
            cp {input} {output.fa};
            cp {input}.accessions.yaml {output.yaml};
            # make sure files are (re)writable by snakemake
            chmod ug+w {output.fa};
            chmod ug+w {output.yaml};
        fi
        """


if MultiReferenceProvider.EXTERIOR_FASTA:
    # explicitly provided EXTERIOR_FASTA species reference
    rule copy_external_species_genome:
        input:
            config["species_reference"]
        output:
            os.path.join(OUT, "reference/species", os.path.basename(config["species_reference"]))
        params:
            # Custom metadata flags to whitelist missing conda/container for check_rule_environments validation function
            conda_not_needed=True,
            container_not_needed=True
        message:
            "copy externally provided reference genome"
        shell:
            # TODO: generate {output}.accessions.yaml file; I assume pipeline will crash if not existing
            """
            cp {input} {output};
            # make sure files are (re)writable by snakemake
            chmod ug+w {output};
            """

if MultiReferenceProvider.EXTERIOR_FASTA and config["clade_reference"] != "None":
    # explicitly provided EXTERIOR_FASTA (multi)clade reference
    rule copy_external_strain_genome:
        input:
            config["clade_reference"]
        output:
            os.path.join(OUT, "reference/strains", os.path.basename(config["clade_reference"]))
        message:
            "copy externally provided multiclade strain genome"
        shell:
            # TODO: generate {output}.accessions.yaml file; I assume pipeline will crash if not existing
            """
            cp {input} {output};
            # make sure files are (re)writable by snakemake
            chmod ug+w {output};
            """
