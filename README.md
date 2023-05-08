<div align="center">
    <h1>Apollo Mapping</h1>
    <br />
    <h2>Reference-based mapping analysis of fungal genomes</h2>
    <br />
    <img src="https://via.placeholder.com/150" alt="pipeline logo">
</div>

## Pipeline information
* **Author(s):**            Boas van der Putten, Roxanne Wolthuis
* **Organization:**         Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
* **Department:**           Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
* **Start date:**           07 - 04 - 2023
* **Commissioned by:**      Thijs Bosch

## About this project
Apollo-mapping is the first pipeline created in the Apollo pipeline series. The Goal of these pipelines is to set up a routine surveillance for fungi (A.fumigatus, Candida). The apollo-mapping pipeline is created with the juno-template and juno-library.

The input of the pipeline is raw Illumina paired-end data  in the form of two fastq files (with extension .fastq, .fastq.gz, .fq or .fq.gz), containing the forward and the reversed reads ('R1' and 'R2' must be part of the file name, respectively).

The pipeline uses the following tools(NOT COMPLETE):
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (Andrews, 2010) is used to assess the quality of the raw Illumina reads
2. [FastP](https://github.com/OpenGene/fastp) (Chen, Zhou, Chen and Gu, 2018) is used to remove poor quality data and adapter sequences 
3. [Picard](https://broadinstitute.github.io/picard/) determines the library fragment lengths
4. [MultiQC](https://multiqc.info/) (Ewels, Magnusson, Lundin, & Käller, 2016) is used to summarize analysis results and quality assessments in a single report for dynamic visualization.
5. [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Bracken](http://ccb.jhu.edu/software/bracken/) for identification of fungal species.  

## Prerequisities
* Linux environment
* (mini)conda
* Python 3.11

## Installation
1. Clone the repository.
```
git clone https://github.com/RIVM-bioinformatics/apollo-mapping.git
```

2. Go to the pipeline directory.
```
cd apollo-mapping
```

3. Create & activate mamba environment.
```
conda env update -f envs/mamba.yaml
```
```
conda activate mamba
```

4. Create & activate apollo environment.
```
mamba env update -f envs/apollo_mapping.yaml
```
```
conda activate apollo_mapping
```

5. Example of run:
```
python3 apollo_mapping.py -i [input] -o [output] -s [species]
```

## Parameters & Usage
### Command for help
* ```-h, --help``` Shows the help of the pipeline

### Required parameters
* ```-i, --input``` Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders)
* ```-s, --species``` Species to use, choose from: ['candida_auris', 'aspergillus_fumigatus']

### Optional parameters
* ```-o --output``` Relative or absolute path to the output directory. If none is given, an 'output' directory will be created in the current directory
* ```-w, --workdir``` Relative or absolute path to the working directory. If none is given, the current directory is used.
* ```-ex, --exclusionfile``` Path to the file that contains samplenames to be excluded.
* ```-p, --prefix``` Conda or singularity prefix. Basically a path to the place where you want to store the conda environments or the singularity images.
* ```-l, --local``` If this flag is present, the pipeline will be run locally (not attempting to send the jobs to an HPC cluster**). The default is to assume that you are working on a cluster. **Note that currently only LSF clusters are supported.
* ```-tl, --time-limit``` Time limit per job in minutes (passed as -W argument to bsub). Jobs will be killed if not finished in this time.
* ```-u, --unlock``` Unlock output directory (passed to snakemake).
* ```-n, --dryrun``` Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake).
* ```-q, --queue``` Name of the queue that the job will be submitted to if working on a cluster.
* ```-mpt, --mean-quality-treshold``` Phred score to be used as threshold for cleaning (filtering) fastq files.
* ```-ws, --window-size``` Window size to use for cleaning (filtering) fastq files.
* ```-ml, --minimum-lenth``` Minimum length for fastq reads to be kept after trimming.
* ```--no-containers``` Use conda environments instead of containers.
* ```--snakemake-args``` Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html).
* ```--reference``` Reference genome to use default is chosen based on species argument, defaults per species can be found in: /mnt/db/apollo/mapping/[species]
* ```--db-dir``` Kraken2 database directory (should include fungi!)                 

### The base command to run this program. 
```
python3 apollo-mapping.py -i [dir/to/fasta_or_fastq_files] -s [species] 
```

### An example on how to run the pipeline.
```
python3 apollo-mapping.py -i [dir/to/fasta_or_fastq_files] -o [/path/to/output/location] -s aspergillus_fumigatus 
```

Detailed information about the pipeline can be found in the [documentation](link to other docs). This documentation is only suitable for users that have access to the RIVM Linux environment.

## Explanation of the output
* **audit_trail:** Logs of conda, git and the pipeline, a sample sheet, the used parameters and a snakemake report.
* **clean_fastq:** cleaned fastq files.
* **identify_species:** Output of kraken and bracken for species identification.
* **log:** Log with output and error file from the cluster for each Snakemake rule/step that is performed.
* **mapped_reads:** Mapping output.
* **multiqc:** Multiqc output and multiqc html report.
* **qc_clean_fastq:** Quality control of clean fastq reads.
* **qc_mapping:** Quality control of mapping.
* **reference:** Reference genome used.
* **variant:** Variant calling results.

## Issues
* This pipeline only works on the RIVM cluster.

## Future ideas for this pipeline
* Make this pipeline available and user friendly for users outside RIVM.

## License
This pipeline is licensed with a AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       IDS-Bioinformatics
* **Email:**                ids-bioinformatics@rivm.nl  

## Acknowledgements


## Contribution guidelines
Apollo pipelines use a [feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow). To work on features, create a branch from the `main` branch to make changes to. This branch can be merged to the main branch via a pull request. Hotfixes for bugs can be committed to the `main` branch.

Please adhere to the [conventional commits](https://www.conventionalcommits.org/) specification for commit messages. These commit messages can be picked up by [release please](https://github.com/googleapis/release-please) to create meaningful release messages.
