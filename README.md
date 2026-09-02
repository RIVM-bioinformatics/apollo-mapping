<div align="center">
    <h1>Apollo Mapping</h1>
    <br />
    <h2>Reference-based mapping analysis of fungal genomes</h2>
    <br />
    <img src="https://via.placeholder.com/150" alt="pipeline logo">
</div>

## Pipeline information
* **Author(s):**            Boas van der Putten, Roxanne Wolthuis, Sofie Hofman, Ate van der Burgt
* **Organization:**         Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
* **Department:**           Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
* **Start date:**           07 - 04 - 2023
* **Refactoring episode:**  2026/Q2
* **Commissioned by:**      Thijs Bosch and Auke de Jong

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
* hatch
* Python 3.11

## Installation

### 1. clone repo (and activate the correct branch)
```bash
git clone https://github.com/RIVM-bioinformatics/apollo-mapping.git
cd apollo-mapping
git checkout v0.5.0
```

### 2a. verify if hatch is already installed system-wide; if so, skip steps 2a-b-c
```bash
hatch --version
```

### 2b. prepare installation using hatch on the RIVM's HPC (where default conda env is interfering with hatch)
```bash
# On "academische werkplek" there's interference with (mini)conda, giving this error message:
# '''error: failed to remove file `/mnt/miniconda/lib/python3.12/site-packages/idna-3.3.dist-info/INSTALLER`: Read-only file system (os error 30)'''
export PYTHONNOUSERSITE=1
export HATCH_DATA_DIR="$HOME/.hatch_isolated"
export HATCH_CACHE_DIR="$HOME/.hatch_isolated/cache"
# please verify the correct pyproject.toml project.requires-python version is used!
conda create -y -p ./local_conda python=3.10
conda activate ./local_conda
pip install hatch hatch-conda
```

### 2c. install hatch (on vanilla linux/ubuntu servers)
```bash
pip install --user hatch hatch-conda
```

### 3. create and activate hatch environment
```bash
hatch env create
hatch shell
```

In case you have still/already an hatch environment (Environment `default` already exists), remove it first
```bash
hatch env remove default
```

### Why not using the vanilla conda?

Conda installation - vanilla for most RIVM snakemake repo's - is **not** the preferred installation way anymore for apollo-mapping.
It's possible that conda installation will work, but due to increased usage of inter-repo dependencies it might be obsolete or even broken.

## Parameters & Usage

### Command for help
* 
* ``` apollo-mapping.py -h, --help``` Shows the help of the pipeline; checkout the various --help-*** submenu's
* ```-i, --input``` Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders)
* ```-o --output``` Relative or absolute path to the output directory. If none is given, an 'output' directory will be created in the current directory
          

### The base command(s) to run this program. 
```
# Get a detailed and explained overview of basic pipeline usage

python3 apollo-mapping.py -h 

# Get a detailed overview of currently supported (fungal) species

python3 apollo-mapping.py --help-species

# Get some brief examples on how to parameterize input for custom reference data

python3 apollo-mapping.py --help-customrefs

# Minimal example on how to run the pipeline (on a/the cluster):

python3 apollo-mapping.py -i [dir/to/PE/fastq_files] -o [/path/to/output/location]  

# Minimal example(s) on how to run locally:

python3 apollo_mapping.py -i [dir/to/PE/fastq_files] -o [/path/to/output/location]
    --local --no-containers --skip-kraken --species candida_auris
    --snakemake-args "cores=1" "nodes=1"

python3 apollo_mapping.py -i [dir/to/PE/fastq_files] -o [/path/to/output/location]
    --local --no-containers --skip-kraken --skip-reference-selection --species candida_auris
    --snakemake-args "cores=1" "nodes=1"

python3 apollo_mapping.py -i [dir/to/PE/fastq_files] -o [/path/to/output/location]
    --local --no-containers --skip-kraken --custom-reference /path/to/external.fa
    --snakemake-args "cores=1" "nodes=1"

- Pipeline can be run on a local machine, but expects considerable resources
- For full performance, some additional database files should me made available (e.g. kraken database)
- Make sure you've write permission in designated output dirs.
- Preferably, know a little bit on snakemake (in order to paralellize and thereby speedup your workflow)

# Example for multiclade analyses worflow (assuming you've selected the appropriate input data yourself)

rootdatadir=$HOME/RIVM
input=$rootdatadir/data-apollo-reference/test-fastq-input-cauris-clades
output=$rootdatadir/output-cauris-clades-vanilla
#input=$rootdatadir/data-apollo-reference/test-fastq-input-single-cauris
#output=$rootdatadir/output-cauris-single-sample
python3 apollo_mapping.py -i $input -o $output \
    --local --no-containers \
    --clg cauris --skip-kraken \
    --snakemake-args "cores=1" "nodes=1" --dryrun


# Example pipeline run based on testdata (on the RIVM cluster)

input=/mnt/scratch_dir/hofmansj/projects/apollo_pipelines/testdata_cauris/251121_VH01799_343_AAHKCJYM5_0004
output=$HOME/my_scratch_dir/test-apollo-output
mkdir -p $output
python3 apollo_mapping.py -i $input -o $output \
    --no-containers \
    --clg cauris --skip-kraken \
    --snakemake-args "cores=1" "nodes=1" --dryrun

# Example pipeline run based on testdata (testing --local and --custom-reference-dataset, but on the RIVM cluster)
# TODO: "isolate" a small, custom reference dataset directory for this test case

input=/mnt/scratch_dir/hofmansj/projects/apollo_pipelines/testdata_cauris/251121_VH01799_343_AAHKCJYM5_0004
output=$HOME/my_scratch_dir/test-apollo-output-local
mkdir -p $output
python3 apollo_mapping.py -i $input -o $output \
    --no-containers \
    --clg cauris --skip-kraken \
    --local --custom-reference-dataset /mnt/db/apollo/reference_new \
    --snakemake-args "cores=1" "nodes=1" --dryrun


# Example for multiclade masking (sub)worflow (assuming you've selected the appropriate input data yourself)

rootdatadir=$HOME/RIVM
input=$rootdatadir/data-apollo-reference/test-fastq-input-cauris-clades
output=$rootdatadir/output-cauris-clades
python3 apollo_mapping.py -i $input -o $output \
    --local --no-containers \
    --trigger-multiclade-masking-workflow \
    --species candida_auris \
    --snakemake-args "cores=1" "nodes=1" --dryrun
    
# after the workflow did finish, various /softclipped-species/{sampleId}.species.softclipped.bw need
# to get merged ino a final blacklist-softclipped.bed file
for bw in `find $output/softclipped-species -name "*.species.access-softclipped.bw"`
do
    echo $bw 1>&2
    ~/software/ucsc/bigWigToBedGraph $bw stdout
done | bedtools sort | bedtools merge > /tmp/softclipped-I-II-IV-V-VI.bed


```

In the (near) future, detailed information about the pipeline can be found in the [documentation](link to other docs). This documentation is mostly suitable for users that have access to the RIVM Linux environment.

## Explanation of the output
* **audit_trail:** Logs of conda, git and the pipeline, a sample sheet, the used parameters and a snakemake report.
* **clean_fastq:** cleaned fastq files.
* **identify_species:** Output of kraken and bracken for impurity detection a/o species identification.
* **log:** Log with output and error file (from the cluster) for each Snakemake rule/step that is performed.
* **mapped_reads:** Mapping output.
* **multiqc:** Multiqc output and multiqc html report.
* **qc_clean_fastq:** Quality control of clean fastq reads [juno-mapping].
* **qc_mapping:** Quality control of mapping [juno-mapping].
* **qc_raw_fastq:** Quality control of raw fastq reads [juno-mapping].
* **reference:** Reference genome used & reference species identification intermediates
* **simulated:** Fitted coverage and SNP Quality score distributions based on (depth of the) per-sample fastq 
* **variant:** Variant calling results.

## License
This pipeline is licensed with a AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       IDS-Bioinformatics
* **Email:**                ids-bioinformatics@rivm.nl  

## Acknowledgements

## TODO's
During the last upgrade of this pipeline, several (obvious) improvements were considered or became evident, that didn't make it into current stable release. See these as recommendations or a TODO list for the next wave of improvements.

- (further) harmonization & compartimentalization of apollo-mapping in relation to other RIVM CIB snakemake pipelines
  - most noticeably update juno-library to the newest concepts (a.o. argparse) used in apollo-mapping
    - see the bottom of Snakefile, using QC functions from workflow/helpers/generic_workflow_methods.py
    - and all "generic" argparse stuff now making apollo-mapping.py overly heavy
    - read the todo's in config/pipeline_params.yaml
  - refactor mostly apollo-mapping.py by moving (shared) code to other, generically reuseable repo(s)
    - rivm-ids-swc-snakeutils
    - rivm-ids-swc-name_to_be_defined
    - etc ...
  - generally, streamline apollo-mapping further with juno-mapping
    - in both directions, i.e. define submodules)
    - expand on re-using submodules
    - most noticeably upgrade "clean_fastq.smk" into "clean-fastq" submodule
    - most noticeably the future "identify-species" submodule
    - qc_mapping.smk parse_bbtools* rules should get full paths to workflow/scripts/parse_bbtools*.py
      - and next delete workflow/scripts/parse_bbtools*.py from apollo-mapping/workflow/scripts directory
  - general/harmonization: let each snakemake submodule make its own multiqc report
  - simplify rule definitions
  - generalize/automate definition of output folders and log files
- generalize (serious sample quality related) warnings into a single place, preferably the multiqc report
  - and generalize/harmonize this with other snakemake pipelines and submodules
- make management of genomic coordinate blacklists more mature
  - as soon as the 2th multi-clade species becomes a priority, fully automate the generation of blacklists (apollo-reference)
  - make rule copy_blacklist_bed species-aware (not just copy them all into the output folders of any pipeline)
- do research on the diversity of genomic properties of the reference sequences this pipeline can tackle
  - mostly work in apollo-reference
  - examples: repetitive sequences, diploids and other freaks of nature
  - named `species` with very low ANI differences (might compromise matchreference.py) 
- and have a look at minor other TODO's left in the code

## Contribution guidelines
Apollo pipelines use a [feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow). To work on features, create a branch from the `main` branch to make changes to. This branch can be merged to the main branch via a pull request. Hotfixes for bugs can be committed to the `main` branch.

Please adhere to the [conventional commits](https://www.conventionalcommits.org/) specification for commit messages. These commit messages can be picked up by [release please](https://github.com/googleapis/release-please) to create meaningful release messages.
