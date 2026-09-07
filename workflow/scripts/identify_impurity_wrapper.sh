#!/bin/bash
#####################################################################################################################
### This is a wrapper for kraken2 and bracken used in the snakefile.                                                           
###
### Usage: {scripts_path}/identify_impurity_wrapper.sh {wildcards.sample} {input.r1} {input.r2}                   ###
###        {params.kraken_db_dir} {params.kraken2_kreport} {params.bracken_s} {threads} \                         ###
###        {output.kraken_trigger_txt} {log} > {log} 2>&1                                                         ###
### Reason: Kraken2 and bracken are skipped if the flag --skip-kraken is activated. This functionality was        ###
###         implemented to allow users to bypass kraken due to the large database needed (100 GB). I had to       ###
###         put this simple bash code into a script because of the Snakemake-bash strict settings, see:           ###
### https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-with-errors-about-an-unbound-variable-what-s-wrong
#####################################################################################################################

set -euo pipefail

WILDCARDS_SAMPLE="$1"
R1="$2"
R2="$3"
KRAKEN_DB_DIR="$4"
KRAKEN2_KREPORT="$5"
BRACKEN_S="$6"
THREADS="$7"
OUTPUT_KRAKEN_TRIGGER_TXT="$8"
LOG="$9"

if [[ -z "$KRAKEN_DB_DIR" || "$KRAKEN_DB_DIR" == "None" ]]; then
    # No Kraken database available; skip this rule completely
    # write the output.kraken_trigger_txt file -> SKIPPED
    echo "kraken2 and bracken impurity analysis of ${WILDCARDS_SAMPLE}: SKIPPED" > "$OUTPUT_KRAKEN_TRIGGER_TXT"
else
    # Kraken database available!
    # Adding --confidence 0.05 causes 0 kmers assigned to species in some samples
    # That breaks bracken
    kraken2 --db "$KRAKEN_DB_DIR" \
        --threads "$THREADS" \
        --report "$KRAKEN2_KREPORT" \
        --gzip-compressed \
        --paired \
        --output - \
        "$R1" "$R2" &> "$LOG"
    bracken -d "$KRAKEN_DB_DIR" \
        -i "$KRAKEN2_KREPORT" \
        -o "$BRACKEN_S" \
        -r 150 \
        -l S \
        -t 0 &>> "$LOG"
    # write the output.kraken_trigger_txt file -> DONE
    echo "kraken2 and bracken impurity analysis of ${WILDCARDS_SAMPLE} [on ${KRAKEN_DB_DIR}]: DONE" > "$OUTPUT_KRAKEN_TRIGGER_TXT"
fi