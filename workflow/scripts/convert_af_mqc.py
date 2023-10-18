#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import shutil
import argparse


def main(args):
    shutil.copy(args.mqc, args.output)
    df = pd.read_csv(args.input, sep="\t")
    # Remove multi-allelic counts
    df_monoallelic = df[~df["MULTI-ALLELIC"]]
    # Remove multi-allelic columns
    df_clean = df_monoallelic[["Sample", "AF", "Count"]]
    # Pivot to wide format and clean column names
    df_clean_pivoted = pd.pivot_table(
        df_clean, index="Sample", columns="AF", values=["Count"], fill_value=0
    )
    df_clean_pivoted.columns = [
        "_".join(col).strip() for col in df_clean_pivoted.columns.values
    ]
    df_clean_pivoted.reset_index(inplace=True)
    df_clean_pivoted.columns = df_clean_pivoted.columns.map(
        {"Sample": "Sample", "Count_0.500": "Heterozygous", "Count_1.00": "Homozygous"}
    )
    # Calculate ratio heterozygous (mono-allelic) variants from total number of (mono-allelic) variants
    df_clean_pivoted["Heterozygosity_ratio"] = df_clean_pivoted["Heterozygous"] / (
        df_clean_pivoted["Heterozygous"] + df_clean_pivoted["Homozygous"]
    )
    df_clean_pivoted.to_csv(args.output, mode="a", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--mqc", help="Append data to this MultiQC header", type=Path, required=True
    )
    args = parser.parse_args()

    main(args)
