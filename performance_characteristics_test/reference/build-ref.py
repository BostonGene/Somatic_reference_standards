#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import scipy.stats
from tqdm import tqdm


# ---------------------------
# Functions
# ---------------------------

def add_vaf_coef_columns(df, column_list, fractions):
    """
    Compute Pearson correlation and linear regression slope
    between observed VAFs and expected fractions.
    """

    df_t = df[column_list].T

    vaf_r_coef = []
    vaf_l_coef = []

    for col in tqdm(df_t.columns, desc="Computing VAF coefficients"):
        try:
            x = np.array(df_t[col], dtype=float)
            y = np.array(fractions, dtype=float)

            m, _ = np.polyfit(x, y, 1)
            r, _ = scipy.stats.pearsonr(x, y)

        except Exception:
            m = 0
            r = 0

        vaf_r_coef.append(r)
        vaf_l_coef.append(m)

    df["VAF_r_coefficient"] = vaf_r_coef
    df["VAF_l_coefficient"] = vaf_l_coef

    return df


def read_and_prepare_maf(path, sample_name):
    """
    Read MAF and rename Tumor_VAF column.
    """

    cols = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Tumor_VAF",
    ]

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=cols,
    )

    df = df.rename(columns={"Tumor_VAF": f"Tumor_VAF_{sample_name}"})

    return df


# ---------------------------
# Main
# ---------------------------

def main(args):

    # Expected fractions
    fractions = args.fractions

    # Initialize empty reference table
    ref_df = pd.DataFrame(
        columns=[
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        ]
    )

    # Merge all MAFs
    for sample_name, maf_path in zip(args.samples, args.mafs):
        maf_df = read_and_prepare_maf(maf_path, sample_name)

        ref_df = ref_df.merge(
            maf_df,
            on=[
                "Hugo_Symbol",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            how="outer",
        )

    # Check VAF columns
    vaf_columns = args.vaf_columns

    missing = set(vaf_columns) - set(ref_df.columns)
    if missing:
        raise ValueError(f"Missing VAF columns: {missing}")

    # Compute coefficients
    ref_df = add_vaf_coef_columns(ref_df, vaf_columns, fractions)

    # True somatic variant logic
    ref_df["True_vars"] = np.where(
        (
            (ref_df["VAF_r_coefficient"] >= args.r_threshold)
            & (ref_df["VAF_l_coefficient"] > args.slope_min)
        )
        | (
            (ref_df["VAF_r_coefficient"] >= args.r_threshold)
            & (ref_df["VAF_l_coefficient"] < args.slope_max)
        ),
        "PASS",
        "",
    )

    # Final compact output
    out_df = ref_df[
        [
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        ]
    ].copy()

    out_df["IS_non_FP"] = np.where(ref_df["True_vars"] == "PASS", "YES", "NO")

    out_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Build somatic reference variants from dilution series"
    )

    parser.add_argument(
        "--mafs",
        nargs="+",
        required=True,
        help="List of input MAF files",
    )

    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample names matching MAF files order",
    )

    parser.add_argument(
        "--vaf-columns",
        nargs="+",
        required=True,
        help="VAF column names to use for regression",
    )

    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        required=True,
        help="Expected dilution fractions (e.g. 0.1 0.2 0.3 0.5 0.75 1.0)",
    )

    parser.add_argument(
        "--r-threshold",
        type=float,
        default=0.8,
        help="Minimum Pearson r",
    )

    parser.add_argument(
        "--slope-min",
        type=float,
        default=0.7,
        help="Minimum regression slope",
    )

    parser.add_argument(
        "--slope-max",
        type=float,
        default=1.3,
        help="Maximum regression slope",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output reference TSV",
    )

    args = parser.parse_args()
    main(args)
