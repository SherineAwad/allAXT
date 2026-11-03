#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Find genes that are statistically significant and changing in the Base comparison "
            "(Osteo-nonReg vs. Uninjured) and flat/non-significant in the Group comparison "
            "(Osteo-Reg vs. Uninjured)."
        )
    )

    parser.add_argument(
        "--base",
        required=True,
        help="Base CSV file (Osteo-nonReg vs. Uninjured; expects significant changes)."
    )
    parser.add_argument(
        "--group",
        required=True,
        help="Group CSV file (Osteo-Reg vs. Uninjured; expects flat/non-significant behavior)."
    )
    parser.add_argument(
        "--output",
        default="rescued_genes.csv",
        help="Output CSV path for the final filtered gene list. Default: rescued_genes.csv."
    )
    parser.add_argument(
        "--gene-col",
        default="names",
        help="Column name containing gene symbols or IDs. Default: 'names'."
    )

    # Base Filters (The "Changing" Criteria)
    parser.add_argument(
        "--base-adj-pvalue-cutoff",
        type=float,
        default=0.05,
        help="Base max adjusted p-value threshold (pvals_adj < cutoff). Default: 0.05."
    )
    parser.add_argument(
        "--base-min-logfc",
        type=float,
        default=0.25,
        help="Base minimum absolute log-fold change (abs(logfoldchanges) >= cutoff). Default: 0.25."
    )

    # Group Filters (The "Flat/Unchanged" Criteria)
    parser.add_argument(
        "--group-adj-pvalue-cutoff",
        type=float,
        default=0.05,
        help="Group minimum adjusted p-value threshold (pvals_adj > cutoff). Default: 0.05."
    )
    parser.add_argument(
        "--group-max-logfc",
        type=float,
        default=0.1,
        help="Group maximum absolute log-fold change (abs(logfoldchanges) <= cutoff). Default: 0.1."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    path_base = Path(args.base)
    path_group = Path(args.group)
    output_path = Path(args.output)

    try:
        df_base = pd.read_csv(path_base)
        df_group = pd.read_csv(path_group)
    except Exception as exc:
        sys.exit(f"Error reading input files: {exc}")

    # Validate columns exist in both dataframes
    for name, df in [("Base", df_base), ("Group", df_group)]:
        for col in [args.gene_col, "logfoldchanges", "pvals_adj"]:
            if col not in df.columns:
                sys.exit(f"Error: {name} file is missing required column '{col}'.")

    # Coerce numeric columns to handle potential string formatting issues
    for df in [df_base, df_group]:
        df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
        df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")

    # 1. Filter Base: Keep significant genes with a minimum effect size
    mask_base = (
        (df_base["pvals_adj"] < args.base_adj_pvalue_cutoff) &
        (df_base["logfoldchanges"].abs() >= args.base_min_logfc)
    )
    filtered_base = df_base.loc[mask_base].copy()

    # 2. Merge with Group on the gene identifier column to align common genes
    merged = pd.merge(
        filtered_base,
        df_group,
        on=args.gene_col,
        suffixes=("_base", "_group")
    )

    # 3. Filter Group side: Keep only genes that are flat / non-significant
    mask_group = (
        (merged["pvals_adj_group"] > args.group_adj_pvalue_cutoff) &
        (merged["logfoldchanges_group"].abs() <= args.group_max_logfc)
    )
    final_result = merged.loc[mask_group].copy()

    try:
        final_result.to_csv(output_path, index=False)
    except Exception as exc:
        sys.exit(f"Error writing {output_path}: {exc}")

    print(f"Base total rows: {len(df_base)}")
    print(f"Base significant subset: {len(filtered_base)}")
    print(f"Group total rows: {len(df_group)}")
    print(f"Final retained genes (Rescued genes): {len(final_result)}")
    print(f"Output saved to: {output_path}")


if __name__ == "__main__":
    main()
