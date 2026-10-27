#!/usr/bin/env python3

import argparse
import os
import warnings

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")


def read_focus_file(path):
    """
    Read one cell type per line, with comma-separated samples.

    Example:
        Osteosarcoma.cluster0: Reg, nonReg, Reg_4wpa
        Mesenchymal.cluster1: Regen_14_DPA, Uninjured
    """
    selected = []
    seen = set()

    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            if ":" not in line:
                raise ValueError(
                    f"Invalid focus file line {line_number}: "
                    "expected 'Celltype: sample1, sample2'."
                )

            celltype, sample_text = line.split(":", 1)
            celltype = celltype.strip()

            if not celltype:
                raise ValueError(
                    f"Missing cell type on focus file line {line_number}."
                )

            samples = [
                sample.strip()
                for sample in sample_text.split(",")
                if sample.strip()
            ]

            if not samples:
                raise ValueError(
                    f"No samples specified for '{celltype}' "
                    f"on focus file line {line_number}."
                )

            for sample in samples:
                key = (celltype, sample)

                if key not in seen:
                    selected.append(key)
                    seen.add(key)

    if not selected:
        raise ValueError(f"No valid selections found in {path}.")

    return selected


def read_input_list(path):
    """
    Read list of .h5ad file paths from a text file (one path per line).
    """
    h5ad_paths = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                if not os.path.exists(line):
                    raise FileNotFoundError(f"Input .h5ad file not found: {line}")
                h5ad_paths.append(line)

    if not h5ad_paths:
        raise ValueError(f"No valid .h5ad file paths found in input text file: {path}")

    return h5ad_paths


def plot_stacked_barplot(
    prop_df,
    output_path,
    title,
):
    """
    Plot a 100% stacked bar chart of subtype proportions across samples.
    """
    if prop_df.empty:
        print("Warning: No data to plot.")
        return

    fig_width = max(8, min(24, len(prop_df) * 0.6))
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    # Plot stacked bar chart
    prop_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        colormap="tab20",
        edgecolor="none",
    )

    ax.set_ylabel("Percentage", fontsize=10)
    ax.set_xlabel("Sample", fontsize=10)
    ax.set_title(title, fontsize=13, y=1.02)
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))

    # Rotate x tick labels for readability
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    # Place legend outside to the right
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
        frameon=False,
        title="Subtypes",
        title_fontsize=9,
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    print(f"Saved stacked bar plot: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute compositional proportions across single-cell "
            "celltype-sample groups and plot a stacked bar chart."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Text file containing paths to input .h5ad files (one per line).",
    )
    parser.add_argument("--prefix", required=True, help="Output prefix")
    parser.add_argument("--celltype", default="celltype", help="Column for cell type (default: celltype)")
    parser.add_argument("--cluster", default="leiden", help="Column for cluster labels (default: leiden)")
    parser.add_argument("--sample", default="sample", help="Column for sample labels (default: sample)")
    parser.add_argument("--layer", default="log1p", help="Expression layer")
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--max_cells", type=int, default=10000)
    parser.add_argument("--min_cells", type=int, default=3)
    parser.add_argument("--min_cells_per_celltype", type=int, default=10)
    parser.add_argument("--hvg", type=int, default=5000)
    parser.add_argument("--reg", type=float, default=1.0, help="Sinkhorn regularization parameter (fallback)")
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument(
        "--focus",
        default=None,
        help=(
            "Path to focus text file. Each line must be "
            "'Celltype.clusterN: sample1, sample2, ...'."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()
    np.random.seed(args.seed)

    print(f"Reading input file list from {args.input}...")
    h5ad_paths = read_input_list(args.input)

    adatas = []
    for path in h5ad_paths:
        print(f"Loading {path}...")
        ad = sc.read_h5ad(path)
        for col in [args.celltype, args.cluster, args.sample]:
            if col not in ad.obs.columns:
                raise ValueError(f"Required column '{col}' not found in {path}")
        adatas.append(ad)

    # Concatenate all loaded AnnData objects into one
    print("Concatenating AnnData objects...")
    adata = sc.concat(adatas, join="inner")

    # Dynamic composite subtype grouping: 'Osteosarcoma.cluster0'
    adata.obs["subtype_label"] = (
        adata.obs[args.celltype].astype(str)
        + ".cluster"
        + adata.obs[args.cluster].astype(str)
    )

    # Include all groups without filtering by cell count threshold
    counts = adata.obs.groupby(["subtype_label", args.sample], observed=True).size()
    valid_groups = counts.index.tolist()

    # Apply focus selection if provided
    if args.focus:
        focus_list = read_focus_file(args.focus)
        valid_groups = [g for g in valid_groups if g in focus_list]

    if not valid_groups:
        raise ValueError("No matching groups remained after applying focus selection.")

    # Filter adata to valid groups
    mask_valid = np.zeros(adata.n_obs, dtype=bool)
    for stype, smp in valid_groups:
        mask_valid |= (adata.obs["subtype_label"] == stype) & (adata.obs[args.sample] == smp)

    adata_sub = adata[mask_valid].copy()

    # Build contingency table: rows = samples, columns = subtypes, values = cell counts
    contingency = pd.crosstab(
        adata_sub.obs[args.sample],
        adata_sub.obs["subtype_label"],
    )

    # Normalize across each sample row to get proportions summing to 1.0 (100%)
    row_sums = contingency.sum(axis=1)
    row_sums[row_sums == 0] = 1.0  # Prevent division by zero
    prop_df = contingency.div(row_sums, axis=0)

    # -------------------------
    # Save Outputs
    # -------------------------
    csv_path = f"{args.prefix}_composition_matrix.csv"
    prop_df.to_csv(csv_path)
    print(f"Saved CSV matrix: {csv_path}")

    os.makedirs("figures", exist_ok=True)
    figure_path = f"figures/{args.prefix}_composition_barplot.png"

    plot_stacked_barplot(
        prop_df,
        figure_path,
        title=f"Subtype Composition Proportions ({args.prefix})",
    )
    print("Done!")


if __name__ == "__main__":
    main()
