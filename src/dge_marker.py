#!/usr/bin/env python3

import argparse
import os
import sys

import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description="Perform DGE between marker-positive and marker-negative cells."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--marker", required=True, help="Name of the marker (e.g., EGFP, CD8A)")
    parser.add_argument("--layer", default="log1p")
    parser.add_argument("--pvalue", type=float, default=0.05)

    args = parser.parse_args()

    # Define the dynamic column name created by the extraction script
    status_col = f"{args.marker}_status"

    # Load AnnData
    print(f"Loading: {args.input}", flush=True)
    adata = ad.read_h5ad(args.input)

    # Check required annotation
    if status_col not in adata.obs.columns:
        sys.exit(f"Error: '{status_col}' not found in adata.obs. Did you run the extraction script with --marker {args.marker}?")

    # Check expression layer
    if args.layer not in adata.layers:
        sys.exit(f"Error: Layer '{args.layer}' not found in adata.layers.")

    # Check groups
    status = adata.obs[status_col].astype(str)

    positive_count = (status == "positive").sum()
    negative_count = (status == "negative").sum()

    if positive_count == 0 or negative_count == 0:
        sys.exit(
            "Error: Both positive and negative cells are required for DGE."
        )

    print(f"{args.marker}-Positive cells: {positive_count:,}")
    print(f"{args.marker}-Negative cells: {negative_count:,}")

    # Perform differential gene expression
    print("\nRunning Wilcoxon DGE...", flush=True)

    sc.tl.rank_genes_groups(
        adata,
        groupby=status_col,
        groups=["positive"],
        reference="negative",
        method="wilcoxon",
        layer=args.layer,
        use_raw=False,
        key_added="dge_marker_status"
    )

    # Extract DGE results
    dge = sc.get.rank_genes_groups_df(
        adata,
        group="positive",
        key="dge_marker_status"
    )

    # Save FULL DGE results (unfiltered)
    output_csv = f"{args.prefix}.csv"
    dge.to_csv(output_csv, index=False)
    print(f"Full DGE results saved: {output_csv}")

    # ---------------------------------------------------------
    # Generate filtered significant genes for UCell (.sig.csv)
    # ---------------------------------------------------------
    sig_mask = dge["pvals_adj"] < args.pvalue

    gene_upper = dge["names"].astype(str).str.upper()
    
    # Exclude mitochondrial, ribosomal, rat, GSM, mouse Gm, and mouse Rik genes
    exclude_mask = (
        gene_upper.str.startswith("MT-", na=False)
        | gene_upper.str.startswith("MT.", na=False)
        | gene_upper.str.startswith("RPS", na=False)
        | gene_upper.str.startswith("RPL", na=False)
        | gene_upper.str.startswith("MRPS", na=False)
        | gene_upper.str.startswith("MRPL", na=False)
        | gene_upper.str.startswith("GSM", na=False)     # Human/Rat GSM genes
        | gene_upper.str.startswith("RN", na=False)      # Rat genes (e.g., Rn18s)
        | gene_upper.str.startswith("RT1-", na=False)    # Rat MHC genes
        | gene_upper.str.startswith("LOC", na=False)     # Uncharacterized/LOC genes
        | gene_upper.str.startswith("GM", na=False)      # Mouse Gene Models (Gm12606)
        | gene_upper.str.contains("RIK", na=False)       # Mouse RIKEN clones (8030451A03Rik)
    )

    sig = dge[sig_mask & ~exclude_mask].copy()

    up = sig[sig["logfoldchanges"] > 0].sort_values(
        "logfoldchanges", ascending=False
    )
    down = sig[sig["logfoldchanges"] < 0].sort_values(
        "logfoldchanges", ascending=True
    )

    sig_ordered = np.concatenate([up.to_numpy(), down.to_numpy()], axis=0) \
        if len(up) + len(down) > 0 else np.empty((0, sig.shape[1]))

    sig_out = pd.DataFrame(sig_ordered, columns=sig.columns)

    output_sig_csv = f"{args.prefix}.sig.csv"
    sig_out.to_csv(output_sig_csv, index=False)

    print(f"Significant genes (filtered for UCell) saved: {output_sig_csv}")
    # ---------------------------------------------------------

    # Check UMAP
    if "X_umap" not in adata.obsm:
        sys.exit(
            "Error: 'X_umap' not found in adata.obsm. "
            "DGE results were saved successfully."
        )

    # Prepare UMAP
    print("\nGenerating UMAP...", flush=True)

    umap = adata.obsm["X_umap"]

    positive_mask = (status == "positive").to_numpy()
    negative_mask = (status == "negative").to_numpy()

    # Plot UMAP
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        umap[negative_mask, 0],
        umap[negative_mask, 1],
        s=2,
        label=f"{args.marker}-negative (n={negative_count:,})",
        alpha=0.7
    )

    ax.scatter(
        umap[positive_mask, 0],
        umap[positive_mask, 1],
        s=2,
        label=f"{args.marker}-positive (n={positive_count:,})",
        alpha=0.7
    )

    ax.set_title(args.prefix)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(markerscale=4, frameon=False)

    fig.tight_layout()

    # Save UMAP
    os.makedirs("figures", exist_ok=True)

    output_plot = os.path.join(
        "figures",
        f"{args.prefix}_umap.png"
    )

    fig.savefig(output_plot, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"UMAP saved: {output_plot}")

    print("\nAnalysis complete.")
    print(f"Input: {args.input}")
    print(f"Marker: {args.marker}")
    print(f"Layer: {args.layer}")
    print(f"Positive cells: {positive_count:,}")
    print(f"Negative cells: {negative_count:,}")


if __name__ == "__main__":
    main()
