#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyucell as uc

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate UCell signature scores from an AnnData object, generate a unified violin plot, and a multi-panel feature plot grid with fixed scales."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the unified input .h5ad file.",
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Path to a text file containing genes of interest (one gene per line, no header).",
    )
    parser.add_argument(
        "--comparisons",
        required=True,
        help="Path to a comparisons file specifying celltype and sample pairs.",
    )
    parser.add_argument(
        "--prefix",
        default="ucell_analysis",
        help="Prefix for output files and figures (default: ucell_analysis).",
    )
    parser.add_argument(
        "--missing_genes",
        default="impute",
        choices=["impute", "skip"],
        help="How to handle genes from the signature absent in adata.var_names. "
             "'impute' (default) treats them as zero counts. 'skip' removes them.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=0.0,
        help="Fixed minimum value for UCell plot scales (default: 0.0).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=0.5,
        help="Fixed maximum value for UCell plot scales (default: 1.0).",
    )
    return parser.parse_args()

def main():
    args = parse_arguments()

    os.makedirs("figures", exist_ok=True)

    print(f"Loading AnnData from {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Robust metadata cleaning: convert to string and strip whitespace
    if "celltype" in adata.obs.columns:
        adata.obs["celltype"] = adata.obs["celltype"].astype(str).str.strip()
    if "sample" in adata.obs.columns:
        adata.obs["sample"] = adata.obs["sample"].astype(str).str.strip()

    # Diagnostic prints to verify what is actually inside your .h5ad file
    unique_celltypes = adata.obs['celltype'].unique().tolist() if "celltype" in adata.obs.columns else []
    unique_samples = adata.obs['sample'].unique().tolist() if "sample" in adata.obs.columns else []
    print(f"-> Available unique celltypes in adata.obs['celltype']: {unique_celltypes}")
    print(f"-> Available unique samples in adata.obs['sample']: {unique_samples}")

    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray()
    adata.X = np.array(adata.X, dtype=np.float32)

    if np.any(adata.X < 0):
        print(
            "WARNING: adata.X contains negative values. UCell scores should be computed "
            "on raw or size-factor-normalized (non-negative) expression data."
        )

    if "celltype" not in adata.obs.columns or "sample" not in adata.obs.columns:
        raise KeyError("AnnData .obs must contain 'celltype' and 'sample' columns.")

    print(f"Reading gene list from {args.genes}...")
    gene_list = []
    if os.path.exists(args.genes):
        with open(args.genes, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    parts = [p.strip() for p in line.replace(",", " ").split() if p.strip()]
                    gene_list.extend(parts)
    else:
        gene_list = [g.strip() for g in args.genes.split(",") if g.strip()]

    print(f"Parsed {len(gene_list)} genes from input file.")

    valid_genes = [g for g in gene_list if g in adata.var_names]
    missing_genes_set = set(gene_list) - set(valid_genes)
    if missing_genes_set:
        print(f"Warning: The following genes were not found in adata.var_names: {list(missing_genes_set)}")

    if not valid_genes:
        raise ValueError(f"None of the target genes {gene_list} were found in the AnnData var_names.")

    print(f"Matched {len(valid_genes)} valid genes present in AnnData.")

    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray()
    n_detected = (adata.X > 0).sum(axis=1)
    median_detected = int(np.median(n_detected))
    max_rank = min(median_detected, adata.n_vars)
    print(f"Median detected genes per cell: {median_detected}. Using max_rank = {max_rank}.")

    signatures = {"Custom_Signature": valid_genes}

    print("Computing global UCell scores across all cells...")
    uc.compute_ucell_scores(
        adata,
        signatures=signatures,
        chunk_size=1000,
        max_rank=max_rank,
        missing_genes=args.missing_genes,
        n_jobs=-1,
    )

    score_col = "Custom_Signature_UCell"

    global_vmin = args.vmin
    global_vmax = args.vmax
    print(f"Using fixed UCell score scale for all plots: [{global_vmin}, {global_vmax}]")

    out_csv = f"{args.prefix}_all_cells_ucell_scores.csv"
    adata.obs[["celltype", "sample", score_col]].to_csv(out_csv)
    print(f"Saved full score table to {out_csv}")

    print(f"Reading comparisons from {args.comparisons}...")
    comparisons = []
    with open(args.comparisons, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                parts = line.split()
                if len(parts) >= 2:
                    comparisons.append((parts[0].strip(), parts[1].strip()))

    print(f"Loaded {len(comparisons)} comparisons to process.")

    combined_plot_data = []
    comparison_labels_order = []
    valid_comparisons = []

    for idx, (celltype, sample) in enumerate(comparisons):
        print(f"Processing comparison {idx+1}: celltype='{celltype}', sample='{sample}'")
        
        subset_adata = adata[(adata.obs["celltype"] == celltype) & (adata.obs["sample"] == sample)].copy()
        print(f"-> Found {subset_adata.n_obs} cells matching celltype='{celltype}' and sample='{sample}'")

        if subset_adata.n_obs == 0:
            raise ValueError(
                f"ERROR: 0 cells found for celltype='{celltype}' and sample='{sample}'. "
                f"Check your comparisons file against your AnnData metadata. "
                f"Available celltypes: {unique_celltypes}, Available samples: {unique_samples}"
            )
        
        valid_comparisons.append((celltype, sample, subset_adata))

        label = f"{celltype}\n({sample})"
        obs_copy = subset_adata.obs.copy()
        obs_copy["comparison_label"] = label
        combined_plot_data.append(obs_copy)
        if label not in comparison_labels_order:
            comparison_labels_order.append(label)

    # Generate multi-panel feature plot grid
    if valid_comparisons:
        n_plots = len(valid_comparisons)
        ncols = min(3, n_plots)
        nrows = (n_plots + ncols - 1) // ncols
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4), squeeze=False)
        axes_flat = axes.flatten()

        for idx, (celltype, sample, subset_adata) in enumerate(valid_comparisons):
            ax = axes_flat[idx]
            sc.pl.umap(
                subset_adata,
                color=score_col,
                color_map="viridis",
                show=False,
                ax=ax,
                vmin=global_vmin,
                vmax=global_vmax,
            )
            ax.set_title(f"{celltype} | {sample}", fontsize=10, fontweight="bold")

        for j in range(len(valid_comparisons), len(axes_flat)):
            fig.delaxes(axes_flat[j])

        plt.suptitle(f"UCell Feature Plots ({args.prefix})", fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()

        out_umap = os.path.join("figures", f"{args.prefix}_featureplot.png")
        plt.savefig(out_umap, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved multi-panel feature plot to {out_umap}")

    # Generate unified violin plot
    if combined_plot_data:
        unified_df = pd.concat(combined_plot_data)
        unified_df["comparison_label"] = pd.Categorical(
            unified_df["comparison_label"], 
            categories=comparison_labels_order, 
            ordered=True
        )
        
        plt.figure(figsize=(max(8, len(comparison_labels_order) * 1.5), 6))
        sns.violinplot(
            data=unified_df,
            x="comparison_label",
            y=score_col,
            order=comparison_labels_order,
            palette="Set2",
            inner="box"
        )
        plt.title("Unified UCell Score Distribution Across Comparisons", fontsize=12, fontweight="bold")
        plt.xlabel("Comparison Group (Celltype | Sample)", fontsize=10)
        plt.ylabel("UCell Signature Score", fontsize=10)
        plt.xticks(rotation=45, ha="right")
        
        plt.ylim(global_vmin, global_vmax)
        plt.tight_layout()

        out_violin = os.path.join("figures", f"{args.prefix}_violin.png")
        plt.savefig(out_violin, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved unified violin plot to {out_violin}")

if __name__ == "__main__":
    main()
