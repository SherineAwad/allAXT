#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot gene expression per sample with a shared "
            "colour scale across all samples for each gene."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input AnnData h5ad file."
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output figure filenames."
    )

    parser.add_argument(
        "--genes",
        required=True,
        help="Text file containing one gene per row."
    )

    return parser.parse_args()


def safe_filename(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def main():
    args = parse_args()

    output_dir = Path("figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading AnnData: {args.input}")
    adata = sc.read_h5ad(args.input)

    if "sample" not in adata.obs.columns:
        raise KeyError("'sample' column not found in adata.obs.")

    if "X_umap" not in adata.obsm:
        raise KeyError("'X_umap' not found in adata.obsm.")

    with open(args.genes, "r") as handle:
        genes = list(dict.fromkeys(
            line.strip()
            for line in handle
            if line.strip()
        ))

    if not genes:
        raise ValueError(f"No genes found in {args.genes}.")

    valid_genes = [
        gene for gene in genes
        if gene in adata.var_names
    ]

    missing_genes = [
        gene for gene in genes
        if gene not in adata.var_names
    ]

    if missing_genes:
        print(
            "WARNING: Genes not found and skipped: "
            + ", ".join(missing_genes)
        )

    if not valid_genes:
        raise ValueError(
            "None of the requested genes were found in adata.var_names."
        )

    samples = adata.obs["sample"].dropna().unique().tolist()

    if not samples:
        raise ValueError("No samples found in adata.obs['sample'].")

    embedding = np.asarray(adata.obsm["X_umap"])

    if embedding.shape[0] != adata.n_obs or embedding.shape[1] < 2:
        raise ValueError(
            "X_umap must contain at least two coordinates for every cell."
        )

    x_coord = embedding[:, 0]
    y_coord = embedding[:, 1]

    print(f"Cells: {adata.n_obs:,}")
    print(f"Samples: {len(samples)}")
    print(f"Genes: {len(valid_genes)}")

    for gene in valid_genes:
        print(f"Plotting {gene}...")

        expression = adata[:, gene].X

        if hasattr(expression, "toarray"):
            expression = expression.toarray().ravel()
        else:
            expression = np.asarray(expression).ravel()

        if not np.all(np.isfinite(expression)):
            raise ValueError(
                f"Non-finite expression values found for {gene}."
            )

        # Shared colour scale across all samples for this gene.
        vmin = 0
        vmax = max(float(expression.max()), 1e-9)

        n_samples = len(samples)
        ncols = min(3, n_samples)
        nrows = int(np.ceil(n_samples / ncols))

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(5 * ncols + 0.8, 5 * nrows),
            squeeze=False,
            sharex=True,
            sharey=True
        )

        axes_flat = axes.ravel()
        scatter = None

        for ax, sample in zip(axes_flat, samples):
            mask = (
                adata.obs["sample"]
                .eq(sample)
                .fillna(False)
                .to_numpy()
            )

            indices = np.flatnonzero(mask)
            order = indices[
                np.argsort(expression[indices], kind="stable")
            ]

            scatter = ax.scatter(
                x_coord[order],
                y_coord[order],
                c=expression[order],
                cmap="Reds",
                vmin=vmin,
                vmax=vmax,
                s=8,
                linewidths=0,
                rasterized=True
            )

            ax.set_title(f"{sample} (n={mask.sum():,})")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            ax.set_aspect("equal", adjustable="box")

        for ax in axes_flat[n_samples:]:
            ax.set_visible(False)

        # Reserve a dedicated area on the right for the colour bar.
        fig.subplots_adjust(
            left=0.08,
            right=0.86,
            bottom=0.08,
            top=0.90,
            wspace=0.25,
            hspace=0.30
        )

        # Place the colour bar outside the sample panels.
        cbar_ax = fig.add_axes([0.89, 0.20, 0.018, 0.60])

        fig.colorbar(
            scatter,
            cax=cbar_ax,
            label="Expression"
        )

        fig.suptitle(
            f"{gene} expression across samples",
            fontsize=14,
            y=0.97
        )

        output_file = output_dir / (
            f"{safe_filename(args.prefix)}_"
            f"{safe_filename(gene)}_featureplot.png"
        )

        fig.savefig(
            output_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close(fig)

        print(f"Saved: {output_file}")

    print("Done.")


if __name__ == "__main__":
    main()
