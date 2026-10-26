#!/usr/bin/env python3

import argparse
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import ot

from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

warnings.filterwarnings("ignore")


def get_X(adata, layer):
    """Return a dense expression matrix from the requested layer or adata.X."""
    X = adata.layers[layer] if layer and layer in adata.layers else adata.X
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def pot_cost(A, B, min_cells, reg):
    """
    Compute pairwise exact Earth Mover's Distance (EMD) cost.
    Returns the raw OT cost without converting it to similarity.
    No Sinkhorn fallback is used.
    """
    if A is None or B is None:
        return np.nan

    if len(A) < min_cells or len(B) < min_cells:
        return np.nan

    # Compute Euclidean distance matrix in PCA space
    M = ot.dist(A, B, metric="euclidean")

    a = np.ones(len(A)) / len(A)
    b = np.ones(len(B)) / len(B)

    # Exact Earth Mover's Distance (EMD)
    res = ot.emd2(a, b, M)

    cost = float(res)

    return cost if np.isfinite(cost) else np.nan


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


def cluster_order(distance_matrix):
    """
    Get a hierarchical clustering order from complete, finite OT costs.
    If any pairwise distances are missing, skip clustering rather than
    imputing distances that could distort the dendrogram.
    """
    n = distance_matrix.shape[0]

    if n <= 1:
        return np.arange(n), None

    D = np.array(distance_matrix, dtype=float, copy=True)
    D = (D + D.T) / 2
    np.fill_diagonal(D, 0.0)

    if not np.all(np.isfinite(D)):
        print(
            "Warning: Missing or non-finite OT costs detected; "
            "dendrogram unavailable. Using original label order."
        )
        return np.arange(n), None

    D = np.maximum(D, 0.0)

    condensed = squareform(D, checks=False)
    tree = linkage(condensed, method="average")
    leaves = dendrogram(tree, no_plot=True)["leaves"]

    return np.asarray(leaves), tree


def plot_triangular_heatmap(
    matrix,
    labels,
    output_path,
    title,
    annotate=True,
):
    """
    Plot the lower triangle of a symmetric OT cost matrix,
    preserving numerical OT costs in every visible cell.
    """
    n = len(labels)

    if n == 0:
        print("Warning: No states to plot.")
        return

    order, tree = cluster_order(matrix)

    ordered_matrix = matrix[np.ix_(order, order)]
    ordered_labels = [labels[i] for i in order]

    mask = np.triu(np.ones((n, n), dtype=bool), k=0)
    plotted_matrix = np.ma.array(ordered_matrix, mask=mask)

    finite = ordered_matrix[np.isfinite(ordered_matrix)]
    # Anchor vmin at 0.0 so absolute distance zero is dark purple
    vmin = 0.0
    vmax = float(np.max(finite)) if len(finite) else 1.0

    if vmax == vmin:
        vmax = vmin + 1e-12

    fig_size = max(7, min(22, n * 0.55))

    if tree is not None:
        fig = plt.figure(figsize=(fig_size + 3, fig_size + 2.5))

        grid = fig.add_gridspec(
            2,
            2,
            width_ratios=[1.4, 5],
            height_ratios=[1.7, 5],
            left=0.08,
            right=0.87,
            bottom=0.12,
            top=0.82,
            wspace=0.04,
            hspace=0.08,
        )

        ax_top = fig.add_subplot(grid[0, 1])
        ax_left = fig.add_subplot(grid[1, 0])
        ax_heat = fig.add_subplot(grid[1, 1])

        dendrogram(
            tree,
            ax=ax_top,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="black",
        )
        ax_top.set_axis_off()

        dendrogram(
            tree,
            ax=ax_left,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color="black",
        )
        ax_left.set_axis_off()

    else:
        fig, ax_heat = plt.subplots(figsize=(fig_size, fig_size))
        fig.subplots_adjust(
            left=0.25,
            bottom=0.25,
            top=0.85,
            right=0.88,
        )

    cmap = plt.colormaps["viridis"].copy()
    cmap.set_bad(color="white", alpha=0)

    ax_heat.set_facecolor("white")

    image = ax_heat.imshow(
        plotted_matrix,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        origin="upper",
    )

    fig.colorbar(
        image,
        ax=ax_heat,
        label="Optimal Transport cost (lower = more similar)",
        fraction=0.046,
        pad=0.04,
    )

    ax_heat.set_xticks(range(n))
    ax_heat.set_yticks(range(n))

    ax_heat.set_xticklabels(
        ordered_labels,
        rotation=90,
        fontsize=8,
    )
    ax_heat.set_yticklabels(
        ordered_labels,
        fontsize=8,
    )

    if annotate:
        midpoint = (vmin + vmax) / 2

        for i in range(n):
            for j in range(i):
                value = ordered_matrix[i, j]

                if np.isfinite(value):
                    color = (
                        "white" if value > midpoint else "black"
                    )

                    ax_heat.text(
                        j,
                        i,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=6 if n > 40 else 8,
                        color=color,
                    )

    fig.suptitle(
        title,
        y=0.97,
        fontsize=13,
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    print(f"Saved triangular heatmap: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute pairwise Optimal Transport costs across single-cell "
            "celltype-sample groups and optionally plot a focused "
            "triangular heatmap."
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

    # Concatenate all loaded AnnData objects into one for joint feature space PCA
    print("Concatenating AnnData objects...")
    adata = sc.concat(adatas, join="inner")

    # Dynamic composite subtype grouping: 'Osteosarcoma.cluster0'
    adata.obs["subtype_label"] = (
        adata.obs[args.celltype].astype(str)
        + ".cluster"
        + adata.obs[args.cluster].astype(str)
    )

    # Filter out groups below cell count threshold
    counts = adata.obs.groupby(["subtype_label", args.sample], observed=True).size()
    valid_groups = counts[counts >= args.min_cells_per_celltype].index.tolist()

    if not valid_groups:
        raise ValueError("No groups met the minimum cell threshold.")

    # Apply focus selection if provided
    if args.focus:
        focus_list = read_focus_file(args.focus)
        valid_groups = [g for g in valid_groups if g in focus_list]

    if not valid_groups:
        raise ValueError("No matching groups remained after applying focus selection.")

    # Highly variable genes and PCA computation using robust 'seurat' flavor
    if args.hvg and adata.n_vars > args.hvg:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=args.hvg,
            flavor="seurat",
            layer=args.layer if args.layer in adata.layers else None
        )
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    else:
        adata_hvg = adata.copy()

    X_mat = get_X(adata_hvg, args.layer)
    pca = PCA(n_components=min(args.n_pcs, X_mat.shape[1]), random_state=args.seed)
    X_pca = pca.fit_transform(X_mat)

    # Extract PCA coordinates per distinct subtype-sample group
    group_data = {}
    for stype, smp in valid_groups:
        mask = (adata.obs["subtype_label"] == stype) & (adata.obs[args.sample] == smp)
        indices = np.where(mask)[0]

        if len(indices) > args.max_cells:
            indices = np.random.choice(indices, size=args.max_cells, replace=False)

        group_data[f"{stype}.{smp}"] = X_pca[indices]

    labels = sorted(group_data.keys())
    n = len(labels)
    cost_matrix = np.full((n, n), np.nan)

    print(f"Computing pairwise OT costs across {n} state groups...")
    for i in range(n):
        cost_matrix[i, i] = 0.0
        for j in range(i):
            try:
                c = pot_cost(
                    group_data[labels[i]],
                    group_data[labels[j]],
                    args.min_cells,
                    args.reg,
                )
                if not np.isfinite(c):
                    print(
                        f"Warning: OT cost unavailable for "
                        f"{labels[i]} vs {labels[j]}."
                    )
                cost_matrix[i, j] = c
                cost_matrix[j, i] = c

            except Exception as exc:
                print(
                    f"Warning: Exact EMD failed for "
                    f"{labels[i]} vs {labels[j]}: {exc}"
                )
                cost_matrix[i, j] = np.nan
                cost_matrix[j, i] = np.nan

    # -------------------------
    # Save Outputs
    # -------------------------
    # CSV goes to current directory
    csv_path = f"{args.prefix}_ot_cost_matrix.csv"
    df_cost = pd.DataFrame(cost_matrix, index=labels, columns=labels)
    df_cost.to_csv(csv_path)
    print(f"Saved CSV matrix: {csv_path}")

    # Figures go to figures/ directory
    os.makedirs("figures", exist_ok=True)
    figure_path = f"figures/{args.prefix}_ot_cost_heatmap.png"

    plot_triangular_heatmap(
        cost_matrix,
        labels,
        figure_path,
        title=f"Subtype OT Cost Matrix ({args.prefix})",
    )
    print("Done!")


if __name__ == "__main__":
    main()
