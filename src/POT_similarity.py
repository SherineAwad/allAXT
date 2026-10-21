#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import warnings
import ot

warnings.filterwarnings("ignore")


def get_X(adata, layer):
    X = adata.layers[layer] if layer in adata.layers else adata.X
    return X.toarray() if hasattr(X, "toarray") else X


def pot_score(A, B, min_cells, reg, kernel_scale):
    """
    Computes regularized Optimal Transport (Sinkhorn) similarity between two cell populations.
    """
    try:
        if A is None or B is None:
            return np.nan
        if len(A) < min_cells or len(B) < min_cells:
            return np.nan

        # Compute true Euclidean distance matrix (explicit metric)
        M = ot.dist(A, B, metric='euclidean')

        # Uniform weights for distributions
        a = np.ones(len(A)) / len(A)
        b = np.ones(len(B)) / len(B)

        # Use Sinkhorn algorithm and safely extract scalar cost
        res = ot.sinkhorn2(a, b, M, reg=reg)
        cost = float(np.ravel(res)[0])

        # Convert OT cost to similarity using the median pairwise OT cost
        similarity = np.exp(-cost / kernel_scale)

        return float(np.clip(similarity, 0.0, 1.0))

    except Exception as e:
        return np.nan


def pot_cost(A, B, min_cells, reg):
    """
    Computes regularized Optimal Transport (Sinkhorn) cost between two cell populations.
    """
    try:
        if A is None or B is None:
            return np.nan
        if len(A) < min_cells or len(B) < min_cells:
            return np.nan

        # Compute true Euclidean distance matrix (explicit metric)
        M = ot.dist(A, B, metric='euclidean')

        # Uniform weights for distributions
        a = np.ones(len(A)) / len(A)
        b = np.ones(len(B)) / len(B)

        # Use Sinkhorn algorithm and safely extract scalar cost
        res = ot.sinkhorn2(a, b, M, reg=reg)
        cost = float(np.ravel(res)[0])

        return cost

    except Exception as e:
        return np.nan


def main():
    parser = argparse.ArgumentParser(description="Compute scalable Optimal Transport similarity across single-cell groups.")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--prefix", required=True, help="Output prefix")
    parser.add_argument("--layer", default="log1p", help="Layer to use for expression")
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of PCs")
    parser.add_argument("--max_cells", type=int, default=10000, help="Maximum cells to subsample per group")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum cells required per group")
    parser.add_argument("--min_cells_per_celltype", type=int, default=10, help="Minimum total cells required per cell type across all samples")
    parser.add_argument("--hvg", type=int, default=5000, help="Number of highly variable genes to use (default: 2000)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--focus", default=None, help="Comma-separated list of celltypes or states to plot in a separate focused heatmap (e.g., 'Osteosarcoma,Mesenchymal')")
    args = parser.parse_args()

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    # =====================================================
    # FILTER OUT CELL TYPES WITH TOO FEW CELLS
    # =====================================================
    if "celltype" not in adata.obs or "sample" not in adata.obs:
        raise ValueError("adata.obs must contain 'celltype' and 'sample' columns")

    print(f"\nFiltering cell types with fewer than {args.min_cells_per_celltype} total cells across all samples...")

    celltype_counts = adata.obs["celltype"].value_counts()
    keep_celltypes = celltype_counts[
        celltype_counts >= args.min_cells_per_celltype
    ].index

    removed_celltypes = celltype_counts[
        celltype_counts < args.min_cells_per_celltype
    ]

    if len(removed_celltypes) > 0:
        print("Removing these cell types:")
        for ct, n in removed_celltypes.items():
            print(f"  {ct}: {n} cells")
    else:
        print("No cell types removed.")

    adata = adata[adata.obs["celltype"].isin(keep_celltypes)].copy()

    if adata.n_obs == 0:
        raise ValueError(
            f"No cells remain after filtering cell types with fewer than "
            f"{args.min_cells_per_celltype} total cells."
        )

    print(f"Remaining cell types: {adata.obs['celltype'].nunique()}")
    print(f"Remaining cells: {adata.n_obs}")

    # =====================================================
    # FILTER TO HIGHLY VARIABLE GENES
    # =====================================================
    if args.hvg > 0:
        print(f"Filtering to top {args.hvg} highly variable genes...")
        batch_key = 'sample' if 'sample' in adata.obs else None
        sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg, batch_key=batch_key, flavor='seurat')
        adata = adata[:, adata.var.highly_variable]
        print(f"Kept {adata.shape[1]} genes")
    else:
        print("Warning: Using all genes (not recommended for single-cell PCA/OT).")

    # =====================================================
    # FILTER OUT ZERO-CELL COMBINATIONS
    # =====================================================
    print("Filtering out celltype-sample combinations with zero cells...")

    valid_groups = adata.obs.groupby(['celltype', 'sample']).size()
    valid_groups = valid_groups[valid_groups > 0]

    print("Keeping these combinations:")
    for (ct, s), n in valid_groups.items():
        print(f"  {ct}.{s}: {n} cells")

    celltypes = sorted(set([ct for ct, s in valid_groups.index]))
    samples = sorted(set([s for ct, s in valid_groups.index]))

    # =====================================================
    # BUILD PCA SPACE
    # =====================================================
    print("\nBuilding global PCA space on HVGs...")

    all_cells = []
    for ct in celltypes:
        for s in samples:
            sub = adata[(adata.obs["celltype"] == ct) & (adata.obs["sample"] == s)]
            if sub.n_obs > 0:
                all_cells.append(get_X(sub, args.layer))

    if len(all_cells) == 0:
        raise ValueError("No valid cell groups found. Check filters or input data.")

    X_all = np.vstack(all_cells)

    pca = PCA(n_components=min(args.n_pcs, X_all.shape[1], X_all.shape[0]), random_state=args.seed)
    X_pca = pca.fit_transform(X_all)
    print(f"PCA shape: {X_pca.shape}")

    # =====================================================
    # CREATE STATES DICTIONARY
    # =====================================================
    states = {}
    idx = 0
    for ct in celltypes:
        for s in samples:
            sub = adata[(adata.obs["celltype"] == ct) & (adata.obs["sample"] == s)]
            n = sub.n_obs
            if n > 0:
                states[(ct, s)] = X_pca[idx:idx + n]
                idx += n

    keys = sorted(list(states.keys()))
    n_states = len(keys)
    print(f"\nTotal states to compare: {n_states}")

    # =====================================================
    # PRE-SUBSAMPLE STATES DETERMINISTICALLY
    # =====================================================
    subsampled_states = {}
    for idx_i, key in enumerate(keys):
        group_data = states[key]
        if len(group_data) > args.max_cells:
            group_seed = args.seed + idx_i
            group_rng = np.random.default_rng(group_seed)
            sub_idx = group_rng.choice(len(group_data), args.max_cells, replace=False)
            group_data = group_data[sub_idx]
        subsampled_states[key] = group_data

    # =====================================================
    # COMPUTE ADAPTIVE GLOBAL SCALE
    # =====================================================
    print("\nEstimating adaptive global geometric scale from data...")
    pooled_sample_data = np.vstack(list(subsampled_states.values()))
    if len(pooled_sample_data) > 3000:
        sample_indices = np.random.default_rng(args.seed).choice(len(pooled_sample_data), 3000, replace=False)
        pooled_sample_data = pooled_sample_data[sample_indices]

    sample_dist_matrix = ot.dist(pooled_sample_data, pooled_sample_data, metric='euclidean')
    off_diag_indices = np.triu_indices_from(sample_dist_matrix, k=1)
    global_median_dist = float(np.median(sample_dist_matrix[off_diag_indices]))

    reg = max(0.01, 0.05 * global_median_dist)

    print(f"Adaptive Global Median Distance (Euclidean): {global_median_dist:.4f}")
    print(f"Sinkhorn Regularization (reg): {reg:.4f}")

    # =====================================================
    # COMPUTE PAIRWISE OT COST MATRIX
    # =====================================================
    cost_matrix = np.full((n_states, n_states), np.nan)

    print("\nComputing pairwise Sinkhorn OT costs...")
    print(f"Parameters: max_cells={args.max_cells}, min_cells={args.min_cells}, hvg={args.hvg}\n")

    for i, k1 in enumerate(keys):
        for j, k2 in enumerate(keys):
            if i == j:
                cost_matrix[i, j] = 0.0
            elif i < j:
                cost_matrix[i, j] = pot_cost(
                    subsampled_states[k1],
                    subsampled_states[k2],
                    args.min_cells,
                    reg=reg
                )
                cost_matrix[j, i] = cost_matrix[i, j]
        print(f"  done: {k1[0]}.{k1[1]}")

    # =====================================================
    # COMPUTE SIMILARITY MATRIX USING MEDIAN OT COST
    # =====================================================
    matrix = np.full((n_states, n_states), np.nan)

    valid_costs = cost_matrix[np.triu_indices(n_states, k=1)]
    valid_costs = valid_costs[np.isfinite(valid_costs)]

    if len(valid_costs) > 0:
        median_ot_cost = float(np.median(valid_costs))
    else:
        median_ot_cost = 1.0

    kernel_scale = median_ot_cost if median_ot_cost > 0 else 1.0

    print(f"Median Pairwise OT Cost (reference scale): {median_ot_cost:.4f}")
    print(f"Kernel Scale: {kernel_scale:.4f}")

    for i in range(n_states):
        for j in range(n_states):
            if i == j:
                matrix[i, j] = 1.0
            elif np.isfinite(cost_matrix[i, j]):
                matrix[i, j] = float(
                    np.clip(np.exp(-cost_matrix[i, j] / kernel_scale), 0.0, 1.0)
                )

    # =====================================================
    # PLOT FULL MATRIX
    # =====================================================
    os.makedirs("figures", exist_ok=True)
    labels = [f"{ct}.{s}" for (ct, s) in keys]
    n = len(labels)

    plt.figure(figsize=(max(8, n * 0.4), max(8, n * 0.4)))
    plt.imshow(matrix, cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="Sinkhorn POT similarity")
    plt.xticks(range(n), labels, rotation=90, fontsize=8)
    plt.yticks(range(n), labels, fontsize=8)

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            if not np.isnan(val):
                color = "white" if val < 0.5 else "black"
                plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color=color)

    hvg_text = f"_hvg{args.hvg}" if args.hvg > 0 else ""
    plt.title(f"Sinkhorn OT similarity (Full)\n{args.prefix}")
    plt.tight_layout()

    out = f"figures/{args.prefix}_sinkhorn_pot{hvg_text}_matrix.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\nSaved full matrix: {out}")

    # =====================================================
    # PLOT FOCUSED SUBSET (Optional)
    # =====================================================
    if args.focus:
        focus_list = [f.strip() for f in args.focus.split(",")]
        focus_indices = []
        focus_keys = []

        for idx, (ct, s) in enumerate(keys):
            if ct in focus_list or f"{ct}.{s}" in focus_list:
                focus_indices.append(idx)
                focus_keys.append((ct, s))

        if len(focus_keys) > 0:
            sub_matrix = matrix[np.ix_(focus_indices, focus_indices)]
            labels_focus = [f"{ct}.{s}" for (ct, s) in focus_keys]
            n_focus = len(labels_focus)

            plt.figure(figsize=(max(6, n_focus * 0.5), max(6, n_focus * 0.5)))
            plt.imshow(sub_matrix, cmap="viridis", vmin=0, vmax=1)
            plt.colorbar(label="Sinkhorn POT similarity")
            plt.xticks(range(n_focus), labels_focus, rotation=90, fontsize=8)
            plt.yticks(range(n_focus), labels_focus, fontsize=8)

            for i in range(n_focus):
                for j in range(n_focus):
                    val = sub_matrix[i, j]
                    if not np.isnan(val):
                        color = "white" if val < 0.5 else "black"
                        plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color=color)

            plt.title(f"Sinkhorn OT similarity (Focused: {args.focus})\n{args.prefix}")
            plt.tight_layout()

            out_focus = f"figures/{args.prefix}_sinkhorn_pot{hvg_text}_focus_matrix.png"
            plt.savefig(out_focus, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"Saved focused matrix: {out_focus}")
        else:
            print(f"Warning: None of the focus items {focus_list} matched available keys.")

    print("Done.")


if __name__ == "__main__":
    main()
