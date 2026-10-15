#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import palantir


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def save_trajectory_plots(
    adata,
    branch_name,
    pseudotime_key,
    masks_key,
    figures_dir,
    prefix,
    root_safe,
):
    """Save one trajectory figure showing only the specified branch."""

    fig, ax = plt.subplots(figsize=(8, 8))

    palantir.plot.plot_trajectory(
        adata,
        branch=branch_name,
        pseudo_time_key=pseudotime_key,
        masks_key=masks_key,
        embedding_basis="X_umap",
        cell_color=pseudotime_key,
        n_arrows=10,
        ax=ax,
        figsize=(8, 8),
    )

    ax.set_title(
        f"Palantir trajectory: root {root_safe} ({branch_name})"
    )

    plot_path = figures_dir / (
        f"{prefix}_palantir_root_{root_safe}.png"
    )

    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] Saved Palantir trajectory plot: {plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run Palantir pseudotime analysis on an AnnData dataset."
    )

    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument(
        "--output",
        required=True,
        help="Output .h5ad filename, saved in the current directory",
    )
    parser.add_argument("--prefix", required=True, help="Output filename prefix")
    parser.add_argument(
        "--cluster_key",
        default="celltype",
        help="obs column defining candidate root groups (default: celltype)",
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Leiden cluster to use as the root, e.g. 0, 1, 4, or 5",
    )

    args = parser.parse_args()

    output_dir = Path.cwd()
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    output_h5ad = output_dir / Path(args.output).name

    print(
        f"[INFO] Palantir version: "
        f"{getattr(palantir, '__version__', 'unknown')}"
    )
    print(f"[INFO] Reading: {args.input}")

    adata = ad.read_h5ad(args.input)

    print(
        f"[INFO] Dataset: "
        f"{adata.n_obs:,} cells x {adata.n_vars:,} genes"
    )

    if args.cluster_key not in adata.obs.columns:
        raise KeyError(
            f"'{args.cluster_key}' is not present in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Validate the requested root cluster.
    root_name = str(args.root)

    available_roots = (
        adata.obs[args.cluster_key]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if root_name not in available_roots:
        raise ValueError(
            f"Root cluster '{root_name}' was not found in "
            f"adata.obs['{args.cluster_key}']. "
            f"Available clusters: {available_roots}"
        )

    original_X = adata.X.copy() if hasattr(adata.X, "copy") else adata.X

    if "log1p" in adata.layers:
        adata.X = adata.layers["log1p"].copy()
        print("[INFO] Using adata.layers['log1p'] as adata.X")
    else:
        print(
            "[INFO] No 'log1p' layer found; "
            "retaining the existing adata.X"
        )

    if "X_pca" not in adata.obsm:
        print("[INFO] X_pca not found; computing PCA")
        sc.pp.pca(adata)
    else:
        print("[INFO] Using existing adata.obsm['X_pca']")

    n_components = min(10, adata.obsm["X_pca"].shape[1])

    if n_components < 3:
        raise ValueError(
            f"At least 3 PCA components are required; found {n_components}."
        )

    print("[INFO] Computing Palantir diffusion maps")
    palantir.utils.run_diffusion_maps(
        adata,
        n_components=n_components,
    )

    print("[INFO] Determining Palantir multiscale space")
    palantir.utils.determine_multiscale_space(adata)

    if "DM_EigenVectors_multiscaled" not in adata.obsm:
        raise RuntimeError(
            "Palantir did not create "
            "adata.obsm['DM_EigenVectors_multiscaled']."
        )

    print(
        "[INFO] Palantir multiscale space shape: "
        f"{adata.obsm['DM_EigenVectors_multiscaled'].shape}"
    )

    if "X_umap" not in adata.obsm:
        print("[INFO] X_umap not found; computing UMAP for plotting")
        sc.pp.neighbors(
            adata,
            use_rep="DM_EigenVectors_multiscaled",
        )
        sc.tl.umap(adata, random_state=42)
    else:
        print("[INFO] Using existing adata.obsm['X_umap'] for plotting")

    # Process only the requested root cluster.
    root_safe = safe_name(root_name)

    print("\n" + "=" * 70)
    print(f"[INFO] Processing root group: {root_name}")
    print("=" * 70)

    root_mask = adata.obs[args.cluster_key].astype(str) == root_name
    root_cells = adata.obs_names[root_mask]

    if len(root_cells) == 0:
        raise ValueError(
            f"No cells found for root '{root_name}' "
            f"in '{args.cluster_key}'."
        )

    print(f"[INFO] Root group cells: {len(root_cells):,}")

    # Identify an early cell within the requested cluster.
    try:
        early_cell = palantir.utils.early_cell(
            adata,
            celltype=root_name,
            celltype_column=args.cluster_key,
        )

        if early_cell not in root_cells:
            raise ValueError(
                f"Selected early cell {early_cell} is not in "
                f"root cluster {root_name}."
            )

        print(f"[INFO] Automatically determined early cell: {early_cell}")

    except Exception as e:
        print(
            f"[WARNING] Automated early_cell detection failed ({e}); "
            "falling back to the first cell in the requested root cluster."
        )
        early_cell = root_cells[0]

    print(f"[INFO] Starting cell: {early_cell}")
    print(
        f"[INFO] Starting cell Leiden cluster: "
        f"{adata.obs.loc[early_cell, args.cluster_key]}"
    )

    # Explicitly force Palantir to use the selected early cell as its start.
    pr_res = palantir.core.run_palantir(
        adata,
        early_cell=early_cell,
        num_waypoints=min(1200, max(1, adata.n_obs // 10)),
        max_iterations=25,
        seed=42,
        use_early_cell_as_start=True,
    )

    pseudotime = pd.Series(pr_res.pseudotime).reindex(adata.obs_names)
    entropy = pd.Series(pr_res.entropy).reindex(adata.obs_names)
    branch_probs = pr_res.branch_probs.reindex(adata.obs_names)

    if branch_probs.shape[1] == 0:
        raise RuntimeError(
            f"Palantir returned no fate-probability columns for root "
            f"'{root_name}'."
        )

    pseudotime_key = f"palantir_pseudotime_{root_safe}"
    entropy_key = f"palantir_entropy_{root_safe}"
    branch_key = f"palantir_branch_selected_{root_safe}"
    fate_prob_key = f"palantir_fate_probabilities_{root_safe}"
    masks_key = f"branch_masks_{root_safe}"

    adata.obs[pseudotime_key] = pseudotime
    adata.obs[entropy_key] = entropy
    adata.obsm[fate_prob_key] = branch_probs

    print(
        f"[INFO] Pseudotime range: "
        f"{pseudotime.min():.4f}–{pseudotime.max():.4f}"
    )
    print(
        f"[INFO] Entropy range: "
        f"{entropy.min():.4f}–{entropy.max():.4f}"
    )
    print(f"[INFO] Fate probability matrix shape: {branch_probs.shape}")
    print(f"[INFO] Fate columns: {list(branch_probs.columns)}")
    print("[INFO] Fate probability summary:")
    print(branch_probs.describe().to_string())

    row_sums = branch_probs.sum(axis=1)
    print(
        "[INFO] Fate probability row sums: "
        f"min={row_sums.min():.6f}, "
        f"median={row_sums.median():.6f}, "
        f"max={row_sums.max():.6f}"
    )

    print("[INFO] Selecting branch cells")
    branch_masks = palantir.presults.select_branch_cells(
        adata,
        pseudo_time_key=pseudotime_key,
        fate_prob_key=fate_prob_key,
        q=0.01,
        eps=0.01,
        masks_key=masks_key,
        save_as_df=False,
    )

    branch_masks = np.asarray(branch_masks, dtype=bool)

    if branch_masks.ndim == 1:
        branch_masks = branch_masks.reshape(-1, 1)

    if branch_masks.shape[0] != adata.n_obs:
        raise RuntimeError(
            "Branch mask cell count does not match AnnData: "
            f"{branch_masks.shape[0]} versus {adata.n_obs}."
        )

    all_branch_names = list(branch_probs.columns)

    # Select the dominant fate by mean probability within the requested root.
    dominant_branch = branch_probs.loc[root_cells].mean(axis=0).idxmax()

    print(
        f"[INFO] Most dominant branch for root '{root_name}': "
        f"{dominant_branch}"
    )

    # Store selection for the dominant branch only.
    dom_idx = all_branch_names.index(dominant_branch)

    if dom_idx >= branch_masks.shape[1]:
        raise RuntimeError(
            f"Dominant branch index {dom_idx} is outside the branch-mask "
            f"matrix with {branch_masks.shape[1]} columns."
        )

    adata.obs[branch_key] = branch_masks[:, dom_idx]

    result_df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "root_group": root_name,
            "early_cell": str(early_cell),
            "pseudotime": pseudotime.to_numpy(),
            "entropy": entropy.to_numpy(),
            "branch_selected": adata.obs[branch_key].to_numpy(),
        }
    ).set_index("cell_id")

    for i, branch_name in enumerate(all_branch_names):
        result_df[f"fate_prob_{branch_name}"] = (
            branch_probs[branch_name].to_numpy()
        )
        result_df[f"branch_mask_{branch_name}"] = (
            branch_masks[:, i]
            if i < branch_masks.shape[1]
            else False
        )

    csv_path = output_dir / (
        f"{args.prefix}_palantir_root_{root_safe}.csv"
    )
    result_df.to_csv(csv_path)
    print(f"[INFO] Saved per-root results: {csv_path}")

    mask_columns = all_branch_names[:branch_masks.shape[1]]
    mask_df = pd.DataFrame(
        branch_masks,
        index=adata.obs_names,
        columns=mask_columns,
    )

    mask_path = output_dir / (
        f"{args.prefix}_branch_masks_{root_safe}.csv"
    )
    mask_df.to_csv(mask_path)
    print(f"[INFO] Saved branch masks: {mask_path}")

    save_trajectory_plots(
        adata,
        dominant_branch,
        pseudotime_key,
        masks_key,
        figures_dir,
        args.prefix,
        root_safe,
    )

    summary = result_df[
        [
            "root_group",
            "early_cell",
            "pseudotime",
            "entropy",
            "branch_selected",
        ]
    ].copy()

    summary.index.name = "cell_id"
    combined_df = summary.reset_index()

    combined_path = output_dir / (
        f"{args.prefix}_palantir_comparison.csv"
    )
    combined_df.to_csv(combined_path, index=False)
    print(f"[INFO] Saved combined results: {combined_path}")

    # Restore original adata.X before writing the output file.
    adata.X = original_X

    adata.write_h5ad(output_h5ad)
    print(f"[INFO] Saved annotated AnnData: {output_h5ad}")


if __name__ == "__main__":
    main()
