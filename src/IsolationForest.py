#!/usr/bin/env python3
# lineage_fidelity_condition_subplots.py
# One-class reference model: train on Osteosarcoma in nonReg, score every cell
# by how reference-like it is. Outputs condition-faceted boxplots.

import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import IsolationForest
from scipy.sparse import issparse


def to_dense(X):
    return X.toarray() if issparse(X) else X


def safe_copy_adata(adata, idx):
    try:
        return adata[idx].copy()
    except Exception:
        X = adata.X[idx]
        X = X.copy() if issparse(X) else np.array(X)
        obs = adata.obs.iloc[idx].copy()
        var = adata.var.copy()
        obsm = {}
        if hasattr(adata, 'obsm'):
            for key, value in adata.obsm.items():
                try:
                    obsm[key] = value[idx]
                except Exception:
                    obsm[key] = value
        return sc.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=adata.uns)


def main():
    parser = argparse.ArgumentParser(description="Reference-likeness scores and condition-faceted boxplots")
    parser.add_argument("--input", required=True)
    parser.add_argument("--reference", default="Osteosarcoma")
    parser.add_argument("--condition", default="nonReg")
    parser.add_argument("--prefix", default="")
    parser.add_argument("--hvg_genes", type=int, default=3000)
    parser.add_argument("--n_estimators", type=int, default=100)
    parser.add_argument("--max_cells", type=int, default=3000)
    parser.add_argument("--contamination", type=float, default=0.1)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--chunk_size", type=int, default=5000)
    parser.add_argument("--n_jobs", type=int, default=4)

    args = parser.parse_args()

    os.makedirs("figures", exist_ok=True)
    prefix_str = f"{args.prefix}_" if args.prefix else ""

    print("Loading data...")
    adata = sc.read_h5ad(args.input)

    # Reference cells
    base_mask = (adata.obs["celltype"] == args.reference) & (adata.obs["sample"] == args.condition)
    base_idx = np.where(base_mask)[0]
    if len(base_idx) < 10:
        raise RuntimeError(f"Not enough reference cells for '{args.reference}' in '{args.condition}'.")

    # Build feature space (HVGs) from reference cells only
    print("Selecting highly variable genes from reference cells...")
    ref_anchor = safe_copy_adata(adata, base_idx)
    Xc = to_dense(ref_anchor.X)
    keep = np.var(Xc, axis=0) > 0
    ref_anchor = ref_anchor[:, keep].copy()

    hvg_n = min(args.hvg_genes, ref_anchor.n_vars)
    try:
        sc.pp.highly_variable_genes(ref_anchor, n_top_genes=hvg_n, flavor="seurat", subset=True)
    except Exception:
        var_idx = np.argsort(np.var(to_dense(ref_anchor.X), axis=0))[-hvg_n:]
        ref_anchor = ref_anchor[:, var_idx].copy()

    var_names = ref_anchor.var_names

    # Train one-class model on reference cells only
    rng = np.random.default_rng(args.random_state)
    ref_sub = safe_copy_adata(adata, base_idx)[:, var_names]
    if ref_sub.n_obs > args.max_cells:
        ref_sub = ref_sub[rng.choice(ref_sub.n_obs, args.max_cells, replace=False)].copy()

    X_ref = np.nan_to_num(to_dense(ref_sub.X))
    print(f"Training one-class model on {X_ref.shape[0]} reference cells...")
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("iforest", IsolationForest(
            n_estimators=max(args.n_estimators, 100),
            contamination=args.contamination,
            random_state=args.random_state,
            n_jobs=args.n_jobs,
        )),
    ])
    model.fit(X_ref)

    # Score all cells (chunked)
    print("Scoring all cells against reference model (chunked)...")
    n = adata.n_obs
    ref_likeness = np.empty(n, dtype=np.float32)
    chunk = args.chunk_size
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        idx = np.arange(start, stop)
        Xc = adata[idx, var_names].X
        Xc = Xc.toarray() if issparse(Xc) else np.asarray(Xc)
        Xc = np.nan_to_num(Xc).astype(np.float32)
        
        # raw score_samples: ~ -0.5 (anomaly/different) to 0.0 (normal/blends)
        # Convert to 0-1 scale: 1 = blends with reference (hard to isolate), 0 = distinct (easy to isolate)
        raw_scores = model.score_samples(Xc).astype(np.float32)
        ref_likeness[idx] = (raw_scores + 0.5) * 2
        
        if start % 50000 == 0:
            print(f"  scored {stop}/{n}")

    fidelity_df = pd.DataFrame(index=adata.obs_names)
    fidelity_df["sample"] = adata.obs["sample"].values
    fidelity_df["celltype"] = adata.obs["celltype"].values
    fidelity_df["Fidelity_Reference"] = ref_likeness

    melted_df = fidelity_df.melt(
        id_vars=["sample", "celltype"],
        value_vars=["Fidelity_Reference"],
        var_name="Lineage",
        value_name="Fidelity",
    )

    g = sns.catplot(
        data=melted_df,
        x="celltype",
        y="Fidelity",
        hue="celltype",
        col="sample",
        kind="box",
        sharey=True,
        height=5,
        aspect=1.2,
        palette="tab20",
        legend=False,
    )
    g.set_axis_labels("Cell type", "Reference-likeness (high = blends with reference)")
    g.set_titles(col_template="Condition: {col_name}")
    g.set_xticklabels(rotation=45, horizontalalignment='right')
    g.fig.subplots_adjust(top=0.85)
    g.fig.suptitle(f"Reference-likeness to {args.reference} ({args.condition}) per cell type and condition")

    box_plot_path = os.path.join("figures", f"{prefix_str}osteosarcoma_lineage_condition_subplots.png")
    g.savefig(box_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved condition-faceted boxplot to {box_plot_path}")
    print("Done!")


if __name__ == "__main__":
    main()
