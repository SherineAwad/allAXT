
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import argparse
import os

# ----------------------------
# ARGS
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--N", type=int, required=True)  # heatmap size
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--group", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--sample", default=None)
args = parser.parse_args()

os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
adata = sc.read_h5ad(args.input)

# ----------------------------
# OPTIONAL SAMPLE FILTER
# ----------------------------
if args.sample is not None:
    if "sample" not in adata.obs.columns:
        raise ValueError("'sample' column not found in adata.obs")

    if args.sample not in adata.obs["sample"].astype(str).unique():
        raise ValueError(
            f"Sample '{args.sample}' not found. "
            f"Available samples: {adata.obs['sample'].astype(str).unique()}"
        )

    adata = adata[
        adata.obs["sample"].astype(str) == args.sample
    ].copy()

    print(f"Using sample: {args.sample}")
else:
    print("Using all cells.")

# ----------------------------
# FIX: TYPE MATCHING FOR GROUP
# ----------------------------
col_dtype = adata.obs[args.celltype].dtype
is_cat_int = (
    isinstance(col_dtype, pd.CategoricalDtype)
    and pd.api.types.is_integer_dtype(adata.obs[args.celltype].cat.categories)
)

if args.group.isdigit() and (
    pd.api.types.is_integer_dtype(col_dtype) or is_cat_int
):
    args.group = int(args.group)

# ----------------------------
# GLOBAL DGE: GROUP VS REST
# ----------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby=args.celltype,
    groups=[args.group],
    reference="rest",
    method="wilcoxon",
    layer="log1p",
    use_raw=False
)

df = sc.get.rank_genes_groups_df(adata, group=args.group)

# Save full DGE (all genes), including Wilcoxon scores
full_df = df[["names", "logfoldchanges", "pvals_adj", "scores"]].rename(
    columns={"names": "gene"}
)

full_out = f"{args.prefix}_{args.group}_vs_rest_dge.csv"
full_df.to_csv(full_out, index=False)

# Keep the original columns for all downstream steps
df = df[["names", "logfoldchanges", "pvals_adj"]].rename(
    columns={"names": "gene"}
)

# ----------------------------
# FILTER SIGNIFICANT
# ----------------------------
sig = df[df["pvals_adj"] < args.pvalue].copy()
sig_out = f"{args.prefix}_{args.group}_vs_rest_dge_{args.pvalue}.csv"
sig.to_csv(sig_out, index=False)

if sig.empty:
    raise ValueError("No significant genes found")

# ============================================================
# HEATMAP (TOP N BY LOGFC + DIRECTION)
# ============================================================

top_up_hm = sig[sig["logfoldchanges"] > 0].sort_values(
    "logfoldchanges", ascending=False
).head(args.N)

top_down_hm = sig[sig["logfoldchanges"] < 0].sort_values(
    "logfoldchanges", ascending=True
).head(args.N)

top_heat = pd.concat([top_up_hm, top_down_hm])

if top_heat.empty:
    raise ValueError("No genes found to plot in heatmap after filtering.")

heatmap_df = top_heat.set_index("gene")[["logfoldchanges"]]

if heatmap_df.empty:
    raise ValueError("Heatmap dataframe is empty.")

plt.figure(figsize=(4, max(3, len(heatmap_df) * 0.12)))

sns.heatmap(
    heatmap_df,
    cmap="coolwarm",
    center=0,
    annot=False,
    cbar_kws={"label": "logFC"}
)

plt.title(f"Global DE: {args.group} vs rest")
plt.ylabel("Genes")
plt.xlabel("")

plt.tight_layout()
plt.savefig(
    f"figures/{args.prefix}_{args.group}_vs_rest_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

