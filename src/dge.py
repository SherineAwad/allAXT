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
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--N", type=int, required=True)  # heatmap size
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--group", required=True)
parser.add_argument("--reference", required=True)
args = parser.parse_args()

os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
adata = sc.read_h5ad(args.input)

# ----------------------------
# SUBSET TO GROUP VS REFERENCE
# ----------------------------
adata = adata[adata.obs["sample"].isin([args.group, args.reference])].copy()

# ----------------------------
# GLOBAL DGE (UNCHANGED)
# ----------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="sample",
    reference=args.reference,
    method="wilcoxon",
    layer="log1p",
    use_raw=False
)

df = sc.get.rank_genes_groups_df(adata, group=None)

df = df[["names", "logfoldchanges", "pvals_adj"]].rename(
    columns={"names": "gene"}
)

# Save full DGE (all genes)
full_out = f"{args.group}_{args.reference}_dge.csv"
df.to_csv(full_out, index=False)

# ----------------------------
# FILTER SIGNIFICANT
# ----------------------------
sig = df[df["pvals_adj"] < args.pvalue].copy()
sig_out = f"{args.group}_{args.reference}_dge_{args.pvalue}.csv"
sig.to_csv(sig_out, index=False)

if sig.empty:
    raise ValueError("No significant genes found")

# ============================================================
#  HEATMAP (TOP N BY P-VALUE + DIRECTION)
# ============================================================

top_up_hm = sig[sig["logfoldchanges"] > 0].sort_values(
    ["pvals_adj", "logfoldchanges"], ascending=[True, False]
).head(args.N)

top_down_hm = sig[sig["logfoldchanges"] < 0].sort_values(
    ["pvals_adj", "logfoldchanges"], ascending=[True, True]
).head(args.N)

top_heat = pd.concat([top_up_hm, top_down_hm])

heatmap_df = top_heat.set_index("gene")[["logfoldchanges"]]

plt.figure(figsize=(4, max(3, len(heatmap_df) * 0.12)))


sns.heatmap(
    heatmap_df,
    cmap="coolwarm",
    center=0,
    annot=False,
    cbar_kws={"label": "logFC"}
)

plt.title(f"Global DE: Reg vs {args.reference}")
plt.ylabel("Genes")
plt.xlabel("")

plt.tight_layout()
plt.savefig(
    f"figures/{args.group}_{args.reference}_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ----------------------------
# SAVE H5AD (UNCHANGED)
# ----------------------------
adata.write_h5ad(args.output)

