import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import argparse
import os

# ----------------------------
# ARGS
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--N", type=int, default=25)
parser.add_argument("--pvalue", type=float, default=0.05)
parser.add_argument("--celltype", required=True)  # Two cell types separated by comma
parser.add_argument("--sample", default=None)
args = parser.parse_args()

groups = [x.strip() for x in args.celltype.split(",")]

if len(groups) != 2 or not all(groups):
    raise ValueError("--celltype must contain exactly two cell types separated by a comma")

group1, group2 = groups

os.makedirs("figures", exist_ok=True)

# ----------------------------
# LOAD DATA
# ----------------------------
adata = sc.read_h5ad(args.input)

if "celltype" not in adata.obs.columns:
    raise ValueError("'celltype' column not found in adata.obs")

# ----------------------------
# OPTIONAL SAMPLE FILTER
# ----------------------------
if args.sample is not None:
    if "sample" not in adata.obs.columns:
        raise ValueError("'sample' column not found in adata.obs")

    if args.sample not in adata.obs["sample"].astype(str).unique():
        raise ValueError(f"Sample '{args.sample}' not found")

    adata = adata[
        adata.obs["sample"].astype(str) == args.sample
    ].copy()

    print(f"Using sample: {args.sample}")
else:
    print("Using all samples.")

# ----------------------------
# CHECK CELL TYPES
# ----------------------------
adata.obs["celltype"] = adata.obs["celltype"].astype(str)

available = adata.obs["celltype"].unique()

if group1 not in available:
    raise ValueError(f"Cell type '{group1}' not found")

if group2 not in available:
    raise ValueError(f"Cell type '{group2}' not found")

adata = adata[
    adata.obs["celltype"].isin([group1, group2])
].copy()

print(f"Comparing {group1} vs {group2}")
print(f"Cells in {group1}: {(adata.obs['celltype'] == group1).sum()}")
print(f"Cells in {group2}: {(adata.obs['celltype'] == group2).sum()}")

# ----------------------------
# DGE: GROUP 1 VS GROUP 2
# ----------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="celltype",
    groups=[group1],
    reference=group2,
    method="wilcoxon",
    layer="log1p",
    use_raw=False
)

df = sc.get.rank_genes_groups_df(adata, group=group1)

# ----------------------------
# SAVE FULL DGE
# ----------------------------
full_df = df[
    ["names", "logfoldchanges", "pvals_adj", "scores"]
].rename(columns={"names": "gene"})

comparison = f"{group1}_vs_{group2}"

if args.sample is not None:
    comparison += f"_{args.sample}"

full_df.to_csv(
    f"{args.prefix}_{comparison}_dge.csv",
    index=False
)

# ----------------------------
# FILTER SIGNIFICANT
# ----------------------------
df = df[
    ["names", "logfoldchanges", "pvals_adj"]
].rename(columns={"names": "gene"})

sig = df[df["pvals_adj"] < args.pvalue].copy()

sig.to_csv(
    f"{args.prefix}_{comparison}_dge_{args.pvalue}.csv",
    index=False
)

if sig.empty:
    raise ValueError("No significant genes found")

# ----------------------------
# HEATMAP: TOP N UP AND DOWN
# ----------------------------
top_up = sig[sig["logfoldchanges"] > 0].sort_values(
    "logfoldchanges", ascending=False
).head(args.N)

top_down = sig[sig["logfoldchanges"] < 0].sort_values(
    "logfoldchanges", ascending=True
).head(args.N)

top_heat = pd.concat([top_up, top_down])

if top_heat.empty:
    raise ValueError("No genes found for heatmap")

heatmap_df = top_heat.set_index("gene")[["logfoldchanges"]]

plt.figure(figsize=(4, max(3, len(heatmap_df) * 0.12)))

sns.heatmap(
    heatmap_df,
    cmap="coolwarm",
    center=0,
    annot=False,
    cbar_kws={"label": "logFC"}
)

plt.title(f"DGE: {group1} vs {group2}")
plt.ylabel("Genes")
plt.xlabel("")

plt.tight_layout()
plt.savefig(
    f"figures/{args.prefix}_{comparison}_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("DGE analysis completed.")

