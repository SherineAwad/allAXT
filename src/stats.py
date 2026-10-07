import argparse
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

os.makedirs("figures", exist_ok=True)

# =========================
# COUNTS TABLE
# =========================
ct_counts = pd.crosstab(adata.obs["celltype"], adata.obs["sample"])
ct_counts["Total"] = ct_counts.sum(axis=1)
ct_counts.loc["Total"] = ct_counts.sum(axis=0)

ct_counts.to_csv(f"{args.prefix}_cell_counts.csv")

# =========================
# FRACTIONS FOR PLOT
# =========================
ratios_df = (
    ct_counts.drop("Total")
    .drop(columns="Total")
    .div(ct_counts.drop("Total").drop(columns="Total").sum(axis=0), axis=1)
)

# =========================
# 🔥 COLOR FIX (NO MORE DUPLICATES)
# =========================
celltypes = ratios_df.index.tolist()
n_ct = len(celltypes)

# Use a large distinct palette (tab20 + tab20b + tab20c combined)
palette = (
    list(plt.cm.tab20.colors)
    + list(plt.cm.tab20b.colors)
    + list(plt.cm.tab20c.colors)
)

# If still more celltypes, fallback to evenly spaced HSV
if n_ct > len(palette):
    extra = plt.cm.hsv(np.linspace(0, 1, n_ct - len(palette)))
    palette.extend(extra)

colors = palette[:n_ct]

# =========================
# PLOT
# =========================
ax = ratios_df.T.plot(
    kind='bar',
    stacked=True,
    figsize=(12, 6),
    color=colors,
    width=0.9
)

plt.ylabel("Fraction of cells")
plt.xlabel("Sample")
plt.title("Cell Type Contribution per Sample")
plt.xticks(rotation=45, ha='right')

# cleaner legend (sorted by abundance)
handles, labels = ax.get_legend_handles_labels()
plt.legend(
    handles[::-1],
    labels[::-1],
    title='Cell Type',
    bbox_to_anchor=(1.05, 1),
    loc='upper left'
)

plt.tight_layout()

plt.savefig(
    f"figures/{args.prefix}_cell_ratios.png",
    dpi=600,
    bbox_inches='tight'
)
plt.close()

print(f"Saved: {args.prefix}_cell_counts.csv")
print(f"Saved: figures/{args.prefix}_cell_ratios.png")
