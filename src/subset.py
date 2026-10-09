import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Input h5ad file")
parser.add_argument("--output", required=True, help="Output h5ad file")

parser.add_argument(
    "--celltype",
    default=None,
    help="Celltype(s) to subset (comma-separated)"
)

parser.add_argument(
    "--sample",
    default=None,
    help="Sample(s) to subset (comma-separated)"
)

parser.add_argument("--prefix", default="output", help="Prefix for saved figures")

args = parser.parse_args()

# -------------------------
# Read data
# -------------------------
adata = sc.read_h5ad(args.input)

# -------------------------
# Validation
# -------------------------
if "celltype" not in adata.obs.columns:
    raise ValueError("celltype column not found in adata.obs")

if "sample" not in adata.obs.columns:
    raise ValueError("sample column not found in adata.obs")

# Validate selections against the original data
if args.celltype is not None:
    celltypes = [ct.strip() for ct in args.celltype.split(",")]

    invalid = [ct for ct in celltypes if ct not in adata.obs["celltype"].unique()]
    if invalid:
        available = adata.obs["celltype"].unique()
        raise ValueError(f"Celltype(s) '{invalid}' not found. Available: {available}")

if args.sample is not None:
    samples = [s.strip() for s in args.sample.split(",")]

    invalid = [s for s in samples if s not in adata.obs["sample"].unique()]
    if invalid:
        available = adata.obs["sample"].unique()
        raise ValueError(f"Sample(s) '{invalid}' not found. Available: {available}")

# -------------------------
# Subset by celltype (if provided)
# -------------------------
if args.celltype is not None:
    adata = adata[adata.obs["celltype"].isin(celltypes)].copy()

# -------------------------
# Subset by sample (if provided)
# -------------------------
if args.sample is not None:
    adata = adata[adata.obs["sample"].isin(samples)].copy()

if adata.n_obs == 0:
    raise ValueError("No cells remain after subsetting. Check the celltype and sample selections.")

print(f"Subset result: {adata.n_obs} cells, {adata.n_vars} genes")

# -------------------------
# Create output dir
# -------------------------
os.makedirs("figures", exist_ok=True)

# -------------------------
# GLOBAL UMAP (subset only)
# -------------------------
fig, ax = plt.subplots(figsize=(6, 5))

sc.pl.umap(
    adata,
    color="sample",
    size=20,
    ax=ax,
    show=False
)

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_umap.png", dpi=300)
plt.close()

# -------------------------
# PER-SAMPLE UMAP
# -------------------------
samples = adata.obs["sample"].unique()

fig, axes = plt.subplots(1, len(samples), figsize=(5 * len(samples), 5))

# handle single sample case
if len(samples) == 1:
    axes = [axes]

for i, s in enumerate(samples):
    sc.pl.umap(
        adata[adata.obs["sample"] == s],
        color="sample",
        size=20,
        title=f"{s}",
        ax=axes[i],
        show=False
    )

plt.tight_layout()
plt.savefig(f"figures/{args.prefix}_perSample_umap.png", dpi=300)
plt.close()

# -------------------------
# Save subset
# -------------------------
adata.write_h5ad(args.output)
print(f"Saved: {args.output}")

# -------------------------
# Summary
# -------------------------
print("\nSummary:")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")
print(f"  Samples: {adata.obs['sample'].value_counts().to_dict()}")
