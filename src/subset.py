#!/usr/bin/env python3

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
parser.add_argument(
    "--resolution",
    type=float,
    default=0.5,
    help="Leiden clustering resolution"
)

args = parser.parse_args()

# -------------------------
# Read data
# -------------------------
adata = sc.read_h5ad(args.input)

# -------------------------
# Validation
# -------------------------
for column in ["celltype", "sample"]:
    if column not in adata.obs.columns:
        raise ValueError(f"{column} column not found in adata.obs")

if args.celltype is not None:
    celltypes = [ct.strip() for ct in args.celltype.split(",")]

    invalid = [
        ct for ct in celltypes
        if ct not in adata.obs["celltype"].unique()
    ]

    if invalid:
        raise ValueError(
            f"Celltype(s) {invalid} not found. "
            f"Available: {adata.obs['celltype'].unique()}"
        )

if args.sample is not None:
    selected_samples = [s.strip() for s in args.sample.split(",")]

    invalid = [
        s for s in selected_samples
        if s not in adata.obs["sample"].unique()
    ]

    if invalid:
        raise ValueError(
            f"Sample(s) {invalid} not found. "
            f"Available: {adata.obs['sample'].unique()}"
        )

# -------------------------
# Subset by celltype
# -------------------------
if args.celltype is not None:
    adata = adata[
        adata.obs["celltype"].isin(celltypes)
    ].copy()

# -------------------------
# Subset by sample
# -------------------------
if args.sample is not None:
    adata = adata[
        adata.obs["sample"].isin(selected_samples)
    ].copy()

if adata.n_obs == 0:
    raise ValueError("No cells remain after subsetting.")

print(f"Subset result: {adata.n_obs} cells, {adata.n_vars} genes")

# -------------------------
# Validate existing Harmony PCA
# -------------------------
if "X_pca_harmony" not in adata.obsm:
    raise ValueError("X_pca_harmony not found in adata.obsm.")

if adata.obsm["X_pca_harmony"].shape[0] != adata.n_obs:
    raise ValueError("X_pca_harmony does not match the subset cells.")

# -------------------------
# Recompute neighbors on the subset
# -------------------------
sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=15
)

# -------------------------
# Leiden clustering
# Save cluster labels as 'leiden'
# -------------------------
sc.tl.leiden(
    adata,
    resolution=args.resolution,
    key_added="leiden"
)

print("\nLeiden cluster counts:")
print(adata.obs["leiden"].value_counts().sort_index())

# -------------------------
# UMAP coloured by Leiden clusters
# Do not save the UMAP embedding to the output
# -------------------------
sc.tl.umap(adata, random_state=42)

os.makedirs("figures", exist_ok=True)

sc.pl.umap(
    adata,
    color="leiden",
    size=20,
    show=False
)

plt.savefig(
    f"figures/{args.prefix}_clusters_umap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close("all")

# -------------------------
# Save subset
# -------------------------
adata.write_h5ad(args.output)
print(f"\nSaved: {args.output}")

# -------------------------
# Summary
# -------------------------
print("\nSummary:")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")
print(f"  Samples: {adata.obs['sample'].value_counts().to_dict()}")
print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")
print(f"  Resolution: {args.resolution}")
