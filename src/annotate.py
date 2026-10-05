#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--annotations", required=True)
parser.add_argument("--remove", required=False, help="Comma-separated list of cluster numbers to remove (e.g., '5,7,9')")
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

# Remove specified clusters if provided
if args.remove:
    remove_clusters = [c.strip() for c in args.remove.split(",")]
    original_cells = adata.n_obs
    adata = adata[~adata.obs["leiden"].astype(str).isin(remove_clusters)]
    removed_cells = original_cells - adata.n_obs
    print(f"Removed clusters: {remove_clusters}")
    print(f"Removed {removed_cells} cells, {adata.n_obs} remaining")

annotations = pd.read_csv(args.annotations, header=None, names=["leiden", "celltype"])
annotations["leiden"] = annotations["leiden"].astype(str)

mapping = dict(zip(annotations["leiden"], annotations["celltype"]))

present_clusters = set(adata.obs["leiden"].astype(str).unique())
mapped_clusters = set(mapping.keys())

for cluster in present_clusters:
    if cluster not in mapped_clusters:
        print(f"WARNING: Cluster {cluster} has no annotation - will be labeled as 'unknown'")

adata.obs["celltype"] = adata.obs["leiden"].astype(str).map(mapping).fillna("unknown")

adata.write_h5ad(args.output)

# UMAP with annotation written on the clusters
sc.pl.umap(adata, color="celltype", legend_loc="on data", title="Cell Type Annotation", save=f"_{args.prefix}_celltypeON.png")

# UMAP with annotation in legend only
sc.pl.umap(adata, color="celltype", legend_loc="right margin", title="Cell Type Annotation", save=f"_{args.prefix}_celltype.png")
