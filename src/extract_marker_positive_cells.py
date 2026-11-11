#!/usr/bin/env python3

import argparse
import sys

import anndata as ad
import numpy as np
from scipy.sparse import issparse


def main():

    parser = argparse.ArgumentParser(
        description="Extract marker-positive and marker-negative cells."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--cutoff", type=float, default=3)
    parser.add_argument("--layer", default="log1p")
    parser.add_argument("--prefix", required=True)

    args = parser.parse_args()

    # Load AnnData
    print(f"Loading: {args.input}", flush=True)
    adata = ad.read_h5ad(args.input)

    # Check cell-type column
    if "celltype" not in adata.obs.columns:
        sys.exit("Error: 'celltype' not found in adata.obs.")

    # Check marker gene
    if args.marker not in adata.var_names:
        sys.exit(
            f"Error: Marker '{args.marker}' not found in adata.var_names."
        )

    # Check expression layer
    if args.layer not in adata.layers:
        sys.exit(
            f"Error: Layer '{args.layer}' not found in adata.layers."
        )

    # Subset requested cell type
    celltype_mask = (
        adata.obs["celltype"].astype(str) == args.celltype
    )

    if celltype_mask.sum() == 0:
        sys.exit(
            f"Error: No cells found for cell type '{args.celltype}'."
        )

    subset = adata[celltype_mask].copy()

    # Extract marker expression
    gene_index = subset.var_names.get_loc(args.marker)

    expression = subset.layers[args.layer][:, gene_index]

    if issparse(expression):
        expression = expression.toarray()

    expression = np.asarray(expression, dtype=np.float32).ravel()

    # Classify cells using the cutoff
    # FIX: Dynamically name the column based on the marker
    status_col = f"{args.marker}_status"
    
    subset.obs[status_col] = np.where(
        expression > args.cutoff,
        "positive",
        "negative"
    )

    # Count cells
    positive_count = (
        subset.obs[status_col] == "positive"
    ).sum()

    negative_count = (
        subset.obs[status_col] == "negative"
    ).sum()

    # Save combined AnnData object
    output = f"{args.prefix}.h5ad"

    subset.write_h5ad(output)

    # Report results
    print("\nExtraction complete.")
    print(f"Cell type: {args.celltype}")
    print(f"Marker: {args.marker}")
    print(f"Layer: {args.layer}")
    print(f"Cutoff: {args.cutoff}")
    # FIX: Dynamically print the counts using the marker name
    print(f"{args.marker}-positive cells: {positive_count:,}")
    print(f"{args.marker}-negative cells: {negative_count:,}")
    print(f"Total cells: {subset.n_obs:,}")
    print(f"Output: {output}", flush=True)


if __name__ == "__main__":
    main()
