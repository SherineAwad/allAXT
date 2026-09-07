import scanpy as sc
import sys
import importlib_metadata
import argparse
import numpy as np
import scipy.sparse as sp

sys.modules['importlib.metadata'] = importlib_metadata


def read_samples(file_path):
    samples = {}
    with open(file_path, 'r') as file:
        for line in file:
            sample_id = line.strip()
            if sample_id:
                filename = f"{sample_id}_filtered_feature_bc_matrix.h5"
                samples[sample_id] = filename
    return samples


def main():
    parser = argparse.ArgumentParser()

    # Required args
    parser.add_argument('--input', required=True, help="Samples file")
    parser.add_argument('--output', required=True, help="Output .h5ad file")
    parser.add_argument('--prefix', required=True, help="Prefix for saved figures")

    # QC parameters 
    parser.add_argument('--min_genes', type=int, default=500)
    parser.add_argument('--max_genes', type=int, default=10000)
    parser.add_argument('--min_counts', type=int, default=1000)
    parser.add_argument('--max_counts', type=int, default=50000)
    parser.add_argument('--max_mt', type=float, default=10.0)

    parser.add_argument('--min_cells_gene', type=int, default=3)
    parser.add_argument('--min_genes_cell', type=int, default=500)

    args = parser.parse_args()

    adatas = {}
    samples = read_samples(args.input)
    print(samples)

    for sample_id, filename in samples.items():
        print(f"Reading {filename}...")
        adata = sc.read_10x_h5(filename)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_id
        adatas[sample_id] = adata

    combined_adata = sc.concat(adatas.values(), label='sample', keys=adatas.keys())

    combined_adata.var["mt"] = combined_adata.var_names.str.startswith("mt-")

    sc.pp.calculate_qc_metrics(
        combined_adata, qc_vars=["mt"], inplace=True, log1p=True
    )

    # BEFORE QC plot
    sc.pl.violin(
        combined_adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{args.prefix}_preQC.png"
    )

    # -----------------------------
    # CELL FILTERING (your logic)
    # -----------------------------
    combined_adata = combined_adata[
        (combined_adata.obs['n_genes_by_counts'] > args.min_genes) &
        (combined_adata.obs['n_genes_by_counts'] < args.max_genes) &
        (combined_adata.obs['total_counts'] > args.min_counts) &
        (combined_adata.obs['total_counts'] < args.max_counts) &
        (combined_adata.obs['pct_counts_mt'] < args.max_mt), :
    ]

    sc.pp.filter_cells(combined_adata, min_genes=args.min_genes_cell)
    sc.pp.filter_genes(combined_adata, min_cells=args.min_cells_gene)


    # AFTER QC plot
    sc.pl.violin(
        combined_adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{args.prefix}_AfterQC.png"
    )

    combined_adata.write(args.output)


if __name__ == "__main__":
    main()
