import scanpy as sc
import argparse
import os
import numpy as np
import harmonypy as hm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Run Harmony – input must be (cells × PCs)
    harmony_out = hm.run_harmony(adata.obsm['X_pca'], adata.obs, "experiment")
    adata.obsm['X_pca_harmony'] = harmony_out.Z_corr  # shape (n_cells, n_pcs)

    # UMAP using Harmony-corrected PCs
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)

    # Plot
    os.makedirs('figures', exist_ok=True)
    sc.pl.umap(adata, color='sample', save=f'_{args.prefix}_harmony_sample.png')
    

    # Each sample separately
    for s in adata.obs["sample"].unique():
        sc.pl.umap(
            adata[adata.obs["sample"] == s],
            color='sample',
            title=f"Sample: {s}",
            save=f'_{args.prefix}_harmony_{s}.png')

    # Save
    adata.write(args.output, compression='gzip')
    print(f'SAVED → {args.output}')

if __name__ == '__main__':
    main()
