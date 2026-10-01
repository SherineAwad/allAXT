import scanpy as sc
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", required=True)
args = parser.parse_args()

adata = sc.read(args.input)

if "connectivities" not in adata.obsp:
    raise ValueError("Neighbor graph missing. Run Harmony/neighbors first.")

# Leiden clustering
sc.tl.leiden(adata, resolution=2.5)

# UMAP (Leiden only)
sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    save=f"_{args.prefix}_leiden.png"
)


qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

for qc in qc_metrics:
    if qc in adata.obs:
        # 1. Create custom figure with desired width (e.g., width=16, height=6)
        fig, ax = plt.subplots(figsize=(16, 6))

        # 2. Pass ax=ax into scanpy's violin plot
        sc.pl.violin(
            adata,
            keys=qc,
            groupby="leiden",
            rotation=90,
            stripplot=False,
            jitter=0.2,
            order=sorted(adata.obs["leiden"].unique()),
            ax=ax,
            show=False,
        )

        # 3. Save directly using matplotlib to preserve the layout
        fig.savefig(
            f"figures/{args.prefix}_QC_{qc}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)

adata.write(args.output, compression="gzip")
