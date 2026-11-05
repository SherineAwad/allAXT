import argparse
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pertpy as pt

from pathlib import Path
from scipy.sparse import issparse


def run_global_comparisons(
    h5ad_path,
    comparisons_file,
    prefix,
    condition_col="sample",
    celltype_col="celltype",
    n_neighbors=35,
    n_nhoods=10000
):

    # --------------------------------------------------
    # 1. Load AnnData
    # --------------------------------------------------

    print(f"1. Loading global AnnData object from {h5ad_path}...")

    adata = ad.read_h5ad(h5ad_path)

    print(f"    Loaded {adata.n_obs} cells.")

    for col in [condition_col, celltype_col]:
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not found in adata.obs."
            )

        if adata.obs[col].isna().any():
            raise ValueError(
                f"Missing values found in adata.obs['{col}']."
            )

    # --------------------------------------------------
    # 2. Global neighbour graph
    # --------------------------------------------------

    print("2. Ensuring global neighbour graph...")

    if (
        "neighbors" not in adata.uns
        or "connectivities" not in adata.obsp
    ):

        if "X_pca" not in adata.obsm:
            raise ValueError(
                "X_pca is missing from adata.obsm. "
                "Cannot construct the global neighbour graph."
            )

        sc.pp.neighbors(
            adata,
            n_neighbors=min(n_neighbors, adata.n_obs - 1),
            use_rep="X_pca"
        )

    connectivities = adata.obsp["connectivities"]

    if connectivities.shape != (adata.n_obs, adata.n_obs):
        raise ValueError(
            f"Unexpected connectivity matrix shape: "
            f"{connectivities.shape}"
        )

    print(f"    Connectivity matrix: {connectivities.shape}")

    # --------------------------------------------------
    # 3. Construct global Milo neighbourhoods
    # --------------------------------------------------

    print("3. Defining global Milo neighbourhoods...")

    milo = pt.tl.Milo()

    prop = min(n_nhoods / adata.n_obs, 1.0)

    milo.make_nhoods(
        adata,
        prop=prop,
        seed=42
    )

    if "nhoods" not in adata.obsm:
        raise RuntimeError(
            "make_nhoods did not populate adata.obsm['nhoods']."
        )

    nhood_matrix = adata.obsm["nhoods"]

    # Ensure rows represent cells and columns represent neighbourhoods.
    if nhood_matrix.shape[0] != adata.n_obs:

        if nhood_matrix.shape[1] == adata.n_obs:
            nhood_matrix = nhood_matrix.T

        else:
            raise RuntimeError(
                f"Unexpected neighbourhood matrix shape: "
                f"{nhood_matrix.shape}"
            )

    print(
        f"    Global neighbourhoods: {nhood_matrix.shape[1]}"
    )

    # --------------------------------------------------
    # 4. Read comparisons
    # --------------------------------------------------

    print(f"4. Reading comparisons file: {comparisons_file}...")

    df_comp = pd.read_csv(
        comparisons_file,
        sep=r"\s+",
        header=None,
        comment="#"
    )

    if df_comp.shape[1] != 4:
        raise ValueError(
            "The comparisons file must contain exactly four columns: "
            "celltype1 condition1 celltype2 condition2."
        )

    if df_comp.empty:
        raise ValueError(
            "The comparisons file contains no comparisons."
        )

    # --------------------------------------------------
    # 5. Prepare metadata
    # --------------------------------------------------

    celltypes = adata.obs[celltype_col].astype(str)
    conditions = adata.obs[condition_col].astype(str)

    available_celltypes = set(celltypes.unique())
    available_conditions = set(conditions.unique())

    print("    Available cell types:", sorted(available_celltypes))
    print("    Available conditions:", sorted(available_conditions))

    # --------------------------------------------------
    # 6. Output directory
    # --------------------------------------------------

    output_dir = Path("figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # 7. Process comparisons
    # --------------------------------------------------

    print(
        f"5. Processing {len(df_comp)} comparisons "
        "on the global manifold..."
    )

    for idx, row in df_comp.iterrows():

        ct1, cond1 = str(row.iloc[0]), str(row.iloc[1])
        ct2, cond2 = str(row.iloc[2]), str(row.iloc[3])

        print("\n" + "=" * 65)

        print(
            f"Comparison {idx + 1}: "
            f"{cond2} ({ct2}) vs {cond1} ({ct1})"
        )

        # --------------------------------------------------
        # 7a. Validate comparison
        # --------------------------------------------------

        if ct1 not in available_celltypes:
            print(f"    WARNING: Cell type '{ct1}' not found. Skipping.")
            continue

        if ct2 not in available_celltypes:
            print(f"    WARNING: Cell type '{ct2}' not found. Skipping.")
            continue

        if cond1 not in available_conditions:
            print(f"    WARNING: Condition '{cond1}' not found. Skipping.")
            continue

        if cond2 not in available_conditions:
            print(f"    WARNING: Condition '{cond2}' not found. Skipping.")
            continue

        # --------------------------------------------------
        # 7b. Population masks
        # --------------------------------------------------

        mask1 = (
            (celltypes == ct1) &
            (conditions == cond1)
        ).to_numpy(dtype=float)

        mask2 = (
            (celltypes == ct2) &
            (conditions == cond2)
        ).to_numpy(dtype=float)

        n1 = int(mask1.sum())
        n2 = int(mask2.sum())

        print(f"    Cells in first population: {n1}")
        print(f"    Cells in second population: {n2}")

        if n1 == 0 or n2 == 0:
            print(
                "    WARNING: One or both populations are empty. "
                "Skipping comparison."
            )
            continue

        # --------------------------------------------------
        # 7c. Global neighbourhood abundance ratios
        # --------------------------------------------------

        if issparse(nhood_matrix):

            counts1 = np.asarray(
                nhood_matrix.T @ mask1
            ).ravel()

            counts2 = np.asarray(
                nhood_matrix.T @ mask2
            ).ravel()

        else:

            counts1 = np.asarray(
                nhood_matrix.T @ mask1
            ).ravel()

            counts2 = np.asarray(
                nhood_matrix.T @ mask2
            ).ravel()

        log2_ratios = np.log2(
            (counts2 + 1) / (counts1 + 1)
        )

        # --------------------------------------------------
        # 7d. Output filename
        # --------------------------------------------------

        safe_name = (
            f"{prefix}_{cond2}_{ct2}_vs_{cond1}_{ct1}"
            .replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("/", "_")
        )

        # --------------------------------------------------
        # 7e. Histogram
        # --------------------------------------------------

        plt.figure(figsize=(10, 6), dpi=300)

        plt.hist(
            log2_ratios,
            bins=30,
            edgecolor="black",
            linewidth=0.5
        )

        plt.axvline(
            0,
            linestyle="--",
            color="tab:blue"
        )

        plt.xlabel(
            f"log2((count + 1) in {cond2} [{ct2}] / "
            f"(count + 1) in {cond1} [{ct1}])",
            fontsize=12
        )

        plt.ylabel(
            "Number of neighbourhoods",
            fontsize=12
        )

        plt.title(
            "Global neighbourhood abundance ratio",
            fontsize=14
        )

        plt.tight_layout()

        hist_path = (
            output_dir /
            f"{safe_name}_ratio_histogram.png"
        )

        plt.savefig(hist_path)
        plt.close()

        print(f"    Saved histogram: {hist_path}")

        # --------------------------------------------------
        # 7f. Save neighbourhood ratios
        # --------------------------------------------------

        ratio_results = pd.DataFrame({
            "neighborhood": np.arange(len(log2_ratios)),
            f"count_{cond1}_{ct1}": counts1,
            f"count_{cond2}_{ct2}": counts2,
            "log2_ratio": log2_ratios
        })

        ratio_path = (
            output_dir /
            f"{safe_name}_neighborhood_ratios.csv"
        )

        ratio_results.to_csv(
            ratio_path,
            index=False
        )

        print(f"    Saved neighbourhood ratios: {ratio_path}")

    # --------------------------------------------------
    # 8. Differential abundance limitation
    # --------------------------------------------------

    print("\n" + "=" * 65)
    print("Milo differential abundance testing was not performed.")
    print(
        "The supplied analysis does not identify independent "
        "biological replicates for statistical testing."
    )
    print(
        "Global neighbourhood abundance ratios and their "
        "descriptive outputs were generated for all valid comparisons."
    )

    print("\nAnalysis completed.")


# --------------------------------------------------
# Main
# --------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Global Milo neighbourhood abundance comparisons."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to global AnnData h5ad file"
    )

    parser.add_argument(
        "--comparisons",
        required=True,
        help="Path to comparisons text file"
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files"
    )

    parser.add_argument(
        "--n-cells",
        type=int,
        default=10000,
        help="Number of neighbourhood centres (default: 10000)"
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=35,
        help="Number of neighbours (default: 35)"
    )

    args = parser.parse_args()

    run_global_comparisons(
        h5ad_path=args.input,
        comparisons_file=args.comparisons,
        prefix=args.prefix,
        n_nhoods=args.n_cells,
        n_neighbors=args.n_neighbors
    )
