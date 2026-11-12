#!/usr/bin/env python3

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Infer CNVs from scRNA-seq using infercnvpy."
    )
    parser.add_argument("--input", required=True, help="Input AnnData h5ad file")
    parser.add_argument("--prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--marker",
        required=False,
        default=None,
        help="Optional marker gene to subset non-Uninjured cells (e.g., EGFP)"
    )
    parser.add_argument(
        "--cutoff",
        required=False,
        type=float,
        default=0.0,
        help="Expression cutoff for the marker to consider a cell positive (default: 0.0)"
    )
    return parser.parse_args()


def retrieve_mouse_coordinates():
    """
    Retrieve mouse genomic coordinates from Ensembl BioMart.
    Coordinates are matched using MGI gene symbols.
    """

    print("Retrieving mouse genomic coordinates from Ensembl BioMart...")

    annotations = sc.queries.biomart_annotations(
        "mmusculus",
        attrs=[
            "mgi_symbol",
            "chromosome_name",
            "start_position",
            "end_position",
        ],
        use_cache=True,
    )

    required = [
        "mgi_symbol",
        "chromosome_name",
        "start_position",
        "end_position",
    ]

    missing = set(required) - set(annotations.columns)

    if missing:
        raise ValueError(
            f"BioMart annotation is missing columns: {sorted(missing)}"
        )

    annotations = annotations[required].copy()

    annotations = annotations.dropna(subset=required)

    annotations["mgi_symbol"] = annotations["mgi_symbol"].astype(str).str.strip()
    annotations["chromosome_name"] = (
        annotations["chromosome_name"].astype(str).str.strip()
    )

    annotations["start_position"] = pd.to_numeric(
        annotations["start_position"], errors="coerce"
    )

    annotations["end_position"] = pd.to_numeric(
        annotations["end_position"], errors="coerce"
    )

    annotations = annotations.dropna(
        subset=["start_position", "end_position"]
    )

    annotations = annotations[
        (annotations["mgi_symbol"] != "")
        & (annotations["mgi_symbol"] != "nan")
    ].copy()

    # Retain standard mouse chromosomes only.
    valid_chromosomes = {
        str(i) for i in range(1, 20)
    } | {"X", "Y", "MT"}

    annotations = annotations[
        annotations["chromosome_name"].isin(valid_chromosomes)
    ].copy()

    annotations["start_position"] = (
        annotations["start_position"].astype(np.int64)
    )

    annotations["end_position"] = (
        annotations["end_position"].astype(np.int64)
    )

    # Remove exact duplicate gene-coordinate records.
    annotations = annotations.drop_duplicates(
        subset=[
            "mgi_symbol",
            "chromosome_name",
            "start_position",
            "end_position",
        ]
    )

    # A symbol mapping to multiple distinct genomic locations is ambiguous.
    loci_per_gene = (
        annotations.groupby("mgi_symbol")[
            ["chromosome_name", "start_position", "end_position"]
        ]
        .apply(lambda x: len(x.drop_duplicates()))
    )

    ambiguous_genes = loci_per_gene[loci_per_gene > 1].index

    annotations = annotations[
        ~annotations["mgi_symbol"].isin(ambiguous_genes)
    ].copy()

    # There should now be one genomic location per gene symbol.
    annotations = annotations.drop_duplicates(
        subset=["mgi_symbol"]
    )

    # CRITICAL FIX:
    # infercnvpy requires chromosome names beginning with "chr".
    # BioMart returns names such as 1, 2, 3, X, Y and MT.
    annotations["chromosome"] = (
        "chr" + annotations["chromosome_name"].replace({"MT": "M"})
    )

    annotations = annotations.rename(
        columns={
            "start_position": "start",
            "end_position": "end",
        }
    )

    annotations = annotations[
        ["mgi_symbol", "chromosome", "start", "end"]
    ].copy()

    print(
        f"Unique BioMart gene symbols: {annotations['mgi_symbol'].nunique():,}"
    )
    print(
        f"Ambiguous gene symbols excluded: {len(ambiguous_genes):,}"
    )

    return annotations


def annotate_genes(adata, prefix):
    """
    Match AnnData gene symbols to BioMart coordinates.
    """

    annotations = retrieve_mouse_coordinates()

    gene_coordinates = annotations.set_index("mgi_symbol")

    gene_symbols = pd.Index(adata.var_names.astype(str))

    mapped = gene_coordinates.reindex(gene_symbols)

    valid = mapped["chromosome"].notna().to_numpy()

    unmapped_genes = gene_symbols[~valid]

    if len(unmapped_genes) > 0:
        unmapped_df = pd.DataFrame({
            "gene": unmapped_genes
        })

        output = f"{prefix}_CNV_Unmapped_Genes.csv"

        unmapped_df.to_csv(output, index=False)

        print(f"Genes without unique coordinates: {len(unmapped_genes):,}")
        print(f"Saved: {output}")

    adata = adata[:, valid].copy()

    mapped = mapped.loc[valid].copy()

    # Assign the genomic coordinates to adata.var.
    adata.var["chromosome"] = mapped["chromosome"].to_numpy()
    adata.var["start"] = mapped["start"].to_numpy(dtype=np.int64)
    adata.var["end"] = mapped["end"].to_numpy(dtype=np.int64)

    # Exclude mitochondrial genes (chrM / chrMT) from the analysis dataset
    mito_genes = adata.var["chromosome"].isin(["chrM", "chrMT"])
    if mito_genes.any():
        adata = adata[:, ~mito_genes].copy()

    # Confirm that chromosome names meet infercnvpy requirements.
    chromosomes = adata.var["chromosome"].dropna().unique()

    if len(chromosomes) == 0:
        raise ValueError("No chromosomes with valid genomic coordinates.")

    invalid = [
        chrom for chrom in chromosomes
        if not str(chrom).startswith("chr")
    ]

    if invalid:
        raise ValueError(
            f"Invalid chromosome names detected: {invalid[:10]}"
        )

    usable_chromosomes = [
        chrom for chrom in chromosomes
        if chrom != "chrM"
    ]

    if not usable_chromosomes:
        raise ValueError(
            "No usable chromosomes remain after excluding chrM."
        )

    print(f"Genes retained: {adata.n_vars:,}")
    print(f"Usable chromosomes: {len(usable_chromosomes)}")
    print(
        "Chromosomes:",
        ", ".join(sorted(usable_chromosomes))
    )

    # Check that each chromosome has enough genes for the default window.
    chromosome_counts = adata.var["chromosome"].value_counts()

    too_small = chromosome_counts[
        (chromosome_counts.index != "chrM")
        & (chromosome_counts < 2)
    ]

    if len(too_small):
        warnings.warn(
            "Some chromosomes have fewer than two genes: "
            + ", ".join(
                f"{chrom} ({count})"
                for chrom, count in too_small.items()
            )
        )

    return adata


def filter_by_marker(adata, marker=None, cutoff=0.0):
    """
    Filter non-Uninjured cells based on marker expression.
    Keeps all Uninjured cells.
    """
    if marker is None:
        return adata

    print(f"Filtering non-Uninjured cells for {marker} > {cutoff}...")

    if "sample" not in adata.obs.columns:
        raise ValueError("Missing required adata.obs column: sample")

    uninjured = adata.obs["sample"].astype(str) == "Uninjured"

    # Check for marker in var_names
    if marker in adata.var_names:
        print(f"  Found marker '{marker}' in var_names.")
        marker_expr = adata[:, marker].X
        if sparse.issparse(marker_expr):
            marker_expr = marker_expr.toarray().flatten()
        else:
            marker_expr = np.array(marker_expr).flatten()
        marker_positive = marker_expr > cutoff
        
    # Check for marker in obs
    elif marker in adata.obs.columns:
        print(f"  Found marker '{marker}' in obs columns.")
        marker_positive = adata.obs[marker].to_numpy() > cutoff
        
    else:
        raise ValueError(
            f"Marker '{marker}' not found in adata.var_names or adata.obs. "
            f"Please check the exact spelling."
        )

    # Keep all Uninjured cells + marker positive cells from other samples
    keep_mask = uninjured | marker_positive
    adata = adata[keep_mask].copy()
    
    print(f"Cells kept after marker filtering: {adata.n_obs:,}")
    
    return adata


def prepare_reference(adata):
    """
    Use ALL Uninjured cells as a single, unified reference population.
    """

    print("Preparing Uninjured reference populations...")

    if "sample" not in adata.obs.columns:
        raise ValueError("Missing required adata.obs column: sample")

    uninjured = adata.obs["sample"].astype(str) == "Uninjured"

    if not uninjured.any():
        raise ValueError("No Uninjured cells found.")

    adata.obs["_cnv_reference_group"] = "non_reference"

    # Group all Uninjured cells into one single category
    adata.obs.loc[uninjured, "_cnv_reference_group"] = "Reference"
    
    reference_categories = ["Reference"]

    print("Reference populations:")
    print(f"  Reference: {int(uninjured.sum()):,} cells")
    print(f"Total reference cells: {int(uninjured.sum()):,}")

    return adata, reference_categories


def run_infercnv(adata, reference_categories):
    """
    Run infercnvpy using Uninjured reference populations.
    """

    print("Running infercnvpy...")

    cnv.tl.infercnv(
        adata,
        reference_key="_cnv_reference_group",
        reference_cat=reference_categories,
        key_added="cnv",
        n_jobs=4,
        chunksize=1000,
    )

    print("infercnvpy completed successfully.")

    return adata


def save_outputs(adata, prefix):
    """
    Save CNV results, cell-type summary CSV, chromosome heatmap,
    and CNV score violin plot.
    """

    os.makedirs("figures", exist_ok=True)

    # ---------------------------------------------------------
    # Save AnnData
    # ---------------------------------------------------------

    h5ad_output = f"{prefix}_CNV.h5ad"

    print(f"Saving: {h5ad_output}")

    adata.write_h5ad(h5ad_output)

    # ---------------------------------------------------------
    # Save average CNV signal per cell type
    # ---------------------------------------------------------

    print("Calculating mean CNV signal per cell type...")

    cnv_matrix = adata.obsm["X_cnv"]

    if sparse.issparse(cnv_matrix):
        cnv_matrix = cnv_matrix.toarray()

    celltypes = adata.obs["celltype"].astype(str)

    summary = pd.DataFrame(
        cnv_matrix,
        index=adata.obs_names
    )

    summary["celltype"] = celltypes.to_numpy()

    mean_cnv = summary.groupby("celltype", observed=True).mean()

    mean_cnv.columns = [
        f"region_{i + 1:05d}"
        for i in range(mean_cnv.shape[1])
    ]

    csv_output = f"{prefix}_CNV_Summary.csv"

    mean_cnv.to_csv(csv_output)

    print(f"Saved: {csv_output}")

    # ---------------------------------------------------------
    # Save chromosome heatmap
    # ---------------------------------------------------------

    print("Generating chromosome heatmap...")

    plt.figure(figsize=(16, 8))

    cnv.pl.chromosome_heatmap_summary(
        adata,
        groupby="celltype",
        use_rep="cnv",
        figsize=(16, 8),
        show=False,
        dendrogram=False,
    )

    heatmap_output = f"figures/{prefix}_CNV_Heatmap.png"

    plt.savefig(
        heatmap_output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close("all")

    print(f"Saved: {heatmap_output}")

    # ---------------------------------------------------------
    # CNV score per cell and violin plot by cell type
    # ---------------------------------------------------------

    print("Calculating per-cell CNV score...")

    if sparse.issparse(cnv_matrix):
        cnv_matrix = cnv_matrix.toarray()

    # Per-cell CNV burden: mean absolute CNV signal across bins.
    # This is a robust overall "how abnormal is this cell" score.
    per_cell_score = np.abs(cnv_matrix).mean(axis=1)

    adata.obs["cnv_score"] = per_cell_score

    print("Generating CNV score violin plot...")

    plt.figure(figsize=(12, 6))

    sc.pl.violin(
        adata,
        keys="cnv_score",
        groupby="celltype",
        rotation=90,
        show=False,
    )

    violin_output = f"figures/{prefix}_CNV_Score_Violin.png"

    plt.savefig(
        violin_output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close("all")

    print(f"Saved: {violin_output}")


def main():

    args = parse_args()

    print(f"Reading {args.input}")

    adata = sc.read_h5ad(args.input)

    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")

    if not adata.var_names.is_unique:
        raise ValueError(
            "Gene names are not unique. Resolve duplicate gene symbols first."
        )

    # 1. Filter by marker FIRST (before dropping unmapped genes like EGFP)
    adata = filter_by_marker(adata, marker=args.marker, cutoff=args.cutoff)

    # 2. Gene symbols must be matched to genomic coordinates.
    adata = annotate_genes(adata, args.prefix)

    # 3. Prepare Uninjured reference populations.
    adata, reference_categories = prepare_reference(adata)

    # 4. Run infercnvpy.
    adata = run_infercnv(adata, reference_categories)

    # 5. Save outputs.
    save_outputs(adata, args.prefix)

    print("\nCNV analysis completed.")
    print(f"Output prefix: {args.prefix}")


if __name__ == "__main__":
    main()
