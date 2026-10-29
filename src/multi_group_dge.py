import argparse
import pandas as pd
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run pairwise DGE using log1p layer and save individual per-comparison CSVs."
    )
    parser.add_argument("--input", required=True, help="Path to input h5ad file")
    parser.add_argument(
        "--comparisons",
        required=True,
        help="Path to file with 5 columns: comp_prefix, target_celltype, target_condition, ref_celltype, ref_condition",
    )
    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for significance (default: 0.05)",
    )
    parser.add_argument(
        "--log2fc",
        type=float,
        default=0.25,
        help="Absolute log2 fold change cutoff for significance (default: 0.25)",
    )
    parser.add_argument(
        "--layer",
        type=str,
        default="log1p",
        help="Layer to use for differential expression testing (default: log1p)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # 1. Load dataset
    print(f"Loading AnnData from {args.input}...")
    adata = sc.read_h5ad(args.input)

    # Validate metadata
    if "celltype" not in adata.obs.columns or "sample" not in adata.obs.columns:
        raise ValueError("adata.obs must contain 'celltype' and 'sample' columns.")

    adata.obs["celltype_sample"] = (
        adata.obs["celltype"].astype(str) + "_" + adata.obs["sample"].astype(str)
    )

    # Validate layer selection
    use_layer = None
    if args.layer in adata.layers:
        use_layer = args.layer
        print(f"Using specified expression layer: '{use_layer}'")
    else:
        print(
            f"WARNING: Layer '{args.layer}' not found in adata.layers. Defaulting to main matrix `adata.X`."
        )

    # 2. Parse comparisons file
    sep = "," if args.comparisons.endswith(".csv") else "\t"
    comp_df = pd.read_csv(args.comparisons, sep=sep, header=None).dropna(how="all")

    # Drop header row if present
    first_cell = str(comp_df.iloc[0, 0]).strip().lower()
    if first_cell in ["comp_prefix", "prefix", "comparison_prefix", "name"]:
        comp_df = comp_df.iloc[1:].reset_index(drop=True)

    if comp_df.shape[1] < 5:
        raise ValueError("The comparisons file must have at least 5 columns.")

    comp_df = comp_df.iloc[:, :5]
    comp_df.columns = [
        "comp_prefix",
        "target_celltype",
        "target_sample",
        "ref_celltype",
        "ref_sample",
    ]

    available_groups = set(adata.obs["celltype_sample"].unique())

    # 3. Iterate through comparison tests and save individual files only
    for idx, row in comp_df.iterrows():
        comp_prefix = str(row["comp_prefix"]).strip()
        target_celltype = str(row["target_celltype"]).strip()
        target_sample = str(row["target_sample"]).strip()
        ref_celltype = str(row["ref_celltype"]).strip()
        ref_sample = str(row["ref_sample"]).strip()

        target_group = f"{target_celltype}_{target_sample}"
        ref_group = f"{ref_celltype}_{ref_sample}"

        print(
            f"\n[{idx + 1}/{len(comp_df)}] Running comparison: {comp_prefix} ({target_group} vs {ref_group})"
        )

        if target_group not in available_groups or ref_group not in available_groups:
            print(
                f"  WARNING: Skipping {comp_prefix}. One or both groups missing from dataset."
            )
            print(
                f"  Target in dataset: {target_group in available_groups} | Ref in dataset: {ref_group in available_groups}"
            )
            continue

        # Subset for pairwise test
        mask = adata.obs["celltype_sample"].isin([target_group, ref_group])
        sub_adata = adata[mask].copy()

        # Execute Differential Expression on log1p layer
        sc.tl.rank_genes_groups(
            sub_adata,
            groupby="celltype_sample",
            groups=[target_group],
            reference=ref_group,
            method="wilcoxon",
            layer=use_layer,
            use_raw=False,
        )

        # Structure DGE dataframe
        full_result = sc.get.rank_genes_groups_df(sub_adata, group=target_group)
        full_result["comparison"] = comp_prefix
        full_result["target_group"] = target_group
        full_result["ref_group"] = ref_group

        # Filter significant genes using CLI thresholds
        sig_result = full_result[
            (full_result["pvals_adj"] < args.padj)
            & (full_result["logfoldchanges"].abs() >= args.log2fc)
        ].copy()

        # Save individual Full and Significant CSVs ONLY
        full_csv = f"{comp_prefix}_DGE_full.csv"
        sig_csv = f"{comp_prefix}_DGE_sig.csv"

        full_result.to_csv(full_csv, index=False)
        sig_result.to_csv(sig_csv, index=False)

        print(f"  Saved full DGE ({len(full_result)} genes) to: {full_csv}")
        print(f"  Saved sig DGE ({len(sig_result)} genes) to: {sig_csv}")


if __name__ == "__main__":
    main()
