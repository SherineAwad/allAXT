#!/usr/bin/env python3

import argparse
import os
import textwrap
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import gseapy as gp

warnings.filterwarnings("ignore")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run multi-comparison GSEA using GSEApy with case-insensitive human contamination filtering, numeric ranking sanitization, duplicate gene resolution, and gene identifier overlap verification."
    )
    parser.add_argument(
        "--comparisons",
        required=True,
        help="Path to file specifying comparisons (CSV or whitespace-delimited table: dge_csv, comp_prefix)",
    )
    parser.add_argument(
        "--organism",
        default="mouse",
        help="Organism name (default: mouse)",
    )
    parser.add_argument(
        "--library",
        default="KEGG_2019_Mouse",
        help="Enrichr library name (default: KEGG_2019_Mouse)",
    )
    parser.add_argument(
        "--resources",
        default=None,
        help="Resource/library name shortcut (e.g. kegg, go, bp, cc, mf, reactome)",
    )
    parser.add_argument(
        "--prefix",
        default="multi_gsea",
        help="Prefix for final output files (default: multi_gsea)",
    )
    parser.add_argument(
        "--qvalue",
        type=float,
        default=0.05,
        help="FDR q-value threshold for significance (default: 0.05)",
    )
    parser.add_argument(
        "--pvalue",
        type=float,
        default=0.05,
        help="Nominal p-value threshold for significance (default: 0.05)",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=5,
        help="Number of top significant enriched/depleted pathways per comparison for unified plotting",
    )
    return parser.parse_args()


def parse_gene_percent(value):
    """Convert GSEApy Gene % values to numeric percentages."""
    if pd.isna(value):
        return np.nan

    value = str(value).strip()

    try:
        if "/" in value:
            numerator, denominator = value.split("/", 1)
            denominator = float(denominator)
            return (
                np.nan
                if denominator == 0
                else float(numerator) / denominator * 100
            )

        if value.endswith("%"):
            return float(value.rstrip("%"))

        return float(value)

    except (ValueError, ZeroDivisionError):
        return np.nan


def resolve_library_name(lib_input, organism):
    """Map short resource names (kegg, go, bp, cc, mf, reactome) to valid Enrichr/GSEApy library names."""
    if not lib_input:
        return "KEGG_2019_Mouse"
    
    val = lib_input.strip().lower()
    org = organism.strip().lower()

    is_mouse = org in ["mouse", "mus musculus"]

    if val in ["kegg", "KEGG"]:
        return "KEGG_2019_Mouse" if is_mouse else "KEGG_2019"
    elif val in ["go", "GO", "gene_ontology", "bp"]:
        return "GO_Biological_Process_2021" if is_mouse else "GO_Biological_Process_2021"
    elif val in ["cc", "cellular_component"]:
        return "GO_Cellular_Component_2021" if is_mouse else "GO_Cellular_Component_2021"
    elif val in ["mf", "molecular_function"]:
        return "GO_Molecular_Function_2021" if is_mouse else "GO_Molecular_Function_2021"
    elif val in ["reactome", "Reactome"]:
        return "Reactome_2022" if is_mouse else "Reactome_2016"
    
    return lib_input


def load_filtered_gene_sets(library, organism):
    """Fetch library via GSEApy with explicit fallback logging and strictly purge human contamination."""
    resolved_lib = resolve_library_name(library, organism)
    print(f"Fetching library '{resolved_lib}' (from input '{library}') for organism '{organism}'...")
    
    loaded_library_name = resolved_lib
    try:
        raw_sets = gp.get_library(name=resolved_lib, organism=organism)
    except Exception as e:
        if organism.lower() in ["mouse", "mus musculus"]:
            try:
                alt_lib = resolved_lib.replace("_Mouse", "").replace("_mus_musculus", "")
                print(f"Retrying fetch with alternative library name '{alt_lib}'...")
                raw_sets = gp.get_library(name=alt_lib, organism=organism)
                loaded_library_name = alt_lib
            except Exception as e2:
                if "reactome" in str(library).lower():
                    try:
                        fallback_reactome = "Reactome_2016"
                        print(f"WARNING: Primary Reactome retrieval failed. Explicitly falling back to '{fallback_reactome}'...")
                        raw_sets = gp.get_library(name=fallback_reactome, organism=organism)
                        loaded_library_name = fallback_reactome
                    except Exception as e3:
                        raise RuntimeError(f"Failed to retrieve library '{library}' (resolved as '{resolved_lib}'): {e3}")
                else:
                    raise RuntimeError(f"Failed to retrieve library '{library}' (resolved as '{resolved_lib}'): {e2}")
        else:
            raise RuntimeError(f"Failed to retrieve library '{library}': {e}")

    print(f"Successfully loaded library version: '{loaded_library_name}'")

    if not isinstance(raw_sets, dict):
        raise RuntimeError(f"Failed to retrieve valid gene set dictionary from library '{library}'.")

    filtered_sets = {}
    excluded_count = 0

    for pathway, genes in raw_sets.items():
        pathway_lower = str(pathway).lower()
        
        # Strictly reject any pathway name containing human or hsa in any capitalization format
        if organism.lower() in ["mouse", "mus musculus"] and (
            "human" in pathway_lower or "hsa" in pathway_lower
        ):
            excluded_count += 1
            continue
            
        filtered_sets[pathway] = genes

    print(f"Loaded {len(filtered_sets)} clean pathways (excluded {excluded_count} human-contaminated sets).")
    return filtered_sets, loaded_library_name


def run_single_gsea(
    csv_file,
    comp_prefix,
    gene_sets,
    qval_cutoff,
    pval_cutoff,
    prefix,
):
    """Execute preranked GSEA for a single DGE comparison table."""
    print(f"\nProcessing comparison: {comp_prefix} ({csv_file})")

    df = pd.read_csv(csv_file)

    gene_col = "gene" if "gene" in df.columns else ("names" if "names" in df.columns else df.columns[0])
    score_col = "scores" if "scores" in df.columns else ("stat" if "stat" in df.columns else df.columns[1])

    # 1. Sanitize ranking scores numerically
    df["ranking"] = pd.to_numeric(df[score_col], errors="coerce")

    # 2. Clean gene identifiers safely using string datatype and filter out empty strings / NaNs
    df[gene_col] = df[gene_col].astype("string").str.strip()
    df = df.dropna(subset=[gene_col, "ranking"])
    df = df[df[gene_col].ne("")]

    # 3. Resolve duplicate gene names by keeping the entry with the strongest absolute score
    df["_abs_score"] = df["ranking"].abs()
    df = df.sort_values(by="_abs_score", ascending=False).drop_duplicates(subset=[gene_col], keep="first")
    df = df.drop(columns=["_abs_score"])

    ranking = df[[gene_col, "ranking"]].sort_values("ranking", ascending=False)
    ranking.columns = ["Gene", "Score"]

    # 4. Gene Identifier Match / Overlap Verification against Library
    all_library_genes = set(g for genes in gene_sets.values() for g in genes)
    input_genes = set(ranking["Gene"])
    overlapping_genes = input_genes.intersection(all_library_genes)
    overlap_pct = (len(overlapping_genes) / len(input_genes) * 100) if len(input_genes) > 0 else 0

    print(f"  Gene Identifier Check: {len(overlapping_genes)} / {len(input_genes)} input genes found in library ({overlap_pct:.2f}% overlap)")
    if overlap_pct < 15.0:
        print("  WARNING: Very low gene overlap detected! Verify that your input gene symbols match the mouse library naming conventions.")

    # Handle tied ranking scores cleanly
    dup_count = ranking["Score"].duplicated().sum()
    if dup_count > 0:
        np.random.seed(42)
        ranking["Score"] = ranking["Score"] + np.random.normal(0, 1e-10, len(ranking))

    out_dir = f"{comp_prefix}_gsea_output"

    results = gp.prerank(
        rnk=ranking,
        gene_sets=gene_sets,
        permutation_num=1000,
        outdir=out_dir,
        min_size=15,
        max_size=500,
        seed=42,
    )

    results_df = results.res2d

    numeric_cols = ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %"]
    for col in numeric_cols:
        if col in results_df.columns:
            results_df[col] = pd.to_numeric(results_df[col], errors="coerce")

    results_df.to_csv(f"{prefix}_{comp_prefix}_full.csv", index=False)

    sig_df = results_df[
        (results_df["FDR q-val"] < qval_cutoff)
        & (results_df["NOM p-val"] < pval_cutoff)
    ].copy()

    sig_df.to_csv(f"{prefix}_{comp_prefix}_significant.csv", index=False)

    sig_df["comparison"] = comp_prefix
    results_df["comparison"] = comp_prefix

    print(
        f"  Total pathways: {len(results_df)} | "
        f"Significant (FDR < {qval_cutoff}, p-val < {pval_cutoff}): {len(sig_df)}"
    )

    return results_df, sig_df


def main():
    args = parse_arguments()

    active_lib = args.resources if args.resources else args.library

    os.makedirs("figures", exist_ok=True)

    if args.comparisons.endswith(".csv"):
        comp_list_df = pd.read_csv(args.comparisons, header=None).dropna(how="all")
    else:
        comp_list_df = pd.read_csv(args.comparisons, sep=r"\s+", header=None).dropna(how="all")

    comp_list_df = comp_list_df.dropna(axis=1, how="all")
    if comp_list_df.shape[1] >= 2:
        comp_list_df = comp_list_df.iloc[:, :2]
        comp_list_df.columns = ["dge_csv", "comp_prefix"]
    else:
        raise ValueError("The comparisons file must contain at least 2 columns (file path and label/prefix).")

    gene_sets, loaded_lib_name = load_filtered_gene_sets(active_lib, args.organism)

    all_sig_results = []
    all_full_results = []

    for _, row in comp_list_df.iterrows():
        csv_file = str(row["dge_csv"]).strip()
        comp_prefix = str(row["comp_prefix"]).strip()

        if not os.path.exists(csv_file):
            print(f"  WARNING: File {csv_file} not found. Skipping comparison {comp_prefix}.")
            continue

        full_res, sig_res = run_single_gsea(
            csv_file,
            comp_prefix,
            gene_sets,
            args.qvalue,
            args.pvalue,
            args.prefix,
        )

        all_full_results.append(full_res)
        all_sig_results.append(sig_res)

    if not all_full_results:
        raise RuntimeError("No comparisons completed successfully.")

    master_sig = pd.concat(all_sig_results, axis=0, ignore_index=True)
    master_full = pd.concat(all_full_results, axis=0, ignore_index=True)

    print("\nCompiling unified comparative dotplot across comparisons...")

    top_terms = []
    for comp in master_sig["comparison"].unique():
        sub = master_sig[master_sig["comparison"] == comp]
        top_pos = sub[sub["NES"] > 0].nlargest(args.top_n, "NES")["Term"].tolist()
        top_neg = sub[sub["NES"] < 0].nsmallest(args.top_n, "NES")["Term"].tolist()
        top_terms.extend(top_pos + top_neg)

    selected_terms = list(set(top_terms))

    if not selected_terms:
        print("No significant pathways passed the specified thresholds. Unified plot skipped.")
        return

    plot_df = master_full[master_full["Term"].isin(selected_terms)].copy()
    plot_df["-log10(FDR)"] = -np.log10(plot_df["FDR q-val"].clip(lower=1e-10))
    plot_df["_gene_pct"] = plot_df["Gene %"].apply(parse_gene_percent)
    plot_df["Clean_Term"] = plot_df["Term"].apply(
        lambda x: textwrap.fill(str(x).replace("_", " "), width=45)
    )

    num_comps = len(plot_df["comparison"].unique())
    num_terms = len(selected_terms)

    fig_height = max(10, num_terms * 0.65)
    fig_width = max(12, num_comps * 3.0)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.scatterplot(
        data=plot_df,
        x="comparison",
        y="Clean_Term",
        size="_gene_pct",
        hue="NES",
        palette="vlag",
        sizes=(60, 350),
        hue_norm=(-3, 3),
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
    )

    plt.xticks(rotation=45, ha="right", fontsize=11)
    plt.yticks(fontsize=10)

    plt.title(
        f"Unified GSEA Pathway Analysis ({loaded_lib_name})\nOrganism: {args.organism} (Strict Case-Insensitive Human-Filtered)",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )

    plt.xlabel("Comparison", fontsize=12, labelpad=10)
    plt.ylabel("Pathway / Gene Set", fontsize=12)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=True,
        title="Metrics",
    )

    plt.grid(True, linestyle="--", alpha=0.3, axis="y")
    plt.tight_layout()

    out_plot = os.path.join("figures", f"{args.prefix}_unified_dotplot.png")
    plt.savefig(out_plot, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"\nUnified cross-comparison dotplot saved to: {out_plot}")


if __name__ == "__main__":
    main()
