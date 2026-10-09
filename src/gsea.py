
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import warnings
warnings.filterwarnings('ignore')


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', required=True, help='DGE CSV file')
    parser.add_argument('--organism', default='mouse', help='Organism (default: mouse)')
    parser.add_argument('--resources', default='go', help='Gene set resource (default: go)')
    parser.add_argument('--prefix', default='GSEA_results', help='Prefix for output files')
    parser.add_argument('--qvalue', type=float, default=0.05,
                        help='FDR q-value threshold (default: 0.05)')
    return parser.parse_args()


def parse_gene_percent(value):
    """Convert GSEApy Gene % values, including fractions such as 10/200, to percentages."""
    if pd.isna(value):
        return np.nan

    value = str(value).strip()

    try:
        if '/' in value:
            numerator, denominator = value.split('/', 1)
            denominator = float(denominator)
            if denominator == 0:
                return np.nan
            return float(numerator) / denominator * 100

        if value.endswith('%'):
            return float(value.rstrip('%'))

        return float(value)

    except (ValueError, ZeroDivisionError):
        return np.nan


def main():
    args = parse_arguments()

    # Create figures directory if it doesn't exist
    os.makedirs('figures', exist_ok=True)

    # Load data
    df = pd.read_csv(args.csv)
    print(f"Loaded {len(df)} genes")

    # Rank genes by signed Wilcoxon scores
    df['ranking'] = df['scores']

    # Prepare ranking
    ranking = df[['gene', 'ranking']].dropna().sort_values(
        'ranking', ascending=False
    )
    ranking.columns = ['Gene', 'Score']

    dup_count = ranking['Score'].duplicated().sum()
    print(f"Duplicate scores: {dup_count} ({dup_count/len(ranking)*100:.1f}%)")

    if dup_count > 0:
        np.random.seed(42)
        ranking['Score'] = ranking['Score'] + np.random.normal(
            0, 1e-10, len(ranking)
        )

    # Get all matching gene-set libraries
    from gseapy.parser import get_library_name
    libraries = get_library_name(organism=args.organism.capitalize())

    resource = args.resources.strip().lower()

    if resource in ['go', 'go_bp', 'gobp', 'biological_process']:
        matched = [
            lib for lib in libraries
            if 'GO_Biological_Process' in lib
        ]
    elif resource == 'reactome':
        matched = [
            lib for lib in libraries
            if 'reactome' in lib.lower()
        ]
    elif resource == 'kegg':
        matched = [
            lib for lib in libraries
            if 'kegg' in lib.lower()
        ]
    elif resource == 'hallmark':
        matched = [
            lib for lib in libraries
            if 'hallmark' in lib.lower()
        ]
    else:
        matched = [
            lib for lib in libraries
            if resource in lib.lower()
        ]

    if not matched:
        raise ValueError(
            f"No available gene-set library matches '{args.resources}' "
            f"for organism '{args.organism}'. Available libraries must "
            f"be checked using GSEApy's get_library_name()."
        )

    print(f"Matched {len(matched)} libraries:")
    for lib in matched:
        print(f"  {lib}")

    # Load and combine gene sets from ALL matching libraries
    gene_set = {}

    for lib in matched:
        library_gene_sets = gp.get_library(
            name=lib,
            organism=args.organism.capitalize()
        )

        for term, genes in library_gene_sets.items():
            combined_term = f"{lib}::{term}"
            gene_set[combined_term] = genes

    if not gene_set:
        raise ValueError(
            f"No gene sets could be loaded from matching libraries: {matched}"
        )

    print(f"Combined total: {len(gene_set)} gene sets")

    # Run GSEA
    results = gp.prerank(
        rnk=ranking,
        gene_sets=gene_set,
        permutation_num=1000,
        outdir=f'{args.prefix}_gsea_output',
        min_size=15,
        max_size=500,
        seed=42
    )

    # Convert numeric columns
    results_df = results.res2d
    numeric_cols = [
        'ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val', 'Tag %'
    ]
    for col in numeric_cols:
        if col in results_df.columns:
            results_df[col] = pd.to_numeric(
                results_df[col], errors='coerce'
            )

    # Save all pathway results
    results_df.to_csv(f'{args.prefix}_results.csv', index=False)
    print(f"\nAll pathway results saved to {args.prefix}_results.csv")
    print(f"Found {len(results_df)} pathways")

    # Filter significant pathways by FDR q-value
    significant_df = results_df[
        results_df['FDR q-val'] < args.qvalue
    ].copy()

    significant_df = significant_df.sort_values(
        'NES', ascending=False
    )

    # Save significant pathway results
    significant_df.to_csv(
        f'{args.prefix}_significant_results.csv', index=False
    )
    print(
        f"Significant pathway results saved to "
        f"{args.prefix}_significant_results.csv"
    )
    print(
        f"Found {len(significant_df)} significant pathways "
        f"(FDR q-value < {args.qvalue})"
    )

    # Stop if no pathways pass the FDR threshold
    if significant_df.empty:
        print(
            f"No pathways passed FDR q-value < {args.qvalue}. "
            "Skipping dotplot and top-pathway reporting."
        )
        return

    # Prepare significant pathways for plotting
    results_df = significant_df.dropna(subset=['NES'])

    # Get top 15 enriched and top 15 depleted significant pathways
    top_enriched = results_df[
        results_df['NES'] > 0
    ].sort_values('NES', ascending=False).head(15)

    top_depleted = results_df[
        results_df['NES'] < 0
    ].sort_values('NES', ascending=True).head(15)

    # Combine for plotting
    plot_df = pd.concat([top_enriched, top_depleted]).copy()

    # Stop if there are no pathways with a valid, non-zero NES
    if plot_df.empty:
        print(
            "No significant pathways with positive or negative NES "
            "are available to plot."
        )
        return

    # Add -log10(FDR) for colour
    plot_df['-log10(FDR)'] = -np.log10(
        plot_df['FDR q-val'].clip(lower=1e-10)
    )

    # Parse Gene % for dot size
    plot_df['_gene_pct'] = plot_df['Gene %'].apply(parse_gene_percent)

    # Handle missing Gene % values
    valid_gene_pct = plot_df['_gene_pct'].dropna()

    if valid_gene_pct.empty:
        print(
            "Gene % could not be parsed; using equal dot sizes "
            "because pathway coverage is unavailable."
        )
        sizes = np.full(len(plot_df), 275.0)
        size_legend_values = []
    else:
        min_gene_pct = valid_gene_pct.min()
        max_gene_pct = valid_gene_pct.max()

        plot_df['_gene_pct'] = plot_df['_gene_pct'].fillna(min_gene_pct)

        if max_gene_pct == min_gene_pct:
            sizes = np.full(len(plot_df), 275.0)
            size_legend_values = [min_gene_pct]
        else:
            sizes = 50 + (
                (plot_df['_gene_pct'] - min_gene_pct)
                / (max_gene_pct - min_gene_pct)
            ) * 450

            size_legend_values = [
                min_gene_pct,
                (min_gene_pct + max_gene_pct) / 2,
                max_gene_pct
            ]

    # Create the dotplot
    print("\nGenerating dotplot...")
    fig, ax = plt.subplots(figsize=(12, 10))

    # X-axis = NES; dot size = Gene % (pathway coverage)
    scatter = ax.scatter(
        plot_df['NES'],
        range(len(plot_df)),
        s=sizes,
        c=plot_df['-log10(FDR)'],
        cmap='RdYlBu_r',
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )

    # Add colourbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(FDR)', fontsize=10)

    # Add vertical line at 0
    ax.axvline(
        x=0, color='black', linestyle='--',
        linewidth=0.8, alpha=0.5
    )

    # Set y-axis labels
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels([
        t[:50] + '...' if len(t) > 50 else t
        for t in plot_df['Term'].values
    ], fontsize=8)

    # Labels and title
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    ax.set_ylabel('Pathways', fontsize=12)
    ax.set_title(
        f'GSEA Results: Top 30 Pathways\n'
        f'{args.organism.upper()} {args.resources.upper()}',
        fontsize=14, fontweight='bold'
    )

    # Add grid
    ax.grid(True, alpha=0.3, axis='x')

    # Invert y-axis to show top enriched at top
    ax.invert_yaxis()

    # Add size legend matching the actual Gene % size mapping
    from matplotlib.lines import Line2D

    if size_legend_values:
        legend_elements = []

        for gene_pct in size_legend_values:
            if max_gene_pct == min_gene_pct:
                marker_area = 275.0
            else:
                marker_area = 50 + (
                    (gene_pct - min_gene_pct)
                    / (max_gene_pct - min_gene_pct)
                ) * 450

            legend_elements.append(
                Line2D(
                    [0], [0],
                    marker='o',
                    color='w',
                    label=f'{gene_pct:.1f}%',
                    markerfacecolor='gray',
                    markersize=np.sqrt(marker_area),
                    markeredgecolor='black'
                )
            )

        ax.legend(
            handles=legend_elements,
            loc='lower right',
            title='Gene %',
            fontsize=8
        )

    plt.tight_layout()

    # Save figure
    output_path = 'figures/' + args.prefix + '_dotplot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.show()

    # Print top results to console
    print("\n" + "="*60)
    print(
        f"TOP 10 ENRICHED PATHWAYS "
        f"(by NES, FDR < {args.qvalue}):"
    )
    print("="*60)

    for idx, row in top_enriched.head(10).iterrows():
        fdr = row['FDR q-val']
        print(
            f"{row['Term'][:70]}: "
            f"NES={row['NES']:.3f}, FDR={fdr}"
        )

    print("\n" + "="*60)
    print(
        f"TOP 10 DEPLETED PATHWAYS "
        f"(by NES, FDR < {args.qvalue}):"
    )
    print("="*60)

    for idx, row in top_depleted.head(10).iterrows():
        fdr = row['FDR q-val']
        print(
            f"{row['Term'][:70]}: "
            f"NES={row['NES']:.3f}, FDR={fdr}"
        )


if __name__ == "__main__":
    main()
