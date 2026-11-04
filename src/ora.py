#!/usr/bin/env python3
import argparse
import pandas as pd
import gseapy as gp

def main():
    parser = argparse.ArgumentParser(description="Run ORA pathway enrichment on rescued genes.")
    parser.add_argument("--csv", required=True, help="Path to the rescued genes CSV file.")
    parser.add_argument("--prefix", required=True, help="Output prefix for results.")
    parser.add_argument("--genome", default="mouse", help="Organism/genome for Enrichr (default: mouse).")
    parser.add_argument("--gene-col", default="names", help="Column name containing gene symbols (default: names).")
    parser.add_argument("--cutoff", type=float, default=0.05, help="P-value cutoff for significant pathways (default: 0.05).")
    
    args = parser.parse_args()

    # Load CSV and extract unique genes
    df = pd.read_csv(args.csv)
    if args.gene_col not in df.columns:
        raise ValueError(f"Column '{args.gene_col}' not found in {args.csv}. Columns available: {list(df.columns)}")

    gene_list = df[args.gene_col].dropna().astype(str).unique().tolist()
    print(f"Loaded {len(gene_list)} unique genes from {args.csv}")

    if len(gene_list) == 0:
        print("Error: Gene list is empty. Aborting enrichment.")
        return

    # Standard pathway databases for mouse
    gene_sets = [
        'GO_Biological_Process_2023',
        'KEGG_2021_Mouse',
        'Reactome_2022'
    ]

    print(f"Running Enrichr ORA for organism: {args.genome}...")
    
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=args.genome,
            outdir=None,
            cutoff=args.cutoff
        )
        
        if enr.results is not None and not enr.results.empty:
            out_file = f"{args.prefix}_ora_results.csv"
            enr.results.to_csv(out_file, index=False)
            
            sig_file = f"{args.prefix}_significant_ora_results.csv"
            sig_df = enr.results[enr.results['Adjusted P-value'] < args.cutoff]
            sig_df.to_csv(sig_file, index=False)
            
            print(f"Success! All results saved to {out_file} and significant results saved to {sig_file}")
        else:
            print(f"Enrichment ran, but returned no significant pathways at p < {args.cutoff}.")
            
    except Exception as e:
        print(f"Error executing gseapy: {e}")

if __name__ == "__main__":
    main()
