import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="Inspect AnnData metadata for Harmony integration compatibility.")
parser.add_argument('--input', required=True, help="Path to the merged .h5ad file.")
args = parser.parse_args()

print(f"Loading {args.input}...")
adata = sc.read_h5ad(args.input)

print("\n==========================================")
print("DATASET OVERVIEW")
print("==========================================")
print(f"Total Cells (obs): {adata.n_obs}")
print(f"Total Genes (var): {adata.n_vars}")

print("\n==========================================")
print("METADATA COLUMNS (adata.obs)")
print("==========================================")
for col in adata.obs.columns:
    n_unique = adata.obs[col].nunique()
    print(f"  - '{col}': {n_unique} unique categories/values")

print("\n==========================================")
print("SAMPLE/BATCH COUNTS PER CATEGORICAL COLUMN")
print("==========================================")
for col in adata.obs.select_dtypes(include=['category', 'object']).columns:
    print(f"\nBreakdown for obs['{col}']:")
    print(adata.obs[col].value_counts().head(10))

print("\n==========================================")
print("HARMONY READINESS CHECK")
print("==========================================")
if 'batch' in adata.obs.columns:
    print("✓ FOUND 'batch' column in adata.obs.")
    print(f"  Unique batches detected: {adata.obs['batch'].nunique()}")
else:
    print("✗ NO 'batch' column found in adata.obs.")
    print("  Available columns you can use as batch_key in Harmony:", list(adata.obs.columns))
