import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description="Rename sample categories in AnnData obs metadata.")
parser.add_argument('--input', required=True, help="Path to input .h5ad file.")
parser.add_argument('--output', required=True, help="Path to save updated .h5ad file.")
args = parser.parse_args()

print(f"Loading {args.input}...")
adata = sc.read_h5ad(args.input)

# Check if 'sample' column exists
if 'sample' not in adata.obs.columns:
    raise KeyError("The column 'sample' was not found in adata.obs.")

# Define mapping
mapping = {
    'SITTB10': 'nonReg_4wpa',
    'SITTB2': 'Reg_4wpa'
}

# Ensure string type to avoid categorical assignment issues, then convert back to category
adata.obs['sample'] = adata.obs['sample'].astype(str).replace(mapping).astype('category')

print("\nUpdated sample counts:")
print(adata.obs['sample'].value_counts())

print(f"\nSaving updated dataset to {args.output}...")
adata.write(args.output, compression="gzip")
print("DONE.")
