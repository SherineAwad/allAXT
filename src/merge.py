import argparse
import os
import numpy as np
import scanpy as sc
from scipy import sparse

# 1. Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# 2. Read Input File List
with open(args.input, "r") as f:
    files = [line.strip() for line in f if line.strip()]

adatas = []

# 3. Read AnnDatas and Attach Metadata
for f in files:
    adata = sc.read_h5ad(f)

    filename = os.path.basename(f)
    experiment = filename.replace("_noDoublets.h5ad", "")

    adata.obs["experiment"] = experiment
    adatas.append(adata)

# 4. Merge All Datasets (Outer Join Keeps All Genes, Including EGFP)
merged = sc.concat(adatas, join="outer", index_unique="-")

# 5. Fill Missing Outer-Join Values with 0 (Required for Sparse Matrices)
if sparse.issparse(merged.X):
    merged.X.data = np.nan_to_num(merged.X.data, copy=False)
    merged.X.eliminate_zeros()
else:
    merged.X = np.nan_to_num(merged.X, copy=False)

# 6. Ensure Unique Cell Barcodes
if not merged.obs_names.is_unique:
    merged.obs_names_make_unique()

# 7. Print Cell Counts to Log
if "sample" in merged.obs:
    print("\nUpdated sample counts:")
    print(merged.obs["sample"].value_counts())

if "experiment" in merged.obs:
    print("\nCell counts per experiment:")
    print(merged.obs["experiment"].value_counts())

# 8. Save Compressed Dataset
merged.write(args.output, compression="gzip")
print(f"\nSAVED MERGED DATASET → {args.output}")
