import argparse
import scanpy as sc
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

with open(args.input, 'r') as f:
    files = [line.strip() for line in f if line.strip()]

adatas = []

for f in files:
    adata = sc.read_h5ad(f)

    # Extract experiment name before "_noDoublets.h5ad"
    filename = os.path.basename(f)
    experiment = filename.replace("_noDoublets.h5ad", "")

    adata.obs["experiment"] = experiment
    adatas.append(adata)

merged = sc.concat(
    adatas,
    join="inner",
    index_unique="-"
)

merged.write(args.output, compression="gzip")
