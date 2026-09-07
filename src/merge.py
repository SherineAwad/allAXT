import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

with open(args.input, 'r') as f:
    files = [line.strip() for line in f if line.strip()]

adatas = [sc.read_h5ad(f) for f in files]

merged = sc.concat(adatas, join="outer", fill_value=0, index_unique="-")

merged.write(args.output, compression="gzip")
