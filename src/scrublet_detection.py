import scrublet as scr
import scanpy as sc
import numpy as np
from scipy import sparse
import argparse
import os
import matplotlib.pyplot as plt

# -------------------------
# Args (UNCHANGED)
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--prefix', required=True)
parser.add_argument('--expected_doublet_rate', type=float, default=0.07)
parser.add_argument('--min_gene_variability_pctl', type=float, default=75)

args = parser.parse_args()

# -------------------------
# Load data
# -------------------------
adata = sc.read_h5ad(args.input)

if not adata.obs_names.is_unique:
    adata.obs_names_make_unique()

print(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")

# -------------------------
# Prepare matrix
# -------------------------
X = adata.X
if sparse.issparse(X):
    X = X.toarray().astype(np.float32)
else:
    X = np.array(X, dtype=np.float32)

# -------------------------
# Run Scrublet
# -------------------------
print("Running Scrublet...")

scrub = scr.Scrublet(
    X,
    expected_doublet_rate=args.expected_doublet_rate
)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=args.min_gene_variability_pctl,
    n_prin_comps=30
)

# -------------------------
# USE SCRUBLET THRESHOLD
# -------------------------
final_doublets = predicted_doublets

# -------------------------
# STORE RESULTS (BEFORE FILTERING)
# -------------------------
adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = final_doublets.astype(bool)

print(f"Detected doublets: {final_doublets.sum()} / {adata.n_obs}")

# -------------------------
# PLOT FIRST (IMPORTANT FIX)
# -------------------------
os.makedirs("figures", exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(doublet_scores, bins=50, edgecolor='black')
ax.set_title(f"{args.prefix} Scrublet scores")
ax.set_xlabel("Doublet score")
ax.set_ylabel("Cells")

fig.savefig(f"figures/{args.prefix}_scrublet_scores.png", dpi=300)
plt.close(fig)

# -------------------------
# REMOVE DOUBLETS (AFTER PLOT)
# -------------------------
adata = adata[~final_doublets].copy()

print(f"Remaining cells after filtering: {adata.n_obs}")

# -------------------------
# SAVE FINAL OUTPUT
# -------------------------
adata.write(args.output, compression="gzip")

print(f"SAVED CLEAN DATASET → {args.output}")

