# AXT project - merging different experiments


## Single-cell RNA-seq quality control

Single-cell RNA-seq data from multiple samples were combined and subjected to quality control to remove low-quality cells and genes that could interfere with downstream analysis.

Cells were retained if they had **500–10,000 detected genes**, **1,000–50,000 total RNA counts**, and **less than 10% mitochondrial RNA**. These thresholds were used to exclude cells with very low RNA content, unusually high RNA content, or a high proportion of mitochondrial transcripts, which can indicate poor-quality or stressed cells.

Genes detected in fewer than **3 cells** were also removed because they provide very limited information for downstream analysis.

The resulting dataset was used for subsequent single-cell analysis.

Below are the QC plots showing the distribution of cell-quality metrics before and after filtering for each dataset. 

#### GSE135985

<img src="figures/violin_GSE135985_preQC.png?v=1" width="48%"> <img src="figures/violin_GSE135985_AfterQC.png?v=1" width="48%">

#### AXT

<img src="figures/violin_AXT_preQC.png?v=1" width="48%"> <img src="figures/violin_AXT_AfterQC.png?v=1" width="48%">

#### SLX-28276

<img src="figures/violin_SLX-28276_preQC.png?v=1" width="48%"> <img src="figures/violin_SLX-28276_AfterQC.png?v=1" width="48%">


## Doublet detection

Following quality control, potential **doublets** were identified using Scrublet. Doublets are droplets containing two cells, which can produce mixed expression profiles and potentially create artificial cell populations.

Scrublet assigned each cell a doublet score based on its expression profile, and cells classified as doublets were removed before downstream analysis.

The following parameters were used for **all datasets**:
- Minimum counts: **2**
- Minimum cells: **3**
- Minimum gene variability percentile: **75**
- Number of principal components: **30**

Dataset-specific parameters were:

| Dataset | Expected doublet rate | Doublet score threshold |
|---|---:|---:|
| SLX-28276 | 0.08 (8%) | Scrublet automatic |
| AXT | 0.06 (6%) | 0.25 |
| GSE135985 | 0.03 (3%) | 0.25 |

Cells identified as doublets using these criteria were excluded from the datasets used for subsequent single-cell analysis.

<img src="figures/AXT_scrublet_scores.png?v=1" width="32%"> <img src="figures/GSE135985_scrublet_scores.png?v=1" width="32%"> <img src="figures/SLX-28276_scrublet_scores.png?v=1" width="32%">

| Dataset | Cells Before | Predicted Doublets | Predicted Doublet Rate | Cells After |
|---|---:|---:|---:|---:|
| SLX-28276 | 47,188 | 921 | 2.0% | 46,267 |
| AXT | 32,593 | 446 | 1.37% | 32,147 |
| GSE135985 | 13,784 | 29 | 0.21% | 13,755 |

## Merging datasets

The quality-controlled datasets were combined into a single dataset for joint downstream analysis. Each dataset was labelled with its corresponding experiment of origin so that cells could be tracked back to their source.

Only genes shared across all datasets were retained, ensuring that the combined dataset contained a consistent set of features across experiments. Cell identifiers were also made unique during merging.

The resulting combined dataset was saved for downstream analysis.

The merged dataset contained **92,169 cells** across three experiments:

| Experiment | Cells |
|---|---:|
| SLX-28276 | 46,267 |
| AXT | 32,147 |
| GSE135985 | 13,755 |

The cells were distributed across the following samples:

| Sample | Count |
|---|---:|
| nonReg_4wpa | 23,948 |
| Reg_4wpa | 22,319 |
| nonReg | 20,103 |
| Reg | 12,044 |
| Uninjured2 | 5,638 |
| Non_Regen_14_DPA | 4,086 |
| Uninjured1 | 2,610 |
| Regen_14_DPA | 1,421 |

Only genes shared across all three experiments were retained to provide a consistent feature set for joint analysis.

## Normalisation, feature selection and dimensionality reduction

Following doublet removal, the data were normalised to account for differences in sequencing depth between cells and then log-transformed to reduce the influence of highly abundant transcripts.

The original raw counts were retained separately for subsequent analyses. Highly variable genes were then identified to focus downstream analysis on genes showing the most informative variation between cells.

The following parameters were used:
- Normalisation target: **10,000 counts per cell**
- Highly variable genes: **5,000 genes**
- Highly variable gene selection method: **Seurat v3**
- Number of principal components: **30**
- PCA scaling: **maximum absolute value of 10**
- PCA solver: **ARPACK**

The 30-dimensional PCA representation was then used to construct a cell-to-cell neighbourhood graph and generate a **UMAP embedding**, allowing the overall structure and relationships between cells to be visualised.

The resulting dataset was saved for downstream single-cell analysis.


![](figures/umap_allAXT_umap.png?v=1) 

#### Per sample 
<img src="figures/umap_allAXT_Uninjured1.png?v=1" width="48%"> <img src="figures/umap_allAXT_Uninjured2.png?v=1" width="48%">

<img src="figures/umap_allAXT_nonReg_4wpa.png?v=1" width="48%"> <img src="figures/umap_allAXT_Non_Regen_14_DPA.png?v=1" width="48%">

<img src="figures/umap_allAXT_nonReg.png?v=1" width="48%"> <img src="figures/umap_allAXT_Reg.png?v=1" width="48%">

<img src="figures/umap_allAXT_Reg_4wpa.png?v=1" width="48%"> <img src="figures/umap_allAXT_Regen_14_DPA.png?v=1" width="48%">


## Batch correction using Harmony

To reduce differences between experiments that could arise from technical variation rather than biological differences, **Harmony** was used to integrate the three datasets based on their experiment of origin.

Harmony was applied to the **30-dimensional PCA representation**, using experiment as the batch variable. The Harmony-corrected representation was then used to construct a new cell-to-cell neighbourhood graph and generate a **UMAP embedding**.

The resulting UMAP was visualised by sample to assess the distribution of cells across samples after batch correction.

The Harmony-corrected dataset was saved for subsequent single-cell analysis.

 
![](figures/umap_allAXT_harmony_sample.png?v=1)

#### Per sample: 

<img src="figures/umap_allAXT_harmony_Reg_4wpa.png?v=1" width="48%"> <img src="figures/umap_allAXT_harmony_nonReg_4wpa.png?v=1" width="48%">

<img src="figures/umap_allAXT_harmony_Reg.png?v=1" width="48%"> <img src="figures/umap_allAXT_harmony_nonReg.png?v=1" width="48%">

<img src="figures/umap_allAXT_harmony_Non_Regen_14_DPA.png?v=1" width="48%"> <img src="figures/umap_allAXT_harmony_Regen_14_DPA.png?v=1" width="48%">

<img src="figures/umap_allAXT_harmony_Uninjured1.png?v=1" width="48%"> <img src="figures/umap_allAXT_harmony_Uninjured2.png?v=1" width="48%">

## Clustering

Following batch correction, cells were grouped into clusters based on their **similarity in gene expression profiles**. Leiden clustering was performed using the Harmony-corrected neighbourhood graph.

A clustering **resolution of 2.5** was used to identify relatively fine-grained cell populations.

The resulting clusters were visualised on the UMAP embedding, with cluster identities displayed directly on the plot.

To assess whether clustering was associated with differences in cell quality, the distributions of **detected genes, total RNA counts, and mitochondrial RNA percentage** were also examined across clusters.

![](figures/umap_allAXT_leiden.png?v=1)

#### Quality per cluster 
<img src="figures/violin_allAXT_QC_n_genes_by_counts.png?v=1" width="32%"> <img src="figures/violin_allAXT_QC_total_counts.png?v=1" width="32%"> <img src="figures/violin_allAXT_QC_pct_counts_mt.png?v=1" width="32%">


