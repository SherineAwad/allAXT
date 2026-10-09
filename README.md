# AXT project - merging different experiments


## Single-cell RNA-seq quality control

Single-cell RNA-seq data from multiple samples were combined and subjected to quality control to remove low-quality cells and genes that could interfere with downstream analysis.

Cells were retained if they had **500–10,000 detected genes**, **1,000–50,000 total RNA counts**, and **less than 10% mitochondrial RNA**. These thresholds were used to exclude cells with very low RNA content, unusually high RNA content, or a high proportion of mitochondrial transcripts, which can indicate poor-quality or stressed cells.

Genes detected in fewer than **3 cells** were also removed because they provide very limited information for downstream analysis.

The resulting dataset was used for subsequent single-cell analysis.

Below are the QC plots showing the distribution of cell-quality metrics before and after filtering for each dataset. 

#### GSE135985

<img src="figures/violin_GSE135985_preQC.png?v=3" width="48%"> <img src="figures/violin_GSE135985_AfterQC.png?v=3" width="48%">

#### AXT

<img src="figures/violin_AXT_preQC.png?v=3" width="48%"> <img src="figures/violin_AXT_AfterQC.png?v=3" width="48%">

#### SLX-28276

<img src="figures/violin_SLX-28276_preQC.png?v=3" width="48%"> <img src="figures/violin_SLX-28276_AfterQC.png?v=3" width="48%">


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

<img src="figures/AXT_scrublet_scores.png?v=3" width="32%"> <img src="figures/GSE135985_scrublet_scores.png?v=3" width="32%"> <img src="figures/SLX-28276_scrublet_scores.png?v=3" width="32%">

| Dataset | Cells Before | Predicted Doublets | Predicted Doublet Rate | Cells After |
|---|---:|---:|---:|---:|
| SLX-28276 | 47,188 | 921 | 2.0% | 46,267 |
| AXT | 32,593 | 446 | 1.37% | 32,147 |
| GSE135985 | 13,784 | 29 | 0.21% | 13,755 |

## Merging datasets

The quality-controlled datasets were combined into a single dataset for joint downstream analysis. Each dataset was labelled with its corresponding experiment of origin so that cells could be tracked back to their source.

All genes across all datasets were retained using an outer join (including dataset-specific genes such as *EGFP*), providing a comprehensive feature set for joint analysis.

The resulting combined dataset was saved for downstream analysis.

The merged dataset contained 92,169 cells across three experiments:

| Experiment | Cells |
|---|---:|
| SLX-28276 | 46,267 |
| AXT | 32,147 |
| GSE135985 | 13,755 |


The cells were distributed across the following samples:

| Sample | Cells |
|---|---:|
| nonReg_4wpa | 23,948 |
| Reg_4wpa | 22,319 |
| nonReg | 20,103 |
| Reg | 12,044 |
| Uninjured2 | 5,638 |
| Non_Regen_14_DPA | 4,086 |
| Uninjured1 | 2,610 |
| Regen_14_DPA | 1,421 |


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


![](figures/umap_allAXT_umap.png?v=3) 

#### Per sample 
<img src="figures/umap_allAXT_Uninjured1.png?v=3" width="48%"> <img src="figures/umap_allAXT_Uninjured2.png?v=3" width="48%">

<img src="figures/umap_allAXT_nonReg_4wpa.png?v=3" width="48%"> <img src="figures/umap_allAXT_Non_Regen_14_DPA.png?v=3" width="48%">

<img src="figures/umap_allAXT_nonReg.png?v=3" width="48%"> <img src="figures/umap_allAXT_Reg.png?v=3" width="48%">

<img src="figures/umap_allAXT_Reg_4wpa.png?v=3" width="48%"> <img src="figures/umap_allAXT_Regen_14_DPA.png?v=3" width="48%">


## Batch correction using Harmony

To reduce differences between experiments that could arise from technical variation rather than biological differences, **Harmony** was used to integrate the three datasets based on their experiment of origin.

Harmony was applied to the **30-dimensional PCA representation**, using experiment as the batch variable. The Harmony-corrected representation was then used to construct a new cell-to-cell neighbourhood graph and generate a **UMAP embedding**.

The resulting UMAP was visualised by sample to assess the distribution of cells across samples after batch correction.

The Harmony-corrected dataset was saved for subsequent single-cell analysis.

 
![](figures/umap_allAXT_harmony_sample.png?v=3)

#### Per sample: 

<img src="figures/umap_allAXT_harmony_Reg_4wpa.png?v=3" width="48%"> <img src="figures/umap_allAXT_harmony_nonReg_4wpa.png?v=3" width="48%">

<img src="figures/umap_allAXT_harmony_Reg.png?v=3" width="48%"> <img src="figures/umap_allAXT_harmony_nonReg.png?v=3" width="48%">

<img src="figures/umap_allAXT_harmony_Non_Regen_14_DPA.png?v=3" width="48%"> <img src="figures/umap_allAXT_harmony_Regen_14_DPA.png?v=3" width="48%">

<img src="figures/umap_allAXT_harmony_Uninjured1.png?v=3" width="48%"> <img src="figures/umap_allAXT_harmony_Uninjured2.png?v=3" width="48%">

## Clustering

Following batch correction, cells were grouped into clusters based on their **similarity in gene expression profiles**. Leiden clustering was performed using the Harmony-corrected neighbourhood graph.

A clustering **resolution of 2.5** was used to identify relatively fine-grained cell populations.

The resulting clusters were visualised on the UMAP embedding, with cluster identities displayed directly on the plot.

To assess whether clustering was associated with differences in cell quality, the distributions of **detected genes, total RNA counts, and mitochondrial RNA percentage** were also examined across clusters.

![](figures/umap_allAXT_leiden.png?v=3)

#### Quality per cluster 
<img src="figures/violin_allAXT_QC_n_genes_by_counts.png?v=3" width="32%"> <img src="figures/violin_allAXT_QC_total_counts.png?v=3" width="32%"> <img src="figures/violin_allAXT_QC_pct_counts_mt.png?v=3" width="32%">



#### Feature plots and dotplot

##### Marker genes
```python 
marker_genes = {
    "Fibroblast": ["Prrx1", "Pdgfra", "Col1a1", "Dcn", "Pi16", "Cd34"],
    "Endothelial": ["Esam", "Flt1", "Vwf", "Plvap", "Cdh5", "Pecam1", "Kdr", "Emcn", "Erg", "Cd34"],
    "Macrophage": ["Adgre1", "Csf1r", "Cd68", "Mrc1", "Cd163"],
    "Keratinocyte": ["Krt14", "Krt5", "Epcam", "Cdh1", "Krt17", "Dsg3"],
    "Osteoblast": ["Runx2", "Postn", "Mmp13", "Spp1", "Dmp1", "Sp7", "Bglap", "Alpl", "Ibsp", "Col1a1"],
    "Pericyte": ["Cspg4", "Pdgfrb", "Kcnj8", "Abcc9", "Rgs5"],
    "SMC": ["Myh11", "Tagln", "Acta2", "Cnn1", "Des"],
    "Chondrocyte": ["Col2a1", "Acan", "Sox9", "Matn1", "Col9a1", "Comp"],
    "Schwann": ["Mbp", "Plp1", "Sox10", "S100b", "Pmp22"],
    "T-cell": ["Cd3d", "Cd3e", "Cd4", "Cd8a", "Cd28"],
    "Osteoclast": ["Ctsk", "Acp5", "Calcr", "Oscar", "Nfatc1", "Dcstamp", "Tnfrsf11a"],
    "Synoviocyte": ["Prg4", "Ucma", "Gdf5", "Cilp2", "Frzb"],
    "Neutrophil": ["Ly6g", "S100a8", "S100a9", "Mpo", "Csf3r"],
    "Lymphatic_Endothelial": ["Prox1", "Lyve1", "Pdpn", "Flt4", "Pecam1"],
    "B-cell": ["Cd79a", "Cd19", "Ms4a1", "Ighm", "Pax5", "Cd22"],
    "Rspo3_Col23a1": ["Rspo3", "Col23a1"],
    "MSC": ["Lepr", "Cxcl12", "Ngfr", "Nes", "Cd44", "Scf"],
    "Osteosarcoma": ["EGFP"],
    "Nail_Epithelium": ["Lgr6", "Sp6", "Sp8"],
    "Sweat_glands": ["Aqp5", "Scnn1a", "Scnn1b", "Scnn1g", "Krt19", "Krt7", "Krt8", "Krt18", "Krt5", "Krt14", "Foxa1"],
    "Mast cell": ["Kit", "Cpa3", "Tpsab1", "Tpsb2", "Ms4a2", "Hdc", "Hpgds", "Mcpt8", "Cd200r3", "Ccr3"]
}
```
![](figures/dotplot__allAXT_dotplot.png?v=1)


<img src="figures/umap_allAXT_Matn1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Oscar.png?v=3" width="32%"> <img src="figures/umap_allAXT_Pi16.png?v=3" width="32%">

<img src="figures/umap_allAXT_Epcam.png?v=3" width="32%"> <img src="figures/umap_allAXT_Ngfr.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cdh5.png?v=3" width="32%">

<img src="figures/umap_allAXT_Ighm.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd22.png?v=3" width="32%"> <img src="figures/umap_allAXT_Scnn1a.png?v=3" width="32%">

<img src="figures/umap_allAXT_Lgr6.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cilp2.png?v=3" width="32%"> <img src="figures/umap_allAXT_Dsg3.png?v=3" width="32%">

<img src="figures/umap_allAXT_Krt14.png?v=3" width="32%"> <img src="figures/umap_allAXT_Rspo3.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd68.png?v=3" width="32%">

<img src="figures/umap_allAXT_Ibsp.png?v=3" width="32%"> <img src="figures/umap_allAXT_Ly6g.png?v=3" width="32%"> <img src="figures/umap_allAXT_Krt18.png?v=3" width="32%">

<img src="figures/umap_allAXT_Col1a1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Dcstamp.png?v=3" width="32%"> <img src="figures/umap_allAXT_Foxa1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Pdgfrb.png?v=3" width="32%"> <img src="figures/umap_allAXT_Col23a1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Pmp22.png?v=3" width="32%">

<img src="figures/umap_allAXT_Krt7.png?v=3" width="32%"> <img src="figures/umap_allAXT_Flt1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cxcl12.png?v=3" width="32%">

<img src="figures/umap_allAXT_S100b.png?v=3" width="32%"> <img src="figures/umap_allAXT_Dmp1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cpa3.png?v=3" width="32%">

<img src="figures/umap_allAXT_Krt8.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd163.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cspg4.png?v=3" width="32%">

<img src="figures/umap_allAXT_Ctsk.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd200r3.png?v=3" width="32%"> <img src="figures/umap_allAXT_Nfatc1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Plp1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd3e.png?v=3" width="32%"> <img src="figures/umap_allAXT_Pecam1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Pdpn.png?v=3" width="32%"> <img src="figures/umap_allAXT_Kcnj8.png?v=3" width="32%"> <img src="figures/umap_allAXT_Mrc1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Pdgfra.png?v=3" width="32%"> <img src="figures/umap_allAXT_Krt5.png?v=3" width="32%"> <img src="figures/umap_allAXT_Kdr.png?v=3" width="32%">

<img src="figures/umap_allAXT_Krt19.png?v=3" width="32%"> <img src="figures/umap_allAXT_Flt4.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd3d.png?v=3" width="32%">

<img src="figures/umap_allAXT_Abcc9.png?v=3" width="32%"> <img src="figures/umap_allAXT_Frzb.png?v=3" width="32%"> <img src="figures/umap_allAXT_Ccr3.png?v=3" width="32%">

<img src="figures/umap_allAXT_EGFP.png?v=3" width="32%"> <img src="figures/umap_allAXT_Hdc.png?v=3" width="32%"> <img src="figures/umap_allAXT_S100a9.png?v=3" width="32%">

<img src="figures/umap_allAXT_Acp5.png?v=3" width="32%"> <img src="figures/umap_allAXT_Scnn1b.png?v=3" width="32%"> <img src="figures/umap_allAXT_Lyve1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Erg.png?v=3" width="32%"> <img src="figures/umap_allAXT_Ucma.png?v=3" width="32%"> <img src="figures/umap_allAXT_Adgre1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Emcn.png?v=3" width="32%"> <img src="figures/umap_allAXT_Calcr.png?v=3" width="32%"> <img src="figures/umap_allAXT_Pax5.png?v=3" width="32%">

<img src="figures/umap_allAXT_Sox9.png?v=3" width="32%"> <img src="figures/umap_allAXT_Tnfrsf11a.png?v=3" width="32%"> <img src="figures/umap_allAXT_Sp6.png?v=3" width="32%">

<img src="figures/umap_allAXT_Prox1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd8a.png?v=3" width="32%"> <img src="figures/umap_allAXT_S100a8.png?v=3" width="32%">

<img src="figures/umap_allAXT_Tpsab1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Runx2.png?v=3" width="32%"> <img src="figures/umap_allAXT_Alpl.png?v=3" width="32%">

<img src="figures/umap_allAXT_Cd34.png?v=3" width="32%"> <img src="figures/umap_allAXT_Spp1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Myh11.png?v=3" width="32%">

<img src="figures/umap_allAXT_Des.png?v=3" width="32%"> <img src="figures/umap_allAXT_Sp8.png?v=3" width="32%"> <img src="figures/umap_allAXT_Scnn1g.png?v=3" width="32%">

<img src="figures/umap_allAXT_Cdh1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Sp7.png?v=3" width="32%"> <img src="figures/umap_allAXT_Mcpt8.png?v=3" width="32%">

<img src="figures/umap_allAXT_Mbp.png?v=3" width="32%"> <img src="figures/umap_allAXT_Bglap.png?v=3" width="32%"> <img src="figures/umap_allAXT_Postn.png?v=3" width="32%">

<img src="figures/umap_allAXT_Prg4.png?v=3" width="32%"> <img src="figures/umap_allAXT_Lepr.png?v=3" width="32%"> <img src="figures/umap_allAXT_Ms4a2.png?v=3" width="32%">

<img src="figures/umap_allAXT_Acan.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd44.png?v=3" width="32%"> <img src="figures/umap_allAXT_Csf1r.png?v=3" width="32%">

<img src="figures/umap_allAXT_Dcn.png?v=3" width="32%"> <img src="figures/umap_allAXT_Krt17.png?v=3" width="32%"> <img src="figures/umap_allAXT_Aqp5.png?v=3" width="32%">

<img src="figures/umap_allAXT_Cd28.png?v=3" width="32%"> <img src="figures/umap_allAXT_Col2a1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd4.png?v=3" width="32%">

<img src="figures/umap_allAXT_Cnn1.png?v=3" width="32%"> <img src="figures/umap_allAXT_Acta2.png?v=3" width="32%"> <img src="figures/umap_allAXT_Cd79a.png?v=3" width="32%">

<img src="figures/umap_allAXT_Mmp13.png?v=3" width="32%"> <img src="figures/umap_allAXT_Gdf5.png?v=3" width="32%"> <img src="figures/umap_allAXT_Sox10.png?v=3" width="32%">

<img src="figures/umap_allAXT_Comp.png?v=3" width="32%"> <img src="figures/umap_allAXT_Mpo.png?v=3" width="32%"> <img src="figures/umap_allAXT_Kit.png?v=3" width="32%">

<img src="figures/umap_allAXT_Hpgds.png?v=3" width="32%"> <img src="figures/umap_allAXT_Nes.png?v=3" width="32%"> <img src="figures/umap_allAXT_Col9a1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Tagln.png?v=3" width="32%"> <img src="figures/umap_allAXT_Vwf.png?v=3" width="32%"> <img src="figures/umap_allAXT_Plvap.png?v=3" width="32%">

<img src="figures/umap_allAXT_Csf3r.png?v=3" width="32%"> <img src="figures/umap_allAXT_Tpsb2.png?v=3" width="32%"> <img src="figures/umap_allAXT_Ms4a1.png?v=3" width="32%">

<img src="figures/umap_allAXT_Cd19.png?v=3" width="32%"> <img src="figures/umap_allAXT_Esam.png?v=3" width="32%"> <img src="figures/umap_allAXT_Rgs5.png?v=3" width="32%">

<img src="figures/umap_allAXT_Prrx1.png?v=3" width="32%">


## Annotations 

![](figures/umap_allAXT_celltype.png?v=3)

![](figures/umap_allAXT_celltypeON.png?v=3)


## Cell ratios 


![](figures/allAXT_cell_ratios.png?v=1)

| celltype | Non_Regen_14_DPA | Reg | Reg_4wpa | Regen_14_DPA | Uninjured1 | Uninjured2 | nonReg | nonReg_4wpa | Total |
|---|---|---|---|---|---|---|---|---|---|
| B-Cells | 3 | 260 | 303 | 4 | 8 | 29 | 300 | 135 | 1042 |
| Endothelial | 101 | 1576 | 1477 | 511 | 484 | 1676 | 1112 | 1450 | 8387 |
| Keratinocyte | 803 | 697 | 1193 | 2 | 276 | 97 | 508 | 616 | 4192 |
| Lymphatic_Endothelial | 62 | 149 | 239 | 2 | 35 | 7 | 231 | 370 | 1095 |
| Macrophage/Osteoclast | 504 | 1931 | 4351 | 175 | 69 | 226 | 3858 | 6440 | 17554 |
| Mast Cells | 19 | 31 | 45 | 1 | 0 | 1 | 125 | 75 | 297 |
| Mesenchymal | 609 | 4300 | 4362 | 399 | 1126 | 1362 | 2741 | 6166 | 21065 |
| Neutrophil | 1479 | 259 | 1018 | 3 | 0 | 12 | 513 | 492 | 3776 |
| Osteosarcoma | 7 | 1247 | 6553 | 17 | 13 | 28 | 9009 | 6246 | 23120 |
| Pericyte/MSC | 135 | 574 | 393 | 234 | 534 | 2031 | 581 | 545 | 5027 |
| Schwann | 37 | 190 | 126 | 53 |



## Differential Gene Expression (DGE)

This workflow performs pairwise differential gene expression analysis between a target condition (**Group**) and a baseline control (**Reference**). 

To ensure the statistical comparison is strict and focused:
1. The dataset is first subsetted to isolate only the cells belonging to the target **Group** and the baseline **Reference**, excluding all unrelated sample groups.
2. Pairwise statistical testing is performed across the remaining cells using the **Wilcoxon rank-sum test** on normalized, log-transformed expression counts (`log1p`).


#### Differential Expression Logic & Directionality
The output fold changes reflect expression levels in the target group relative to the reference baseline:
* **Positive Log Fold Change (logFC > 0):** Indicates genes that are **upregulated** in the target group compared to reference.
* **Negative Log Fold Change (logFC < 0):** Indicates genes that are **downregulated** in the target group compared to reference.

`
**Filtered CSV Output:** Contains only genes that pass the significant p-value threshold.

#### Heatmap Selection & Sorting Logic
The top-ranked differential genes selected for display in the heatmap are chosen based on a two-step process prioritizing **biological effect size (magnitude of expression change)** among statistically validated genes:

1. **Filtering by Statistical Significance:** Only genes meeting the adjusted p-value threshold ($p_{\text{adj}}$ < 0.05) are eligible for heatmap inclusion.
2. **Directional Splitting:** Significant genes are split into two distinct groups:
   * Upregulated genes (logFC > 0)
   * Downregulated genes (logFC < 0)
3. **LogFC Ranking:** 
   * **Upregulated Block (Top):** Ranked strictly by **highest positive logFC** in descending order. The top N strongest upregulated genes are selected.
   * **Downregulated Block (Bottom):** Ranked strictly by **most negative logFC** in ascending order. The top N strongest downregulated genes are selected.



###### Reg vs nonReg


[Download Filtered DGE Results (Adjusted p-value < 0.05)](https://docs.google.com/spreadsheets/d/1omm4zeO5QSxG9W-GHoHs8MzxyBSJReD2iRRc1uunU2s/edit?usp=sharing)


###### Reg_4wp vs nonReg_4wpa


[Download Filtered DGE Results (Adjusted p-value < 0.05)](https://docs.google.com/spreadsheets/d/161OVUODgvmB0QYfQ0_YQVAhP0mV6h9t_UKBwZpWo8Ig/edit?usp=sharing)





