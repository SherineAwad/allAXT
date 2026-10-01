#AXT project - merging different experiments


### GSE135985

<img src="figures/violin_GSE135985_preQC.png?v=1" width="48%"> <img src="figures/violin_GSE135985_AfterQC.png?v=1" width="48%">

### AXT

<img src="figures/violin_AXT_preQC.png?v=1" width="48%"> <img src="figures/violin_AXT_AfterQC.png?v=1" width="48%">

### SLX-28276

<img src="figures/violin_SLX-28276_preQC.png?v=1" width="48%"> <img src="figures/violin_SLX-28276_AfterQC.png?v=1" width="48%">

### Now doublet detection for each sample 

<img src="figures/AXT_scrublet_scores.png?v=1" width="32%"> <img src="figures/GSE135985_scrublet_scores.png?v=1" width="32%"> <img src="figures/SLX-28276_scrublet_scores.png?v=1" width="32%">

| Dataset | Cells Before | Predicted Doublets | Predicted Doublet Rate | Cells After |
|---|---:|---:|---:|---:|
| SLX-28276 | 47,188 | 921 | 2.0% | 46,267 |
| AXT | 32,593 | 446 | 1.37% | 32,147 |
| GSE135985 | 13,784 | 29 | 0.21% | 13,755 |

## Now, we did merge the 3 datasets

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

and now per dataset:

| Experiment | Cells |
|---|---:|
| SLX-28276 | 46,267 |
| AXT | 32,147 |
| GSE135985 | 13,755 |



#### Scale, normalise, and UMAP

![](figures/umap_allAXT_Regen_14_DPA.png?v=1)

##### Per sample 
<img src="figures/umap_allAXT_Uninjured1.png?v=1" width="48%"> <img src="figures/umap_allAXT_Uninjured2.png?v=1" width="48%">

<img src="figures/umap_allAXT_nonReg_4wpa.png?v=1" width="48%"> <img src="figures/umap_allAXT_Non_Regen_14_DPA.png?v=1" width="48%">

<img src="figures/umap_allAXT_nonReg.png?v=1" width="48%"> <img src="figures/umap_allAXT_Reg.png?v=1" width="48%">

<img src="figures/umap_allAXT_Reg_4wpa.png?v=1" width="48%"> <img src="figures/umap_allAXT_umap.png?v=1" width="48%">
