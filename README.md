# GSE205506_IL1B_CD8
Post-treatment single-cell RNA-seq analysis of GSE205506 to investigate IL1B⁺ monocyte–CD8⁺ T-cell ligand–receptor interactions associated with tumor response to PD-1 blockade

## The Question

Which ligand–receptor interactions between IL1B+ monocytes and CD8+ T cells are associated with differential tumor response to PD-1 blockade in dMMR/MSI-H colorectal cancer?

## The Data

The dataset was obtained from the published study:

Li J, Wu C, Hu H, et al. Remodeling of the immune and stromal cell compartment by PD-1 blockade in mismatch repair-deficient colorectal cancer. Cancer Cell. 2023.

- GEO accession: GSE205506
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205506

### Download

Download the processed single-cell data from the GSE205506 GEO page. The `GSE205506_RAW.tar` archive contains the processed gene-expression matrices for the samples. Extract the files before running the analysis.

Raw data and large intermediate files are not stored in this GitHub repository.

## What We Did

Cell-level QC was performed on 324,020 cells using paper-based gene/UMI filters, followed by scDblFinder, leaving 183,601 singlets.  
<br>
Merged the 40 samples (Seurat v5, one layer per sample) and log-normalized (scale factor 10,000).
<br>
Selected 2,000 HVGs (vst), computed in batches of layers to limit memory.
<br>
Downsampled to at most 2,060 cells per sample (from totally 183601 cells to 80158) because full scaling exceeded RAM.
<br>
Scaled the HVGs, regressed out nCount_RNA, ran PCA (20 PCs calculated, 15 used from the elbow plot).
<br>
Checked batch effects with a UMAP colored by GSM. Samples mixed, so Harmony was not applied.
<br>
Built the SNN graph (15 PCs) and clustered with Louvain at resolution 1.2 (30 clusters).
<br>
Broad cell compartments were assigned using canonical markers to support mitochondrial filtering. The dataset was reduced from 80,158 to 63,635 cells by compartment-specific mitochondrial filtering.  
<br>
The filtered dataset was then re-clustered for DE-based final annotation including cell subtypes.

## How to Run It

Run the scripts in the following order:

1. `01_data_QC.R` —  Calculates QC metrics, applies the gene/UMI filters, removes predicted doublets with scDblFinder.
2. `02_normalization_PCA`. — Performs normalization, identification of highly variable genes, downsampling, scaling, regression, PCA and clustering.
3. `03-Annotation and Mitocondrial filtration.R` - Assigns broad cell compartments based on canonical markers, providing an initial broad annotation before mitochondrial filtering.
    It then calculates compartment-specific mitochondrial thresholds and removes cells with high mitochondrial content, resulting in a filtered dataset of 63,635 cells.
    Finally, the filtered dataset is re-normalized and re-clustered, followed by identification of final DE markers and assignment of final cell-type annotations based on DE marker evidence.
4. `04-Targeted-analysis.R` — Performs ligand–receptor communication analysis between IL1B+ monocytes and CD8+ T cells.


### R

- R 4.6.1

### Packages used in the QC stage

| Package | Version |
|---|---:|
| Seurat | 5.5.1 |
| dplyr | 1.1.2 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |
| Matrix | 1.7.6 |
| SingleCellExperiment | 1.34.0 |
| scDblFinder | 1.26.7 |

## Results
Figure 1. UMAP colored by GSM showing mixed samples clusters so no need harmony integration.
<img src="figures/Norm_figures/umap_before_harmony.png" width="1637" height="849" alt="umap_GSM_noharmony">
Figure 2. Final cell-type annotation after mitochondrial filtering, showing distinct epithelial, T-cell, B-cell, myeloid, endothelial, and fibroblast/stromal subpopulations.
<img width="1637" height="849" alt="Final Cell-Type Annotation After Mitochondrial Filtering 1" src="https://github.com/user-attachments/assets/2fb41b4f-6440-48c1-9661-1d4c2534f319" />


## Team

- [Aya Ayman Abdullah](https://github.com/AyaAymanAbdullah)
- [Shahd Karam](https://github.com/ShahdKaram)
- [Member 3](https://github.com/USERNAME3)
