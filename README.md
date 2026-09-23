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

- Checked and organized the sample metadata.
- Loaded the single-cell expression matrices into Seurat objects.
- Calculated nFeature_RNA, nCount_RNA, and percent.mt for quality assessment.
- Applied initial cell-quality filtering based on gene and UMI counts.
- Performed doublet detection and removed predicted doublets.
- Continued with normalization, dimensionality reduction, clustering, cell-type annotation, and downstream cell–cell communication analysis.

## How to Run It

Run the scripts in the following order:

1. `01_data_QC.R` — Performs metadata checking, QC, initial filtering, and doublet detection/removal.
2. `02_normalization_PCA.R` — Performs normalization, identification of highly variable genes, scaling, and PCA.
3. `03_clustering_annotation.R` — Performs clustering, UMAP, and cell-type annotation.
4. `04_response_IL1B_CD8.R` — Performs response-group and IL1B+ monocyte / CD8+ T-cell analysis.
5. `05_CellChat_IL1B_CD8.R` — Performs ligand–receptor communication analysis between IL1B+ monocytes and CD8+ T cells.


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

### Key Figures

---

## Team

- [Aya Ayman Abdullah](https://github.com/AyaAymanAbdullah)
- [Member 2](https://github.com/USERNAME2)
- [Member 3](https://github.com/USERNAME3)
- [Member 4](https://github.com/USERNAME4)
- [Member 5](https://github.com/USERNAME5) 