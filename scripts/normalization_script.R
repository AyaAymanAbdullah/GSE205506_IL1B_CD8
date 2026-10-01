# ============================================================
# GSE205506
# POST-QC ANALYSIS
# MERGE -> NORMALIZATION -> HVGs -> HVG PLOT
# -> DOWNSAMPLING -> SCALING + REGRESSION -> PCA
# -> ELBOW PLOT -> UMAP BEFORE HARMONY
# -> HARMONY -> UMAP AFTER HARMONY
# -> NEIGHBORS -> CLUSTERING
#
# START:
#   After the QC script has finished successfully.
#
# STOP:
#   Immediately before cell-type annotation and
#   compartment-specific mitochondrial filtering.
#
# NOTE:
#   The QC script produces 40 doublet-filtered RDS files.
#   These files are merged here into one Seurat object
#   while keeping the RNA expression in separate Seurat v5
#   layers.
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(harmony)
library(dplyr)

set.seed(100)

seed_use <- 100


# ============================================================
# 2. Parameters
# ============================================================

# Metadata column that identifies the original sample for each cell.
#
# In this dataset, each GSM corresponds to one sample.
# We use GSM as the BATCH VARIABLE for the integration decision.
#
# IMPORTANT: this does NOT mean Harmony is run separately on one GSM
# at a time. If Harmony is chosen below, all GSMs are supplied together
# and Harmony uses GSM to account for sample-to-sample technical/batch
# variation.
#
# Do NOT replace this with the biological condition we want to compare.
# The project guide says integration should be by donor/batch, never by
# the condition being tested.
batch_var <- "GSM"

# Number of highly variable genes
n_hvg <- 2000

# Number of PCs to calculate before inspecting the elbow plot
n_pcs_to_calculate <- 20

# Working value for downstream analysis.
#
# IMPORTANT:
# Inspect the saved PCA elbow plot first.
# Change this value if the elbow plot supports another choice.
#
# This is NOT being treated as a fixed biological truth.
n_pcs <- 15

# Clustering resolution
resolution <- 1.2

# Maximum number of cells kept from each GSM.
#
# This is a memory-management step because the full
# post-doublet dataset contains 183,601 cells and the
# full ScaleData + regression step exceeded available RAM.
#
# 1,500 cells x 40 samples = up to 60,000 cells.
max_cells_per_sample <- 1500


# ============================================================
# 3. Paths
# ============================================================

# Folder produced by the QC script
post_qc_dir <- "D:/Downloads/GSE205506_doublet_filtered"

# Optional QC summary produced by the QC script
doublet_summary_file <- "D:/Downloads/GSE205506_doublet_summary.rds"

# Output folders
results_dir <- "D:/Downloads/GSE205506_results"
fig_dir <- "D:/Downloads/GSE205506_figures"
checkpoint_dir <- "D:/Downloads/GSE205506_checkpoints"

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  checkpoint_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 4. Load the 40 post-QC samples
# ============================================================
# The raw GEO matrices are NOT loaded here.
# Only the RDS files produced by the completed QC workflow
# are loaded.
# ============================================================

post_qc_files <- sort(
  list.files(
    post_qc_dir,
    pattern = "_doublet_filtered\\.rds$",
    full.names = TRUE
  )
)

cat(
  "Post-QC files found:",
  length(post_qc_files),
  "\n"
)

stopifnot(
  length(post_qc_files) == 40
)

sample_ids <- sub(
  "_doublet_filtered\\.rds$",
  "",
  basename(post_qc_files)
)

post_qc_list <- lapply(
  post_qc_files,
  readRDS
)

names(post_qc_list) <- sample_ids

cat(
  "40 post-QC samples loaded successfully.\n"
)


# ============================================================
# 5. Merge all 40 post-QC samples
# ============================================================
# Seurat v5 keeps the merged expression data in separate
# layers, which is suitable for the later integration step.
#
# DO NOT run JoinLayers() here.
# The full JoinLayers() step previously caused memory errors.
# ============================================================

message("Merging the 40 post-QC samples...")

merged_data <- merge(
  x = post_qc_list[[1]],
  y = post_qc_list[-1],
  add.cell.ids = names(post_qc_list),
  project = "GSE205506"
)

message("Merge completed.")

# The sample list is no longer needed
rm(post_qc_list)
gc()


# ============================================================
# 6. Pre-normalization checks
# ============================================================
# Nothing biological needs to be done between QC and
# normalization other than merging and checking the input.
# Percent.mt was already calculated during QC.
# ============================================================

# Metadata column identifying the original sample
# Each GSM represents one sample and is used as the batch variable
# for the Harmony integration decision.
batch_var <- "GSM"


DefaultAssay(merged_data) <- "RNA"

cat(
  "\nMerged dataset:\n",
  "Cells: ", ncol(merged_data), "\n",
  "Genes: ", nrow(merged_data), "\n",
  sep = ""
)

stopifnot(
  batch_var %in% colnames(merged_data[[]])
)

sample_count <- length(
  unique(
    merged_data[[batch_var]][, 1]
  )
)

cat(
  "Samples represented:",
  sample_count,
  "\n"
)

stopifnot(
  sample_count == 40
)

rna_layers_before_norm <- Layers(
  merged_data[["RNA"]]
)

cat("\nRNA layers before normalization:\n")
print(rna_layers_before_norm)

counts_layers <- Layers(
  merged_data[["RNA"]],
  search = "^counts"
)

stopifnot(
  length(counts_layers) == 40
)

cat(
  "Counts layers:",
  length(counts_layers),
  "\n"
)

# Optional cell-count checkpoint from the QC summary
if (file.exists(doublet_summary_file)) {
  
  doublet_summary <- readRDS(
    doublet_summary_file
  )
  
  expected_cells <- sum(
    doublet_summary$cells_after_doublet
  )
  
  cat(
    "Expected cells from QC summary:",
    expected_cells,
    "\n"
  )
  
  stopifnot(
    ncol(merged_data) == expected_cells
  )
  
  rm(
    doublet_summary,
    expected_cells
  )
}


# ============================================================
# 7. Normalization
# ============================================================
# LogNormalize with scale factor 10,000.
#
# Seurat v5 processes the separate RNA layers without requiring
# JoinLayers() first.
# ============================================================

message("Starting LogNormalize...")

normalized_data <- NormalizeData(
  merged_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

rm(merged_data)
gc()

message("LogNormalize completed.")

data_layers <- Layers(
  normalized_data[["RNA"]],
  search = "^data"
)

cat(
  "Normalized data layers:",
  length(data_layers),
  "\n"
)

stopifnot(
  length(data_layers) == 40
)


# ============================================================
# 8. Highly Variable Genes
# ============================================================
# The full dataset is large, so HVGs are calculated in batches
# of 10 layers to reduce peak RAM usage.
#
# Method:
#   VST
#   2,000 highly variable genes
# ============================================================

# Number of highly variable genes to select
n_hvg <- 2000

message("Finding highly variable genes...")

hvg_batch_size <- 10

hvg_layer_batches <- split(
  counts_layers,
  ceiling(
    seq_along(counts_layers) /
      hvg_batch_size
  )
)

for (i in seq_along(hvg_layer_batches)) {
  
  current_hvg_layers <- hvg_layer_batches[[i]]
  
  message(
    "HVG batch ",
    i,
    " / ",
    length(hvg_layer_batches),
    ": ",
    paste(
      current_hvg_layers,
      collapse = ", "
    )
  )
  
  normalized_data <- FindVariableFeatures(
    normalized_data,
    assay = "RNA",
    layer = current_hvg_layers,
    selection.method = "vst",
    nfeatures = n_hvg,
    verbose = FALSE
  )
}


# ============================================================
# 9. Store the final 2,000 HVGs
# ============================================================

hvg_genes <- VariableFeatures(
  normalized_data[["RNA"]]
)

cat(
  "\nFinal number of HVGs:",
  length(hvg_genes),
  "\n"
)

stopifnot(
  length(hvg_genes) == n_hvg
)

cat("\nTop 20 HVGs:\n")
print(
  head(
    hvg_genes,
    20
  )
)


# ============================================================
# 10. HVG visualization
# ============================================================
# VariableFeaturePlot() gave a layered-object error in the
# current environment, so the equivalent HVF information is
# plotted directly from HVFInfo().
#
# X-axis:
#   Mean expression
#
# Y-axis:
#   Standardized variance
#
# Selected HVGs are highlighted.
# ============================================================

hvg_info <- HVFInfo(
  normalized_data[["RNA"]],
  method = "vst"
)

hvg_info$gene <- rownames(
  hvg_info
)

hvg_info$variable <- (
  hvg_info$gene %in% hvg_genes
)

hvg_info$mean_plot <- pmax(
  hvg_info$mean,
  1e-6
)

top10_hvg <- head(
  hvg_genes,
  10
)

p_hvg <- ggplot(
  hvg_info,
  aes(
    x = mean_plot,
    y = variance.standardized,
    color = variable
  )
) +
  geom_point(
    size = 0.5
  ) +
  scale_x_log10() +
  ggrepel::geom_text_repel(
    data = subset(
      hvg_info,
      gene %in% top10_hvg
    ),
    aes(
      label = gene
    ),
    size = 3,
    max.overlaps = 10,
    show.legend = FALSE
  ) +
  labs(
    title = "Highly Variable Genes",
    x = "Average Expression",
    y = "Standardized Variance",
    color = "HVG"
  ) +
  theme_classic()

print(p_hvg)

ggsave(
  file.path(
    fig_dir,
    "variable_features_labeled.pdf"
  ),
  p_hvg,
  width = 7,
  height = 5
)

write.csv(
  data.frame(
    rank = seq_along(hvg_genes),
    gene = hvg_genes
  ),
  file.path(
    results_dir,
    "GSE205506_HVGs_2000.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. Downsample cells within each GSM
# ============================================================
# This is a computational memory-management step.
#
# Every GSM remains represented.
# No GSM contributes more than max_cells_per_sample cells.
# ============================================================

# Maximum number of cells retained from each GSM
max_cells_per_sample <- 2060



cells_before_downsampling <- ncol(
  normalized_data
)

cells_by_sample <- split(
  colnames(normalized_data),
  normalized_data[[batch_var]][, 1]
)

cells_to_keep <- unlist(
  lapply(
    cells_by_sample,
    function(cell_ids) {
      
      if (
        length(cell_ids) >
        max_cells_per_sample
      ) {
        
        sample(
          cell_ids,
          size = max_cells_per_sample,
          replace = FALSE
        )
        
      } else {
        
        cell_ids
        
      }
    }
  ),
  use.names = FALSE
)

downsampled_data <- subset(
  normalized_data,
  cells = cells_to_keep
)

rm(
  normalized_data,
  cells_by_sample,
  cells_to_keep
)

gc()

cells_after_downsampling <- ncol(
  downsampled_data
)

cat(
  "\nCells before downsampling:",
  cells_before_downsampling,
  "\n"
)

cat(
  "Cells after downsampling:",
  cells_after_downsampling,
  "\n"
)

stopifnot(
  length(
    unique(
      downsampled_data[[batch_var]][, 1]
    )
  ) == 40
)


# ============================================================
# 12. Restore the final 2,000 HVG labels after subsetting
# ============================================================

VariableFeatures(
  downsampled_data[["RNA"]]
) <- hvg_genes

stopifnot(
  length(
    VariableFeatures(
      downsampled_data[["RNA"]]
    )
  ) == n_hvg
)


# ============================================================
# 13. Scale data + regress nCount_RNA
# ============================================================
# Scaling and regression are performed in the same step.
#
# features:
#   2,000 HVGs
#
# vars.to.regress:
#   nCount_RNA = total UMI count per cell
#
# block.size is reduced to limit peak RAM usage.
# ============================================================

message(
  "Scaling data and regressing nCount_RNA..."
)

scaled_data <- ScaleData(
  downsampled_data,
  assay = "RNA",
  features = hvg_genes,
  vars.to.regress = "nCount_RNA",
  block.size = 100,
  verbose = TRUE
)

rm(downsampled_data)
gc()

message("ScaleData completed.")


# ============================================================
# 14. PCA
# ============================================================
# Calculate 20 PCs first.
# The final number used downstream is selected after the
# elbow plot.
# ============================================================
n_pcs_to_calculate <- 20
message("Running PCA...")

pca_data <- RunPCA(
  scaled_data,
  assay = "RNA",
  features = hvg_genes,
  npcs = n_pcs_to_calculate,
  seed.use = seed_use,
  verbose = FALSE
)

rm(scaled_data)
gc()

message("PCA completed.")

stopifnot(
  "pca" %in% Reductions(pca_data)
)





# ============================================================
# 15. PCA visualization
# ============================================================

p_pca <- DimPlot(
  pca_data,
  reduction = "pca"
) +
  ggtitle(
    "PCA"
  )

print(p_pca)

ggsave(
  file.path(
    fig_dir,
    "pca.pdf"
  ),
  p_pca,
  width = 7,
  height = 5
)


# ============================================================
# 16. PCA Elbow Plot
# ============================================================
# This is the decision point for the number of PCs.
#
# The script currently has n_pcs = 15 as a working value.
# Inspect this plot and update n_pcs above if the elbow
# supports a different number.
# ============================================================

p_elbow <- ElbowPlot(
  pca_data,
  ndims = n_pcs_to_calculate,
  reduction = "pca"
)

print(p_elbow)



n_pcs <- 15

ggsave(
  file.path(
    fig_dir,
    "pca_elbow.pdf"
  ),
  p_elbow,
  width = 6,
  height = 4
)

cat(
  "\nWorking n_pcs value:",
  n_pcs,
  "\n"
)

cat(
  "Inspect figures/pca_elbow.pdf before finalizing n_pcs.\n"
)


# ============================================================
# 17. PCA heatmap
# ============================================================

# ============================================================
# 17. PCA heatmap
# ============================================================

dims_for_heatmap <- min(
  n_pcs,
  n_pcs_to_calculate
)

pdf(
  file.path(
    fig_dir,
    "pca_heatmap.pdf"
  ),
  width = 8,
  height = 10
)

print(
  DimHeatmap(
    pca_data,
    dims = 1:dims_for_heatmap,
    cells = 500,
    balanced = TRUE
  )
)

dev.off()


# ============================================================
# 18. UMAP BEFORE Harmony
# ============================================================
# This is used to inspect sample/batch structure before
# integration.
# ============================================================

message(
  "Running UMAP before Harmony..."
)

umap_before_harmony <- RunUMAP(
  pca_data,
  reduction = "pca",
  dims = 1:n_pcs,
  reduction.name = "umap_before",
  seed.use = seed_use,
  verbose = FALSE
)

message(
  "UMAP before Harmony completed."
)

stopifnot(
  "umap_before" %in%
    Reductions(umap_before_harmony)
)

#visualize

DimHeatmap(
  pca_data,
  dims = 1:n_pcs,
  cells = 500,
  balanced = TRUE
)
# ============================================================
# 19. Plot UMAP before Harmony by GSM
# ============================================================

p_before <- DimPlot(
  umap_before_harmony,
  reduction = "umap_before",
  group.by = batch_var
) +
  NoLegend() +
  ggtitle(
    "Before Harmony Integration"
  )

print(p_before)

ggsave(
  file.path(
    fig_dir,
    "umap_before_harmony.pdf"
  ),
  p_before,
  width = 8,
  height = 6
)


# ============================================================
# 20. INTEGRATION DECISION — FOLLOW THE PROJECT GUIDE
# ============================================================
# The guide says:
#   1) Colour the UMAP by sample FIRST.
#   2) If samples separate, integrate.
#   3) If samples overlap/mix, do NOT integrate.
#
# We have already created and saved:
#   figures/umap_before_harmony.pdf
#
# The script asks you to inspect that figure in RStudio and enter:
#   y = samples are clearly separated -> run Harmony
#   n = samples overlap/mix -> do NOT run Harmony
#
# This is intentionally a manual decision because Harmony does not
# automatically decide whether integration is biologically/methodologically
# necessary.
# ============================================================

integration_answer <- tolower(
  trimws(
    readline(
      paste0(
        "\nInspect figures/umap_before_harmony.pdf. ",
        "Do the GSM samples clearly separate? [y/n]: "
      )
    )
  )
)

if (!(integration_answer %in% c("y", "n"))) {
  stop(
    "Please enter only 'y' or 'n' for the integration decision."
  )
}

use_harmony <- integration_answer == "y"

cat(
  "\nHarmony selected:",
  use_harmony,
  "\n"
)

# ============================================================
# 21. No Harmony: Use PCA Representation
# ============================================================

message(
  "Samples overlap/mix in the pre-Harmony UMAP. Harmony will NOT be run."
)

analysis_data <- umap_before_harmony
analysis_reduction <- "pca"

p_no_harmony <- DimPlot(
  analysis_data,
  reduction = "umap_before",
  group.by = batch_var
) +
  NoLegend() +
  ggtitle(
    "No Harmony: Samples Overlap"
  )

print(p_no_harmony)

ggsave(
  file.path(
    fig_dir,
    "umap_analysis_no_harmony.pdf"
  ),
  p_no_harmony,
  width = 8,
  height = 6
)


# ============================================================
# 24. Neighbour graph
# ============================================================
# Use Harmony embeddings when integration was selected; otherwise use
# the original PCA space.
# ============================================================

message(
  "Finding cell neighbors using ",
  analysis_reduction,
  "..."
)

neighbor_data <- FindNeighbors(
  analysis_data,
  reduction = analysis_reduction,
  dims = 1:n_pcs,
  verbose = FALSE
)

rm(analysis_data)
gc()

message(
  "Neighbor graph completed."
)


# ============================================================
# 25. Clustering
# ============================================================
# First-round clustering.
#
# Annotation is NOT performed here.
# ============================================================
message(
  "Finding clusters..."
)
resolution <- 1.2

clustered_data <- FindClusters(
  neighbor_data,
  resolution = resolution,
  random.seed = seed_use,
  verbose = FALSE
)


rm(neighbor_data)
gc()

message(
  "Clustering completed."
)


# ============================================================
# 26. Check cluster assignments
# ============================================================

cluster_sizes <- as.data.frame(
  table(
    clustered_data$seurat_clusters
  )
)

colnames(cluster_sizes) <- c(
  "cluster",
  "cells"
)

print(
  cluster_sizes
)

stopifnot(
  !anyNA(
    clustered_data$seurat_clusters
  )
)

write.csv(
  cluster_sizes,
  file.path(
    results_dir,
    "GSE205506_cluster_sizes.csv"
  ),
  row.names = FALSE
)

cat(
  "\nNumber of clusters:",
  nrow(cluster_sizes),
  "\n"
)


# ============================================================
# 27. Cluster UMAP
# ============================================================

cluster_umap_reduction <- if (use_harmony) {
  "umap_after"
} else {
  "umap_before"
}

cluster_umap_title <- if (use_harmony) {
  "Clusters after Harmony Integration"
} else {
  "Clusters without Harmony"
}

p_clusters <- DimPlot(
  clustered_data,
  reduction = cluster_umap_reduction,
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle(
    cluster_umap_title
  )

print(p_clusters)

ggsave(
  file.path(
    fig_dir,
    "umap_clusters.pdf"
  ),
  p_clusters,
  width = 8,
  height = 6.5
)


# ============================================================
# 28. Save final post-clustering checkpoint
# ============================================================
# Large RDS files should NOT be committed to GitHub.
# Keep this checkpoint in local/external storage.
# ============================================================

final_checkpoint_file <- file.path(
  checkpoint_dir,
  if (use_harmony) {
    "GSE205506_postQC_harmony_clustered.rds"
  } else {
    "GSE205506_postQC_no_harmony_clustered.rds"
  }
)

saveRDS(
  clustered_data,
  final_checkpoint_file
)

message(
  "Final post-clustering checkpoint saved to: ",
  final_checkpoint_file
)


# ============================================================
# 29. Save clustered cell metadata
# ============================================================

clustered_cell_metadata <- clustered_data[[]]

clustered_cell_metadata$cell <- rownames(
  clustered_cell_metadata
)

write.csv(
  clustered_cell_metadata,
  file.path(
    results_dir,
    "GSE205506_cell_metadata_clustered.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 30. Analysis summary
# ============================================================

analysis_summary <- data.frame(
  dataset = "GSE205506",
  samples = 40,
  cells_before_downsampling =
    cells_before_downsampling,
  cells_after_downsampling =
    cells_after_downsampling,
  genes = nrow(clustered_data),
  highly_variable_genes =
    length(
      VariableFeatures(
        clustered_data[["RNA"]]
      )
    ),
  normalization =
    "LogNormalize, scale factor 10000",
  hvg_method =
    "VST",
  hvg_number =
    n_hvg,
  regression_variable =
    "nCount_RNA",
  pcs_calculated =
    n_pcs_to_calculate,
  pcs_used =
    n_pcs,
  harmony_used =
    use_harmony,
  harmony_batch_variable =
    if (use_harmony) batch_var else NA_character_,
  harmony_method =
    if (use_harmony) "Harmony" else "Not used",
  clustering_resolution =
    resolution,
  number_of_clusters =
    length(
      unique(
        clustered_data$seurat_clusters
      )
    ),
  max_cells_per_sample =
    max_cells_per_sample,
  stringsAsFactors = FALSE
)

print(
  analysis_summary
)

write.csv(
  analysis_summary,
  file.path(
    results_dir,
    "GSE205506_analysis_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 31. Save software versions
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "GSE205506_sessionInfo_clustering.txt"
  )
)


# ============================================================
# 32. Final checks
# ============================================================

check_samples <- (
  length(
    unique(
      clustered_data[[batch_var]][, 1]
    )
  ) == 40
)

check_hvg <- (
  length(
    VariableFeatures(
      clustered_data[["RNA"]]
    )
  ) == n_hvg
)

check_pca <- (
  "pca" %in%
    Reductions(clustered_data)
)

check_harmony <- (
  !use_harmony ||
    "harmony" %in% Reductions(clustered_data)
)

check_umap_before <- (
  "umap_before" %in%
    Reductions(clustered_data)
)

check_umap_after <- (
  !use_harmony ||
    "umap_after" %in% Reductions(clustered_data)
)

check_clusters <- (
  !anyNA(
    clustered_data$seurat_clusters
  )
)

cat("\n========================================\n")
cat("FINAL POST-QC CHECKPOINTS\n")
cat("========================================\n")

cat(
  "40 samples represented:",
  check_samples,
  "\n"
)

cat(
  "2,000 HVGs:",
  check_hvg,
  "\n"
)

cat(
  "PCA present:",
  check_pca,
  "\n"
)

cat(
  "Harmony used:",
  use_harmony,
  "\n"
)

cat(
  "Harmony decision/representation check:",
  check_harmony,
  "\n"
)

cat(
  "UMAP before Harmony present:",
  check_umap_before,
  "\n"
)

cat(
  "Post-integration UMAP present when Harmony is used:",
  check_umap_after,
  "\n"
)

cat(
  "No missing cluster assignments:",
  check_clusters,
  "\n"
)

stopifnot(
  check_samples,
  check_hvg,
  check_pca,
  check_harmony,
  check_umap_before,
  check_umap_after,
  check_clusters
)

cat(
  "\nALL POST-QC CHECKPOINTS PASSED.\n"
)


# ============================================================
# STOP HERE
# ============================================================
# Next stage:
#
#   1. Broad cell-type annotation
#   2. Canonical marker + cluster DE evidence
#   3. Compartment-specific mitochondrial filtering
#   4. Re-clustering after the appropriate filtering
#
# Integration note:
#   Harmony was used only if the UMAP-by-GSM inspection indicated that
#   sample integration was needed.
#
# Do NOT continue into annotation in this script.
# ============================================================
