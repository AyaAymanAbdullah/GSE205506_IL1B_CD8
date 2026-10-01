library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
clustered_data <- readRDS(
  "D:/Downloads/GSE205506_before_JoinLayers.rds"
)



dim(clustered_data)


table(clustered_data$seurat_clusters)

clustered_data <- JoinLayers(clustered_data)

# ============================================================
# Check the RNA layers after JoinLayers
# ============================================================

Layers(
  clustered_data[["RNA"]]
)




# ============================================================
# Check that the Seurat v5 RNA layers were successfully joined
# ============================================================
Layers(clustered_data[["RNA"]])
# ============================================================
# Find marker genes for every cluster
# This identifies genes that are enriched in each cluster
# compared with the remaining cells.
# Parameters match the tutorial and the paper:
# only.pos = TRUE
# min.pct = 0.25
# logfc.threshold = 0.25
DefaultAssay(clustered_data) <- "RNA"

markers <- FindAllMarkers(
  clustered_data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# ============================================================
# Save the complete cluster marker table
# ============================================================

write.csv(
  markers,
  "D:/Downloads/GSE205506_cluster_markers.csv",
  row.names = FALSE
)

# ============================================================
# Select the top 10 marker genes for each cluster
# The genes are ranked by average log2 fold-change.
# ============================================================

library(dplyr)

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()

# ============================================================
# Save the top 10 markers for each cluster
# ============================================================

write.csv(
  top10_markers,
  "D:/Downloads/GSE205506_top10_markers_per_cluster.csv",
  row.names = FALSE
)
# ============================================================
# Count how many top marker genes were selected for each cluster
# ============================================================

table(
  top10_markers$cluster
)
# ============================================================
# Select the top 5 DE marker genes for each cluster
# This follows the tutorial's DEG-based annotation workflow.
# ============================================================

top5_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  ) %>%
  ungroup()

# ============================================================
# Check the number of top DE markers selected per cluster
# ============================================================

table(
  top5_markers$cluster
)
# ============================================================
# Save the top 5 DE marker genes for each cluster
# ============================================================

write.csv(
  top5_markers,
  "D:/Downloads/GSE205506_top5_DE_markers_per_cluster.csv",
  row.names = FALSE
)
# ============================================================
# Broad canonical marker panel for first-round annotation
# These are established markers for the six broad cell
# compartments in the GSE205506 dataset.
# No subtype markers are included at this stage.
# ============================================================

broad_canonical_markers <- c(
  
  # T cells
  "CD3D", "CD3E", "TRAC", "TRBC1",
  
  # B cells
  "MS4A1", "CD79A", "CD79B", "CD19",
  
  # Myeloid cells
  "LYZ", "CD14", "CTSS", "CD68",
  
  # Epithelial cells
  "EPCAM", "KRT8", "KRT18", "KRT19",
  
  # Endothelial cells
  "VWF", "PECAM1", "CLDN5",
  
  # Fibroblast / stromal cells
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"
)

# Keep only genes that are present in our dataset
broad_canonical_markers <- broad_canonical_markers[
  broad_canonical_markers %in%
    rownames(clustered_data[["RNA"]])
]

# Create the broad canonical DotPlot
p_broad_canonical_dotplot <- DotPlot(
  clustered_data,
  features = broad_canonical_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle(
    "Broad Canonical Marker Expression Across Clusters"
  )

print(p_broad_canonical_dotplot)


# ============================================================
# Broad T-cell markers Violin Plot
# ============================================================

p_t_violin <- VlnPlot(
  clustered_data,
  features = c("CD3D", "CD3E", "TRAC", "TRBC1"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

print(p_t_violin)
# ============================================================
# Broad B-cell markers Violin Plot
# ============================================================

p_b_violin <- VlnPlot(
  clustered_data,
  features = c("MS4A1", "CD79A", "CD79B", "CD19"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

print(p_b_violin)
# ============================================================
# Broad Myeloid markers Violin Plot
# ============================================================

p_myeloid_violin <- VlnPlot(
  clustered_data,
  features = c("LYZ", "CD14", "CTSS", "CD68"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

print(p_myeloid_violin)
# ============================================================
# Broad Epithelial markers Violin Plot
# ============================================================

p_epithelial_violin <- VlnPlot(
  clustered_data,
  features = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

print(p_epithelial_violin)
# ============================================================
# Broad Endothelial markers Violin Plot
# ============================================================

p_endothelial_violin <- VlnPlot(
  clustered_data,
  features = c("VWF", "PECAM1", "CLDN5"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

print(p_endothelial_violin)
# ============================================================
# Broad Fibroblast / Stromal markers Violin Plot
# ============================================================

p_fibroblast_violin <- VlnPlot(
  clustered_data,
  features = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM"
  ),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)

# ============================================================
# Split the 24 broad canonical markers into groups
# of 5 markers per FeaturePlot
# ============================================================

feature_markers_1 <- broad_canonical_markers[1:5]
feature_markers_2 <- broad_canonical_markers[6:10]
feature_markers_3 <- broad_canonical_markers[11:15]
feature_markers_4 <- broad_canonical_markers[16:20]
feature_markers_5 <- broad_canonical_markers[21:24]
# ============================================================
# Broad canonical Feature Plot - Group 1
# ============================================================

p_broad_featureplot_1 <- FeaturePlot(
  clustered_data,
  features = feature_markers_1,
  reduction = "umap",
  ncol = 5,
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  patchwork::plot_annotation(
    title = "Broad Canonical Markers - Group 1"
  )

print(p_broad_featureplot_1)
# ============================================================
# Broad canonical Feature Plot - Group 2
# ============================================================

p_broad_featureplot_2 <- FeaturePlot(
  clustered_data,
  features = feature_markers_2,
  reduction = "umap",
  ncol = 5,
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  patchwork::plot_annotation(
    title = "Broad Canonical Markers - Group 2"
  )

print(p_broad_featureplot_2)
# ============================================================
# Broad canonical Feature Plot - Group 3
# ============================================================

p_broad_featureplot_3 <- FeaturePlot(
  clustered_data,
  features = feature_markers_3,
  reduction = "umap",
  ncol = 5,
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  patchwork::plot_annotation(
    title = "Broad Canonical Markers - Group 3"
  )

print(p_broad_featureplot_3)
# ============================================================
# Broad canonical Feature Plot - Group 4
# ============================================================

p_broad_featureplot_4 <- FeaturePlot(
  clustered_data,
  features = feature_markers_4,
  reduction = "umap",
  ncol = 5,
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  patchwork::plot_annotation(
    title = "Broad Canonical Markers - Group 4"
  )

print(p_broad_featureplot_4)
# ============================================================
# Broad canonical Feature Plot - Group 5
# ============================================================

p_broad_featureplot_5 <- FeaturePlot(
  clustered_data,
  features = feature_markers_5,
  reduction = "umap",
  ncol = 5,
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  patchwork::plot_annotation(
    title = "Broad Canonical Markers - Group 5"
  )

print(p_broad_featureplot_5)
# ============================================================
# Canonical marker-based broad cell-type annotation
# ============================================================

canonical_annotation <- c(
  `0`  = "T cells",
  `1`  = "Epithelial cells",
  `2`  = "B cells",
  `3`  = "Epithelial cells",
  `4`  = "Epithelial cells",
  `5`  = "T cells",
  `6`  = "Epithelial cells",
  `7`  = "Epithelial cells",
  `8`  = "Epithelial cells",
  `9`  = "Myeloid cells",
  `10` = "Epithelial cells",
  `11` = "Endothelial cells",
  `12` = "Epithelial cells",
  `13` = "Endothelial cells",
  `14` = "Epithelial cells",
  `15` = "Epithelial cells",
  `16` = "Epithelial cells",
  `17` = "Epithelial cells",
  `18` = "Epithelial cells",
  `19` = "Epithelial cells",
  `20` = "Epithelial cells",
  `21` = "Endothelial cells",
  `22` = "Epithelial cells",
  `23` = "Epithelial cells",
  `24` = "Fibroblast cells",
  `25` = "B cells",
  `26` = "Epithelial cells",
  `27` = "Fibroblast cells",
  `28` = "T cells"
)

# ============================================================
# Assign canonical annotation to cell metadata
# ============================================================

clustered_data$canonical_annotation <- unname(
  canonical_annotation[
    as.character(clustered_data$seurat_clusters)
  ]
)

# Check that every cell received an annotation
stopifnot(
  !anyNA(clustered_data$canonical_annotation)
)

# Check annotation counts
print(
  table(clustered_data$canonical_annotation)
)
# ============================================================
# UMAP colored by canonical broad annotation
# ============================================================

p_canonical_umap <- DimPlot(
  clustered_data,
  reduction = "umap",
  group.by = "canonical_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Canonical Marker-Based Broad Cell-Type Annotation"
  )

print(p_canonical_umap)

ggsave(
  "GSE205506_canonical_broad_annotation_umap.pdf",
  p_canonical_umap,
  width = 10,
  height = 7
)
# ============================================================
# SAVE CHECKPOINT AFTER CANONICAL ANNOTATION
# ============================================================
# This checkpoint allows the analysis to continue directly
# from the canonical-annotation stage without re-running:
# JoinLayers, FindAllMarkers, canonical marker plots,
# or canonical annotation.
# ============================================================

canonical_checkpoint_file <- 
  "D:/Downloads/GSE205506_canonical_annotated.rds"

saveRDS(
  clustered_data,
  canonical_checkpoint_file
)

message(
  "\nCanonical annotation checkpoint saved to:\n",
  canonical_checkpoint_file
)
# ============================================================
# SAVE DE MARKER OBJECTS AS RDS
# ============================================================
# These objects will be needed for the independent
# DE-based annotation stage.
# ============================================================

saveRDS(
  markers,
  "D:/Downloads/GSE205506_cluster_markers.rds"
)

saveRDS(
  top10_markers,
  "D:/Downloads/GSE205506_top10_markers_per_cluster.rds"
)

message(
  "\nDE marker objects saved successfully.\n"
)
# ============================================================
# SAVE CANONICAL ANNOTATION INFORMATION
# ============================================================

canonical_annotation_objects <- list(
  broad_canonical_markers = broad_canonical_markers,
  feature_markers_1 = feature_markers_1,
  feature_markers_2 = feature_markers_2,
  feature_markers_3 = feature_markers_3,
  feature_markers_4 = feature_markers_4,
  feature_markers_5 = feature_markers_5,
  canonical_annotation = canonical_annotation
)

saveRDS(
  canonical_annotation_objects,
  "D:/Downloads/GSE205506_canonical_annotation_objects.rds"
)

message(
  "\nCanonical annotation objects saved successfully.\n"
)


# ============================================================
# Check the broad canonical annotation
# ============================================================

table(
  clustered_data$canonical_annotation
)
# ============================================================
# Calculate mitochondrial RNA percentage
# ============================================================

clustered_data[["percent.mt"]] <- PercentageFeatureSet(
  clustered_data,
  pattern = "^MT-"
)

summary(
  clustered_data$percent.mt
)
# ============================================================
# Inspect mitochondrial content by broad cell compartment
# ============================================================

VlnPlot(
  clustered_data,
  features = "percent.mt",
  group.by = "canonical_annotation",
  pt.size = 0
)
# ============================================================
#centered MAD-variance normal distribution
# followed by Bonferroni correction, as described in the paper

library(dplyr)

mad_bonferroni_filter <- function(x, alpha = 0.05) {
  
  med <- median(x, na.rm = TRUE)
  
  # Robust standard deviation estimated from MAD
  mad_sd <- mad(
    x,
    center = med,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  # Standardized values
  z <- (x - med) / mad_sd
  
  # One-sided upper-tail p-values
  p_value <- pnorm(
    z,
    lower.tail = FALSE
  )
  
  # Bonferroni correction within this compartment
  p_adj <- p.adjust(
    p_value,
    method = "bonferroni"
  )
  
  # Cells considered unusually high
  remove <- p_adj < alpha
  
  # Cutoff = highest value that is not significantly high
  if (any(!remove)) {
    threshold <- max(x[!remove], na.rm = TRUE)
  } else {
    threshold <- med
  }
  
  list(
    median = med,
    mad_sd = mad_sd,
    threshold = threshold,
    p_value = p_value,
    p_adj = p_adj,
    remove = remove
  )
}
# ============================================================
### Apply mitochondrial filtering separately to each broad compartment
# using the canonical_annotation metadata column directly
# ============================================================
compartments <- unique(clustered_data$canonical_annotation)

mito_results <- list()

for (comp in compartments) {
  
  # Get cell names belonging to this compartment directly from metadata
  cells <- rownames(clustered_data@meta.data)[
    clustered_data@meta.data$canonical_annotation == comp
  ]
  
  x <- clustered_data$percent.mt[cells]
  
  # Epithelial cells: fixed 75% cutoff from the paper
  if (comp == "Epithelial cells") {
    
    threshold <- 75
    remove <- x > threshold
    
    mito_results[[comp]] <- list(
      median = median(x, na.rm = TRUE),
      mad_sd = NA,
      threshold = threshold,
      p_value = NA,
      p_adj = NA,
      remove = remove
    )
    
  } else {
    
    mito_results[[comp]] <- mad_bonferroni_filter(x)
  }
}
# Check that all compartments were processed successfully

names(mito_results)
# ============================================================
# Build mitochondrial filtering summary
# ============================================================
mito_summary <- data.frame(
  compartment = character(),
  cells_before = integer(),
  median_percent_mt = numeric(),
  mad_sd = numeric(),
  threshold_percent_mt = numeric(),
  cells_to_remove = integer(),
  cells_after = integer(),
  percent_removed = numeric(),
  stringsAsFactors = FALSE
)

for (comp in compartments) {
  
  cells <- rownames(clustered_data@meta.data)[
    clustered_data@meta.data$canonical_annotation == comp
  ]
  
  result <- mito_results[[comp]]
  
  cells_before <- length(cells)
  cells_to_remove <- sum(result$remove)
  cells_after <- cells_before - cells_to_remove
  
  mito_summary <- rbind(
    mito_summary,
    data.frame(
      compartment = comp,
      cells_before = cells_before,
      median_percent_mt = result$median,
      mad_sd = result$mad_sd,
      threshold_percent_mt = result$threshold,
      cells_to_remove = cells_to_remove,
      cells_after = cells_after,
      percent_removed = 100 * cells_to_remove / cells_before
    )
  )
}

mito_summary
# Count cells that will be removed according to the statistical rule
# and the fixed 75% epithelial cutoff

sum(sapply(mito_results, function(x) sum(x$remove)))
# Save mitochondrial filtering results before removing any cells

saveRDS(
  mito_results,
  "D:/Downloads/GSE205506_mito_filtering_results.rds"
)

write.csv(
  mito_summary,
  "D:/Downloads/GSE205506_mito_filtering_summary.csv",
  row.names = FALSE
)
# Save the current object before mitochondrial filtering
# so the original 80,158-cell dataset remains 

saveRDS(
  clustered_data,
  "D:/Downloads/GSE205506_before_mito_filtering.rds"
)
# ============================================================
# IDENTIFY CELLS TO REMOVE BASED ON MITOCHONDRIAL FILTERING
# ============================================================

cells_to_remove <- unlist(
  lapply(
    mito_results,
    function(x) {
      names(x$remove)[x$remove]
    }
  ),
  use.names = FALSE
)

# Number of cells to remove
length(cells_to_remove)

# Check that all selected cells exist in the object
all(cells_to_remove %in% colnames(clustered_data))

# Check for duplicates
length(unique(cells_to_remove))

## ============================================================
# CREATE THE LIST OF CELLS TO KEEP
# ============================================================

cells_to_keep <- setdiff(
  colnames(clustered_data),
  cells_to_remove
)

length(cells_to_keep)

# ============================================================
# ============================================================
# CHECK CURRENT R MEMORY USAGE
# ============================================================

gc()
# ============================================================
# CLEAN LARGE OBJECTS BEFORE MITOCHONDRIAL SUBSETTING
# ============================================================

# Remove large analysis objects that are no longer needed
rm(
  markers,
  top10_markers,
  top5_markers,
  broad_canonical_markers
)

gc()

# ============================================================
# PREPARE A LIGHTER OBJECT FOR MITOCHONDRIAL FILTERING
# ============================================================

filtered_data <- clustered_data

# Remove old dimensional reductions
filtered_data[["pca"]] <- NULL
filtered_data[["umap"]] <- NULL

# Remove old graphs and neighbors
filtered_data@graphs <- list()
filtered_data@neighbors <- list()

gc()



# ============================================================
# CHECK SEURAT OBJECT AND RNA LAYERS
# ============================================================

packageVersion("Seurat")
packageVersion("SeuratObject")

Layers(clustered_data[["RNA"]])

# ============================================================
# CREATE A LIGHT OBJECT FOR MITOCHONDRIAL FILTERING
# KEEP COUNTS + DATA ONLY
# ============================================================

filtered_data <- DietSeurat(
  clustered_data,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE,
  commands = FALSE
)

gc()

# Check the remaining RNA layers
Layers(filtered_data[["RNA"]])

# ============================================================
# APPLY MITOCHONDRIAL FILTERING
# KEEP ONLY THE 63,635 VALID CELLS
# ============================================================

filtered_data <- filtered_data[, cells_to_keep]

# Check number of cells after filtering
ncol(filtered_data)

# Check that the expected RNA layers are still present
Layers(filtered_data[["RNA"]])

# ============================================================
# SAVE CHECKPOINT AFTER MITOCHONDRIAL FILTERING
# 63,635 CELLS
# ============================================================

saveRDS(
  filtered_data,
  "D:/Downloads/GSE205506_after_mito_filtering_63635cells.rds"
)
# ============================================================
# VERIFY FINAL FILTERED OBJECT BEFORE RE-CLUSTERING
# ============================================================

ncol(filtered_data)

dim(filtered_data)

nrow(filtered_data[["RNA"]])

# ============================================================
# VERIFY CELL COUNT IN FINAL METADATA
# ============================================================

nrow(filtered_data@meta.data)



# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 1 — FINAL NORMALIZATION
# ============================================================

DefaultAssay(filtered_data) <- "RNA"

filtered_data <- NormalizeData(
  filtered_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

cat(
  "Final normalization completed.\n"
)

# ============================================================
# END — FINAL NORMALIZATION
# ============================================================


## ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 2 — FINAL HIGHLY VARIABLE GENES
# ============================================================

filtered_data <- FindVariableFeatures(
  filtered_data,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

final_variable_features <- VariableFeatures(
  filtered_data
)

cat(
  "Number of final highly variable genes:",
  length(final_variable_features),
  "\n"
)

# ============================================================
# END — FINAL HIGHLY VARIABLE GENES
# ============================================================


# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 3 — FINAL SCALING
# ============================================================

filtered_data <- ScaleData(
  filtered_data,
  features = VariableFeatures(filtered_data),
  vars.to.regress = "nCount_RNA",
  block.size = 200,
  verbose = TRUE
)

cat(
  "Final scaling completed.\n"
)

# ============================================================
# END — FINAL SCALING
# ============================================================

# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 4 — FINAL PCA
# ============================================================

final_npcs_to_calculate <- 20

filtered_data <- RunPCA(
  filtered_data,
  assay = "RNA",
  features = VariableFeatures(filtered_data),
  npcs = final_npcs_to_calculate,
  reduction.name = "final_pca",
  verbose = TRUE
)

cat(
  "Final PCA completed.\n"
)

print(
  Reductions(filtered_data)
)

# ============================================================
# END — FINAL PCA
# ============================================================

# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 5 — FINAL PCA ELBOW PLOT
# ============================================================

final_elbow_plot <- ElbowPlot(
  filtered_data,
  reduction = "final_pca",
  ndims = final_npcs_to_calculate
)

print(final_elbow_plot)

# ============================================================
# END — FINAL PCA ELBOW PLOT
# ============================================================


# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 6 — SELECT FINAL NUMBER OF PCs
# ============================================================

final_n_pcs <- 10

cat(
  "Final number of PCs used downstream:",
  final_n_pcs,
  "\n"
)

# ============================================================
# END — FINAL PC SELECTION
# ============================================================

# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 7 — FINAL NEIGHBOR GRAPH
# ============================================================

filtered_data <- FindNeighbors(
  filtered_data,
  reduction = "final_pca",
  dims = 1:final_n_pcs,
  verbose = TRUE
)

cat(
  "Final neighbor graph completed using PCs 1-",
  final_n_pcs,
  ".\n",
  sep = ""
)

# ============================================================
# END — FINAL NEIGHBOR GRAPH
# ============================================================


# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 8 — FINAL CLUSTERING
# ============================================================

final_resolution <- 1.2

filtered_data <- FindClusters(
  filtered_data,
  resolution = final_resolution,
  random.seed = 100,
  verbose = TRUE
)

# Store final cluster identities separately
filtered_data$final_clusters <- Idents(filtered_data)

cat(
  "Final clustering completed.\n"
)

cat(
  "Number of final clusters:",
  length(unique(filtered_data$final_clusters)),
  "\n"
)

# ============================================================
# END — FINAL CLUSTERING
# ============================================================


# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 9 — FINAL UMAP
# ============================================================

filtered_data <- RunUMAP(
  filtered_data,
  reduction = "final_pca",
  dims = 1:final_n_pcs,
  reduction.name = "final_umap",
  seed.use = 100,
  verbose = TRUE
)

cat(
  "Final UMAP completed.\n"
)

# ============================================================
# END — FINAL UMAP
# ============================================================



# ============================================================
# FINAL INTEGRATION ASSESSMENT
# STEP 10 — FINAL UMAP BY SAMPLE (GSM)
# ============================================================

final_umap_by_gsm <- DimPlot(
  filtered_data,
  reduction = "final_umap",
  group.by = "GSM"
) +
  NoLegend() +
  ggtitle(
    "Final UMAP by Sample (GSM) — Integration Assessment"
  )

print(final_umap_by_gsm)

ggsave(
  "D:/Downloads/GSE205506_final_umap_by_GSM.pdf",
  final_umap_by_gsm,
  width = 8,
  height = 6
)

# ============================================================
# END — FINAL INTEGRATION ASSESSMENT
# ============================================================



# ============================================================
# FINAL RE-CLUSTERING AFTER MITOCHONDRIAL FILTERING
# STEP 11 — FINAL CLUSTER UMAP
# NO HARMONY
# ============================================================

final_cluster_umap <- DimPlot(
  filtered_data,
  reduction = "final_umap",
  group.by = "final_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Final Clusters After Mitochondrial Filtering"
  )

print(final_cluster_umap)


# ============================================================
# END — FINAL CLUSTER UMAP
# ============================================================


# ============================================================
# FINAL ANNOTATION AFTER MITOCHONDRIAL FILTERING
# STEP 12 — FINAL DIFFERENTIAL EXPRESSION MARKERS
# ============================================================

DefaultAssay(filtered_data) <- "RNA"

final_markers <- FindAllMarkers(
  filtered_data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cat(
  "Final DE marker analysis completed.\n"
)

cat(
  "Number of marker rows:",
  nrow(final_markers),
  "\n"
)

# ============================================================
# FINAL ANNOTATION AFTER MITOCHONDRIAL FILTERING
# STEP 13 — TOP 10 FINAL MARKERS PER CLUSTER
# ============================================================

final_top10_markers <- final_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()

print(final_top10_markers)

# ============================================================
# END — TOP 10 FINAL MARKERS
# ============================================================

# ============================================================
# END — FINAL DE MARKERS
# ============================================================


# ============================================================
# SAVE FINAL DE MARKER RESULTS
# ============================================================

saveRDS(
  final_markers,
  "D:/Downloads/GSE205506_final_DE_markers.rds"
)

saveRDS(
  final_top10_markers,
  "D:/Downloads/GSE205506_final_top10_markers_per_cluster.rds"
)

write.csv(
  final_markers,
  "D:/Downloads/GSE205506_final_DE_markers.csv",
  row.names = FALSE
)

write.csv(
  final_top10_markers,
  "D:/Downloads/GSE205506_final_top10_markers_per_cluster.csv",
  row.names = FALSE
)

# ============================================================
# END — SAVE FINAL DE MARKERS
# ============================================================



# ============================================================
# FINAL ANNOTATION AFTER MITOCHONDRIAL FILTERING
# STEP 12 — INSPECT TOP DE MARKERS BY CLUSTER
# ============================================================

final_top10_markers %>%
  dplyr::select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  print(
    n = Inf
  )

# ============================================================
# END — TOP DE MARKERS
# ============================================================

# ============================================================
# FINAL DE-BASED ANNOTATION
# STEP 12 — EXTRACT SUBTYPE EVIDENCE FROM ALL FINAL MARKERS
# ============================================================

# Marker panels based on the biological subtypes relevant
# to this project and the published GSE205506 analysis.

final_marker_panels <- list(
  
  # ----------------------------------------------------------
  # T / CD8-related subtypes
  # ----------------------------------------------------------
  
  "T cells" = c(
    "CD3D", "CD3E", "TRAC", "TRBC1"
  ),
  
  "CD8 T cells" = c(
    "CD8A", "CD8B"
  ),
  
  "CD4 T cells" = c(
    "CD4", "IL7R", "LTB", "CCR7", "MAL", "TCF7"
  ),
  
  "CD8 Tem" = c(
    "GZMK", "CMC1", "CST7", "KLRG1"
  ),
  
  "CD8 Trm" = c(
    "ENTPD1", "ITGAE", "PDCD1",
    "CXCL13", "GZMB", "GNLY",
    "HAVCR2", "TNFRSF18"
  ),
  
  "CD8 Trm-mitotic" = c(
    "MKI67", "TOP2A", "STMN1",
    "PDCD1", "HMGB2", "PCLAF"
  ),
  
  "CD8 IEL" = c(
    "NR4A1", "NR4A2", "HOPX",
    "CD69", "CD160", "KLRB1"
  ),
  
  "MAIT" = c(
    "KLRB1", "TRAV1-2", "SLC4A10",
    "KLRD1", "ZBTB16"
  ),
  
  "CD4 Th" = c(
    "CD40LG", "IL7R"
  ),
  
  "CD4 Treg" = c(
    "FOXP3", "IL2RA", "CTLA4",
    "TNFRSF18", "TNFRSF4", "TIGIT"
  ),
  
  "Cytotoxic lymphocytes" = c(
    "NKG7", "GNLY", "PRF1",
    "GZMA", "GZMB", "GZMH",
    "CTSW", "XCL1", "XCL2"
  ),
  
  # ----------------------------------------------------------
  # B-cell subtypes
  # ----------------------------------------------------------
  
  "B cells" = c(
    "MS4A1", "CD79A", "CD79B", "CD22"
  ),
  
  "Naive B cells" = c(
    "TCL1A", "IGHD", "FCER2",
    "IL4R", "CR2"
  ),
  
  "Germinal-center B cells" = c(
    "AICDA", "RGS13", "SERPINA9",
    "GCSAM", "BCL6", "CXCR4", "CD69"
  ),
  
  "Plasma cells" = c(
    "JCHAIN", "MZB1", "XBP1",
    "TNFRSF17", "IGHA1", "IGHA2",
    "IGHG1", "IGHG2", "IGHG3", "IGHG4",
    "IGKC", "IGLC2", "IGLC3"
  ),
  
  # ----------------------------------------------------------
  # Myeloid subtypes
  # ----------------------------------------------------------
  
  "IL1B+ inflammatory monocytes" = c(
    "IL1B", "FCN1",
    "S100A8", "S100A9", "S100A12",
    "CXCL8", "CXCL2", "CXCL1",
    "CCL2", "CCL4"
  ),
  
  "Macrophages" = c(
    "C1QA", "C1QB", "C1QC",
    "TREM2", "FOLR2", "APOC1",
    "APOE", "MS4A4A", "VSIG4"
  ),
  
  "DC cells" = c(
    "FCER1A", "CD1C", "CLEC10A",
    "CST3", "HLA-DRA", "HLA-DPA1"
  ),
  
  # ----------------------------------------------------------
  # Endothelial / stromal
  # ----------------------------------------------------------
  
  "ACKR1+ endothelial" = c(
    "ACKR1", "SELP", "CCL14", "AQP1"
  ),
  
  "Endothelial cells" = c(
    "VWF", "PECAM1", "CLDN5",
    "KDR", "EMCN", "MMRN1"
  ),
  
  "COL15A1+ endothelial" = c(
    "COL15A1", "COL4A1", "KDR",
    "PCDH12", "IGFBP3"
  ),
  
  "Fibroblasts" = c(
    "COL1A1", "COL1A2", "COL3A1",
    "COL6A3", "PDGFRA", "CXCL14",
    "DCN", "LUM"
  ),
  
  "PI16/SFRP2+ fibroblasts" = c(
    "PI16", "SFRP1", "SFRP2",
    "MFAP5", "SCARA5", "ADH1B"
  ),
  
  "Myofibroblasts" = c(
    "MYH11", "ACTA2", "CNN1",
    "LMOD1", "RERGL", "PLN"
  ),
  
  # ----------------------------------------------------------
  # Epithelial subtypes
  # ----------------------------------------------------------
  
  "Enterocyte-like epithelial" = c(
    "ADH1C", "CFTR", "SLC9A2",
    "UGT2A3", "HMGCS2"
  ),
  
  "Goblet cells" = c(
    "MUC2", "SPINK4", "FCGBP",
    "TFF3", "REG4"
  ),
  
  "BEST4+ epithelial" = c(
    "BEST4", "AQP8", "SLC26A3",
    "SLC17A4", "OTOP2", "TMIGD1"
  ),
  
  "OLFM4+ epithelial" = c(
    "OLFM4", "ASCL2", "HMGCS2",
    "ADH1C", "CFTR"
  ),
  
  "Tuft cells" = c(
    "POU2F3", "TRPM5", "GNG13",
    "HTR3E", "HPGDS"
  ),
  
  # ----------------------------------------------------------
  # Proliferating cells
  # ----------------------------------------------------------
  
  "Cycling cells" = c(
    "MKI67", "TOP2A", "STMN1",
    "BIRC5", "UBE2C", "NUSAP1",
    "ASPM", "CENPF", "NEK2",
    "PBK", "CDC20"
  )
)


# ------------------------------------------------------------
# Convert marker panels to a dataframe
# ------------------------------------------------------------

marker_panel_df <- do.call(
  rbind,
  lapply(
    names(final_marker_panels),
    function(cell_type) {
      data.frame(
        cell_type = cell_type,
        gene = final_marker_panels[[cell_type]],
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------
# Search ALL final DE markers
# ------------------------------------------------------------

all_final_marker_evidence <- final_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.58
  ) %>%
  inner_join(
    marker_panel_df,
    by = "gene"
  ) %>%
  arrange(
    cluster,
    cell_type,
    desc(avg_log2FC)
  )


# ------------------------------------------------------------
# Summarize evidence per cluster and subtype
# ------------------------------------------------------------

final_subtype_evidence <- all_final_marker_evidence %>%
  group_by(
    cluster,
    cell_type
  ) %>%
  summarise(
    n_matching_markers = n(),
    markers = paste(
      gene,
      collapse = ", "
    ),
    mean_log2FC = mean(avg_log2FC),
    max_log2FC = max(avg_log2FC),
    max_pct1 = max(pct.1),
    .groups = "drop"
  ) %>%
  arrange(
    cluster,
    desc(n_matching_markers),
    desc(mean_log2FC)
  )


# ------------------------------------------------------------
# Display all subtype evidence
# ------------------------------------------------------------

print(
  final_subtype_evidence,
  n = Inf
)

# ============================================================
# END — ALL-MARKER SUBTYPE EVIDENCE
# ============================================================


# ============================================================
# FINAL DE-BASED ANNOTATION
# STEP 14 — FINAL SUBTYPE ANNOTATION MAP
# BASED ON ALL FINAL DE MARKER EVIDENCE
# ============================================================

final_de_annotation_map <- c(
  "0"  = "Enterocyte-like epithelial",
  "1"  = "GZMK+ CD8 Tem-like T cells",
  "2"  = "IL7R+ CD4 T cells",
  "3"  = "Enterocyte-like epithelial",
  "4"  = "B cells",
  "5"  = "Goblet cells",
  "6"  = "ACKR1+ endothelial cells",
  "7"  = "Goblet cells",
  "8"  = "Cytotoxic CD8 T cells (IEL/Trm-like)",
  "9"  = "Cycling cells — lineage unresolved",
  "10" = "Naive B cells",
  "11" = "Enterocyte-like epithelial",
  "12" = "Endothelial cells",
  "13" = "BEST4+ epithelial cells",
  "14" = "IL1B+ inflammatory monocytes",
  "15" = "Cycling epithelial cells",
  "16" = "OLFM4+ / enterocyte-like epithelial",
  "17" = "BEST4+ epithelial cells",
  "18" = "Plasma cells",
  "19" = "Myofibroblasts",
  "20" = "Cycling epithelial cells",
  "21" = "COL15A1+ endothelial cells",
  "22" = "Goblet cells",
  "23" = "Fibroblast cells",
  "24" = "C1Q+ macrophages",
  "25" = "PI16/SFRP2+ fibroblasts",
  "26" = "Tuft cells",
  "27" = "Cycling germinal-center B cells",
  "28" = "Cycling cytotoxic CD8 T cells / Trm-mitotic-like",
  "29" = "Plasma cells",
  "30" = "Germinal-center B-like cells"
)

filtered_data$final_de_annotation <- unname(
  final_de_annotation_map[
    as.character(filtered_data$final_clusters)
  ]
)

# Check annotation for all clusters
annotation_check <- data.frame(
  cluster = names(final_de_annotation_map),
  annotation = unname(final_de_annotation_map),
  stringsAsFactors = FALSE
)

print(annotation_check)

# Check missing annotations
cat(
  "Missing final annotations:",
  sum(is.na(filtered_data$final_de_annotation)),
  "\n"
)

# ============================================================
# END — FINAL DE-BASED ANNOTATION
# ============================================================


# ============================================================
# FINAL DE-BASED ANNOTATION
# STEP 15 — ANNOTATION CHECK
# ============================================================

table(
  filtered_data$final_de_annotation
)

cat(
  "Missing final annotations:",
  sum(is.na(filtered_data$final_de_annotation)),
  "\n"
)

# ============================================================
# END — ANNOTATION CHECK
# ============================================================

# ============================================================
# FINAL ANNOTATION
# STEP 16 — FINAL ANNOTATED UMAP
# ============================================================

final_annotated_umap <- DimPlot(
  filtered_data,
  reduction = "final_umap",
  group.by = "final_de_annotation",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Final Cell-Type Annotation After Mitochondrial Filtering"
  )

print(final_annotated_umap)

ggsave(
  "D:/Downloads/GSE205506_final_annotated_umap.pdf",
  final_annotated_umap,
  width = 11,
  height = 8
)

# ============================================================
# END — FINAL ANNOTATED UMAP
# ============================================================



# ============================================================
# SAVE FINAL RE-CLUSTERING + FINAL DE-BASED ANNOTATION
# CHECKPOINT
# ============================================================

saveRDS(
  filtered_data,
  "D:/Downloads/GSE205506_final_reclustering_annotated_63635cells.rds"
)

# Save the cluster-to-annotation mapping separately
write.csv(
  annotation_check,
  "D:/Downloads/GSE205506_final_cluster_annotation.csv",
  row.names = FALSE
)

cat(
  "Final annotated checkpoint saved successfully.\n",
  "Cells: ",
  ncol(filtered_data),
  "\n",
  sep = ""
)

# ============================================================
# END — SAVE FINAL ANNOTATED CHECKPOINT
# ============================================================