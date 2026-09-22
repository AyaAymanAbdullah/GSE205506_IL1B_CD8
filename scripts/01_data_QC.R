# GSE205506 - Stage 1 & 2
# Data loading, metadata extraction, QC, and doublet detection/removal
#
# This script:
# 1. Retrieves GEO sample metadata
# 2. Calculates cell-level QC metrics
# 3. Applies the Pass 1 gene/UMI filters used in the original paper
# 4. Visualizes QC before and after Pass 1
# 5. Detects and removes predicted doublets using scDblFinder
# 6. Saves QC and doublet checkpoints for downstream analysis
#
# Note:
# Raw data and RDS files are NOT committed to GitHub.
# For a fresh run, place the downloaded GSE205506_RAW folder under:
#   data/GSE205506_RAW


# ============================================================
# 1. Packages
# ============================================================

library(Seurat)
library(GEOquery)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)


# ============================================================
# 2. Paths
# ============================================================

# Raw GEO matrix files
data_path <- file.path("data", "GSE205506_RAW")

# Output folders
results_dir <- file.path("results")
filtered_dir <- file.path(results_dir, "GSE205506_filtered")
doublet_dir <- file.path(results_dir, "GSE205506_doublet_filtered")

# Create output folders if they do not exist
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doublet_dir, recursive = TRUE, showWarnings = FALSE)

# Check that the raw data folder exists
stopifnot(dir.exists(data_path))


# ============================================================
# 3. Get metadata from GEO
# ============================================================

gse <- getGEO(
  "GSE205506",
  GSEMatrix = FALSE,
  destdir = results_dir
)

gsm_list <- GSMList(gse)

metadata_geo <- do.call(
  rbind,
  lapply(gsm_list, function(x) {
    data.frame(
      GSM = Meta(x)$geo_accession,
      title = Meta(x)$title,
      characteristics = paste(
        Meta(x)$characteristics_ch1,
        collapse = " | "
      ),
      stringsAsFactors = FALSE
    )
  })
)


# ============================================================
# 4. Match metadata to the 40 matrix files
# ============================================================

matrix_files <- list.files(
  data_path,
  pattern = "_matrix\\.mtx\\.gz$"
)

sample_info <- data.frame(
  matrix_file = matrix_files,
  stringsAsFactors = FALSE
)

sample_info$GSM <- sub(
  "_.*",
  "",
  sample_info$matrix_file
)

metadata <- merge(
  sample_info,
  metadata_geo,
  by = "GSM",
  all.x = TRUE,
  sort = FALSE
)

metadata <- metadata[
  match(sample_info$GSM, metadata$GSM),
]

rownames(metadata) <- NULL


# Extract useful sample-level information from GEO characteristics
metadata$patient <- sub(
  ".*subject: ([^|]+).*",
  "\\1",
  metadata$characteristics
)

metadata$sample_status <- sub(
  ".*genotype: ([^|]+).*",
  "\\1",
  metadata$characteristics
)

metadata$anatomical_site <- sub(
  ".*tissue: ([^|]+).*",
  "\\1",
  metadata$characteristics
)

metadata$treatment <- sub(
  ".*treatment: ([^|]+).*",
  "\\1",
  metadata$characteristics
)

metadata$cell_type_GEO <- sub(
  ".*cell type: ([^|]+).*",
  "\\1",
  metadata$characteristics
)

metadata$patient <- trimws(metadata$patient)
metadata$sample_status <- trimws(metadata$sample_status)
metadata$anatomical_site <- trimws(metadata$anatomical_site)
metadata$treatment <- trimws(metadata$treatment)
metadata$cell_type_GEO <- trimws(metadata$cell_type_GEO)

# Basic metadata checks
stopifnot(nrow(metadata) == 40)

print(table(metadata$sample_status))
print(table(metadata$treatment))
print(sum(is.na(metadata$patient)))
print(sum(is.na(metadata$sample_status)))
print(sum(is.na(metadata$treatment)))


# ============================================================
# 5. Collect cell-level QC metrics before filtering
# ============================================================

qc_cells_list <- vector("list", nrow(metadata))

for (i in seq_len(nrow(metadata))) {

  gsm <- metadata$GSM[i]

  message(
    "Collecting cell-level QC: ",
    i, " / ", nrow(metadata),
    " - ", gsm
  )

  matrix_file <- list.files(
    data_path,
    pattern = paste0(
      "^", gsm, ".*_matrix\\.mtx\\.gz$"
    ),
    full.names = TRUE
  )

  features_file <- list.files(
    data_path,
    pattern = paste0(
      "^", gsm, ".*_features\\.tsv\\.gz$"
    ),
    full.names = TRUE
  )

  barcodes_file <- list.files(
    data_path,
    pattern = paste0(
      "^", gsm, ".*_barcodes\\.tsv\\.gz$"
    ),
    full.names = TRUE
  )

  stopifnot(
    length(matrix_file) == 1,
    length(features_file) == 1,
    length(barcodes_file) == 1
  )

  counts <- ReadMtx(
    mtx = matrix_file,
    features = features_file,
    cells = barcodes_file,
    feature.column = 2,
    cell.column = 1
  )

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE205506",
    min.cells = 0,
    min.features = 0
  )

  # Calculate mitochondrial RNA percentage
  obj$percent.mt <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  qc_cells_list[[i]] <- data.frame(
    cell = colnames(obj),
    GSM = gsm,
    patient = metadata$patient[i],
    sample_status = metadata$sample_status[i],
    treatment = metadata$treatment[i],
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA = obj$nCount_RNA,
    percent.mt = obj$percent.mt,
    stringsAsFactors = FALSE
  )

  # Free memory before loading the next sample
  rm(counts, obj)
  gc()
}

qc_cells <- do.call(
  rbind,
  qc_cells_list
)

rownames(qc_cells) <- NULL

print(dim(qc_cells))


# ============================================================
# 6. QC visualization before Pass 1 filtering
# ============================================================

ggplot(qc_cells, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  labs(
    title = "nFeature_RNA across all cells",
    x = "nFeature_RNA",
    y = "Number of cells"
  )

ggplot(qc_cells, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  labs(
    title = "nCount_RNA across all cells",
    x = "nCount_RNA",
    y = "Number of cells"
  )

ggplot(qc_cells, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  labs(
    title = "Mitochondrial percentage across all cells",
    x = "percent.mt",
    y = "Number of cells"
  )

# Violin plots by sample
ggplot(qc_cells, aes(x = GSM, y = nFeature_RNA)) +
  geom_violin(fill = "grey70", scale = "width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "nFeature_RNA before Pass 1 filtering",
    x = "Sample",
    y = "nFeature_RNA"
  )

ggplot(qc_cells, aes(x = GSM, y = nCount_RNA)) +
  geom_violin(fill = "grey70", scale = "width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "nCount_RNA before Pass 1 filtering",
    x = "Sample",
    y = "nCount_RNA"
  )

ggplot(qc_cells, aes(x = GSM, y = percent.mt)) +
  geom_violin(fill = "grey70", scale = "width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "percent.mt before Pass 1 filtering",
    x = "Sample",
    y = "percent.mt"
  )


# ============================================================
# 7. Pass 1 filtering
# ============================================================
# Thresholds are based on the original paper:
#   nFeature_RNA > 500 and < 5000
#   nCount_RNA   > 400 and < 25000
#
# Mitochondrial filtering is intentionally NOT applied here.
# The original paper handled mitochondrial filtering by cell
# compartment after initial clustering/broad cell identification.

qc_filter_summary <- vector(
  "list",
  nrow(metadata)
)

for (i in seq_len(nrow(metadata))) {

  gsm <- metadata$GSM[i]

  message(
    "Pass 1 filtering: ",
    i, " / ", nrow(metadata),
    " - ", gsm
  )

  matrix_file <- list.files(
    data_path,
    pattern = paste0(
      "^", gsm, ".*_matrix\\.mtx\\.gz$"
    ),
    full.names = TRUE
  )

  features_file <- list.files(
    data_path,
    pattern = paste0(
      "^", gsm, ".*_features\\.tsv\\.gz$"
    ),
    full.names = TRUE
  )

  barcodes_file <- list.files(
    data_path,
    pattern = paste0(
      "^", gsm, ".*_barcodes\\.tsv\\.gz$"
    ),
    full.names = TRUE
  )

  stopifnot(
    length(matrix_file) == 1,
    length(features_file) == 1,
    length(barcodes_file) == 1
  )

  counts <- ReadMtx(
    mtx = matrix_file,
    features = features_file,
    cells = barcodes_file,
    feature.column = 2,
    cell.column = 1
  )

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE205506",
    min.cells = 0,
    min.features = 0
  )

  obj$percent.mt <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  cells_before <- ncol(obj)

  # Apply Pass 1 gene/UMI thresholds
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 5000 &
      nCount_RNA > 400 &
      nCount_RNA < 25000
  )

  cells_after <- ncol(obj)

  # Add sample metadata
  obj$GSM <- gsm
  obj$patient <- metadata$patient[i]
  obj$sample_status <- metadata$sample_status[i]
  obj$anatomical_site <- metadata$anatomical_site[i]
  obj$treatment <- metadata$treatment[i]

  # Make cell names unique across samples
  obj <- RenameCells(
    obj,
    add.cell.id = gsm
  )

  # Save Pass 1 filtered object
  saveRDS(
    obj,
    file.path(
      filtered_dir,
      paste0(gsm, "_filtered.rds")
    )
  )

  qc_filter_summary[[i]] <- data.frame(
    GSM = gsm,
    patient = metadata$patient[i],
    sample_status = metadata$sample_status[i],
    treatment = metadata$treatment[i],
    cells_before = cells_before,
    cells_after = cells_after,
    cells_removed = cells_before - cells_after,
    percent_retained =
      (cells_after / cells_before) * 100,
    stringsAsFactors = FALSE
  )

  rm(counts, obj)
  gc()
}

qc_filter_summary <- do.call(
  rbind,
  qc_filter_summary
)

rownames(qc_filter_summary) <- NULL

# Overall cell counts after Pass 1
cat(
  "Cells started with =", sum(qc_filter_summary$cells_before), "\n"
)

cat(
  "Cells after Pass 1 =", sum(qc_filter_summary$cells_after), "\n"
)

cat(
  "Cells removed in Pass 1 =", sum(qc_filter_summary$cells_removed), "\n"
)


# ============================================================
# 8. QC visualization after Pass 1 filtering
# ============================================================

qc_after_list <- vector(
  "list",
  nrow(metadata)
)

for (i in seq_len(nrow(metadata))) {

  gsm <- metadata$GSM[i]

  message(
    "Collecting post-filter QC: ",
    i, " / ", nrow(metadata),
    " - ", gsm
  )

  obj <- readRDS(
    file.path(
      filtered_dir,
      paste0(gsm, "_filtered.rds")
    )
  )

  qc_after_list[[i]] <- data.frame(
    cell = colnames(obj),
    GSM = gsm,
    patient = metadata$patient[i],
    sample_status = metadata$sample_status[i],
    treatment = metadata$treatment[i],
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA = obj$nCount_RNA,
    percent.mt = obj$percent.mt,
    stringsAsFactors = FALSE
  )

  rm(obj)
  gc()
}

qc_after <- do.call(
  rbind,
  qc_after_list
)

rownames(qc_after) <- NULL

ggplot(qc_after, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  labs(
    title = "nFeature_RNA after Pass 1 filtering",
    x = "nFeature_RNA",
    y = "Number of cells"
  )

ggplot(qc_after, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  labs(
    title = "nCount_RNA after Pass 1 filtering",
    x = "nCount_RNA",
    y = "Number of cells"
  )

ggplot(qc_after, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  labs(
    title = "Mitochondrial percentage after Pass 1 filtering",
    x = "percent.mt",
    y = "Number of cells"
  )

# Violin plots by sample
ggplot(qc_after, aes(x = GSM, y = nFeature_RNA)) +
  geom_violin(fill = "grey70", scale = "width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "nFeature_RNA after Pass 1 filtering",
    x = "Sample",
    y = "nFeature_RNA"
  )

ggplot(qc_after, aes(x = GSM, y = nCount_RNA)) +
  geom_violin(fill = "grey70", scale = "width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "nCount_RNA after Pass 1 filtering",
    x = "Sample",
    y = "nCount_RNA"
  )

ggplot(qc_after, aes(x = GSM, y = percent.mt)) +
  geom_violin(fill = "grey70", scale = "width") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "percent.mt after Pass 1 filtering",
    x = "Sample",
    y = "percent.mt"
  )


# ============================================================
# 9. Check sample-level cell retention
# ============================================================

# Identify samples with less than 50% cell retention.
# These samples are documented here but are NOT removed automatically.

low_retention_samples <- qc_filter_summary[
  qc_filter_summary$percent_retained < 50,
  c(
    "GSM",
    "patient",
    "cells_before",
    "cells_after",
    "percent_retained"
  )
]

print(low_retention_samples)


# ============================================================
# 10. Doublet detection and removal with scDblFinder
# ============================================================

# scDblFinder identifies cells predicted to contain two or more
# cells captured together.

set.seed(100)

doublet_summary <- vector(
  "list",
  nrow(metadata)
)

for (i in seq_len(nrow(metadata))) {

  gsm <- metadata$GSM[i]

  message(
    "Doublet detection: ",
    i, " / ", nrow(metadata),
    " - ", gsm
  )

  # Load the Pass 1 filtered object
  obj <- readRDS(
    file.path(
      filtered_dir,
      paste0(gsm, "_filtered.rds")
    )
  )

  cells_before_doublet <- ncol(obj)

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(obj)

  # Run scDblFinder
  sce <- scDblFinder(sce)

  # Add scDblFinder results to the Seurat object
  obj$doublet_score <- colData(sce)$scDblFinder.score
  obj$doublet_class <- colData(sce)$scDblFinder.class

  # Count predicted doublets
  n_doublets <- sum(
    obj$doublet_class == "doublet"
  )

  # Count singlets before saving
  n_singlets <- sum(
    obj$doublet_class == "singlet"
  )

  # Remove predicted doublets
  obj <- subset(
    obj,
    subset = doublet_class == "singlet"
  )

  cells_after_doublet <- ncol(obj)

  # Save object after doublet removal
  saveRDS(
    obj,
    file.path(
      doublet_dir,
      paste0(gsm, "_doublet_filtered.rds")
    )
  )

  doublet_summary[[i]] <- data.frame(
    GSM = gsm,
    cells_before_doublet = cells_before_doublet,
    predicted_doublets = n_doublets,
    singlets_after_doublet = n_singlets,
    cells_after_doublet = cells_after_doublet,
    stringsAsFactors = FALSE
  )

  rm(obj, sce)
  gc()
}

doublet_summary <- do.call(
  rbind,
  doublet_summary
)

rownames(doublet_summary) <- NULL


# ============================================================
# 11. Final doublet-removal checkpoint
# ============================================================

# Verify that all 40 post-doublet objects were created
doublet_files <- list.files(
  doublet_dir,
  pattern = "_doublet_filtered\\.rds$",
  full.names = TRUE
)

stopifnot(
  length(doublet_files) == nrow(metadata)
)

# Calculate totals reported by scDblFinder
total_singlets <- sum(doublet_summary$singlets_after_doublet)
total_predicted_doublets <- sum(doublet_summary$predicted_doublets)

cat(
  "Total singlets after removal =",
  total_singlets,
  "\n"
)

cat(
  "Total predicted doublets removed =",
  total_predicted_doublets,
  "\n"
)


# Verify that no predicted doublets remain in the saved objects
total_remaining_doublets <- 0

for (file in doublet_files) {

  obj <- readRDS(file)

  total_remaining_doublets <- total_remaining_doublets +
    sum(obj$doublet_class == "doublet")

  rm(obj)
  gc()
}

cat(
  "Total doublets remaining after removal =",
  total_remaining_doublets,
  "\n"
)


# ============================================================
# 12. Save QC and analysis checkpoints
# ============================================================

# Save metadata
saveRDS(
  metadata,
  file.path(results_dir, "GSE205506_metadata.rds")
)

write.csv(
  metadata,
  file.path(results_dir, "GSE205506_metadata.csv"),
  row.names = FALSE
)

# Save QC table before Pass 1
saveRDS(
  qc_cells,
  file.path(results_dir, "qc_cells.rds")
)

# Save QC table after Pass 1
saveRDS(
  qc_after,
  file.path(results_dir, "qc_after.rds")
)

write.csv(
  qc_after,
  file.path(results_dir, "qc_after.csv"),
  row.names = FALSE
)

# Save Pass 1 filtering summary
saveRDS(
  qc_filter_summary,
  file.path(results_dir, "qc_filter_summary.rds")
)

write.csv(
  qc_filter_summary,
  file.path(results_dir, "qc_filter_summary.csv"),
  row.names = FALSE
)

# Save doublet summary
saveRDS(
  doublet_summary,
  file.path(results_dir, "GSE205506_doublet_summary.rds")
)

write.csv(
  doublet_summary,
  file.path(results_dir, "GSE205506_doublet_summary.csv"),
  row.names = FALSE
)

# Save final doublet checkpoint
doublet_removal_checkpoint <- data.frame(
  total_singlets_after_removal = total_singlets,
  total_predicted_doublets_removed = total_predicted_doublets,
  total_doublets_remaining = total_remaining_doublets
)

saveRDS(
  doublet_removal_checkpoint,
  file.path(
    results_dir,
    "GSE205506_doublet_removal_checkpoint.rds"
  )
)

write.csv(
  doublet_removal_checkpoint,
  file.path(
    results_dir,
    "GSE205506_doublet_removal_checkpoint.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. Record software versions
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "GSE205506_sessionInfo.txt"
  )
)

message("QC and doublet-removal workflow completed.")
