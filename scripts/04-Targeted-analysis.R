# ============================================================
# GSE205506 — TARGETED IL1B+ MONOCYTE / CD8+ T-CELL ANALYSIS
# ============================================================
#
# Main question:
# Which ligand–receptor interactions between IL1B+ monocytes
# and CD8+ T cells are associated with differential tumor
# response to PD-1 blockade in dMMR/MSI-H colorectal cancer?
#
# Pipeline:
# Stage 3 checkpoint
#      ↓
# Mitochondrial QC
#      ↓
# Broad cell identification
#      ↓
# IL1B+ monocytes + CD8+ T cells
#      ↓
# pCR / non-pCR
#      ↓
# CellChat
#      ↓
# Ligand–receptor interactions
#      ↓
# Self-signaling
#      ↓
# Patient-level summaries
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)

# CellChat is loaded later after the Seurat/QC steps.


# ============================================================
# 1. DIRECTORIES
# ============================================================

results_dir <- "C:/Users/Admin/Documents/results"

fig_dir <- file.path(
  results_dir,
  "GSE205506_figures"
)

analysis_dir <- file.path(
  results_dir,
  "GSE205506_targeted_analysis"
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  analysis_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. FIND STAGE 3 SEURAT CHECKPOINT
# ============================================================

rds_files <- list.files(
  results_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

cat("\nRDS files found:\n")
print(rds_files)


# Search for likely Stage 3 objects
stage3_candidates <- rds_files[
  grepl(
    "GSE205506|cluster|before|stage|integrat",
    basename(rds_files),
    ignore.case = TRUE
  )
]

cat("\nPotential Stage 3 files:\n")
print(stage3_candidates)


# ============================================================
# 3. AUTOMATICALLY IDENTIFY SEURAT OBJECTS
# ============================================================

seurat_candidates <- list()

for (f in stage3_candidates) {
  
  obj_test <- tryCatch(
    readRDS(f),
    error = function(e) NULL
  )
  
  if (inherits(obj_test, "Seurat")) {
    
    seurat_candidates[[f]] <- obj_test
    
  }
}

cat("\nSeurat objects among candidates:\n")
print(names(seurat_candidates))


# ============================================================
# SAFETY CHECK
# ============================================================

if (length(seurat_candidates) == 0) {
  
  stop(
    paste0(
      "\nNo Seurat checkpoint was automatically identified.\n",
      "Do not continue.\n",
      "Check the filenames printed above."
    )
  )
  
}


# ============================================================
# 4. DISPLAY CANDIDATES
# ============================================================

cat("\nCandidate Seurat objects:\n\n")

for (nm in names(seurat_candidates)) {
  
  obj <- seurat_candidates[[nm]]
  
  cat(
    "\nFILE: ",
    nm,
    "\nCells: ",
    ncol(obj),
    "\nGenes: ",
    nrow(obj),
    "\nMetadata columns: ",
    paste(
      colnames(obj@meta.data),
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# IMPORTANT
# ============================================================
#
# If there is only ONE candidate, use it automatically.
#
# If there are several candidates, the script stops here so
# that an old object cannot accidentally be analysed.
# ============================================================

if (length(seurat_candidates) == 1) {
  
  stage3_file <- names(seurat_candidates)[1]
  
  filtered_data <- seurat_candidates[[1]]
  
  cat(
    "\nStage 3 object selected automatically:\n",
    stage3_file,
    "\n"
  )
  
} else {
  
  stop(
    paste0(
      "\nMore than one Seurat candidate was found.\n",
      "Choose the correct Stage 3 checkpoint before continuing."
    )
  )
  
}


# ============================================================
# 5. SAVE A COPY OF THE STAGE 3 CHECKPOINT
# ============================================================

saveRDS(
  filtered_data,
  file.path(
    analysis_dir,
    "GSE205506_stage3_checkpoint_for_targeted_analysis.rds"
  )
)

cat("\nStage 3 checkpoint copied successfully.\n")


# ============================================================
# 6. BASIC OBJECT CHECK
# ============================================================

cat("\n========== OBJECT CHECK ==========\n")

cat(
  "Cells:",
  ncol(filtered_data),
  "\n"
)

cat(
  "Genes:",
  nrow(filtered_data),
  "\n"
)

cat(
  "Metadata columns:\n"
)

print(
  colnames(filtered_data@meta.data)
)


# ============================================================
# 7. CHECK EXISTING REDUCTIONS
# ============================================================

cat("\nExisting reductions:\n")

print(
  Reductions(filtered_data)
)


# ============================================================
# 8. CHECK EXISTING CLUSTERS
# ============================================================

cat("\nCurrent identities:\n")

print(
  head(
    table(
      Idents(filtered_data)
    ),
    20
  )
)


# ============================================================
# 9. MITOCHONDRIAL QC
# ============================================================

cat(
  "\n========== MITOCHONDRIAL QC ==========\n"
)


# Calculate mitochondrial percentage if necessary

if (!"percent.mt" %in% colnames(filtered_data@meta.data)) {
  
  filtered_data[["percent.mt"]] <- PercentageFeatureSet(
    filtered_data,
    pattern = "^MT-"
  )
  
}


summary(
  filtered_data$percent.mt
)


# ============================================================
# 10. MITOCHONDRIAL DISTRIBUTION
# ============================================================

p_mt_vln <- VlnPlot(
  filtered_data,
  features = "percent.mt",
  pt.size = 0
) +
  ggtitle(
    "Mitochondrial RNA percentage"
  )

ggsave(
  file.path(
    fig_dir,
    "Figure2_percent_mt_distribution.pdf"
  ),
  p_mt_vln,
  width = 8,
  height = 6
)


# ============================================================
# 11. MITOCHONDRIAL HISTOGRAM
# ============================================================

p_mt_hist <- ggplot(
  filtered_data@meta.data,
  aes(x = percent.mt)
) +
  geom_histogram(
    bins = 100
  ) +
  theme_classic() +
  xlab(
    "Mitochondrial RNA (%)"
  ) +
  ylab(
    "Number of cells"
  ) +
  ggtitle(
    "Distribution of mitochondrial RNA"
  )

ggsave(
  file.path(
    fig_dir,
    "Figure2_percent_mt_histogram.pdf"
  ),
  p_mt_hist,
  width = 8,
  height = 6
)


# ============================================================
# 12. TEST CANDIDATE MITOCHONDRIAL THRESHOLDS
# ============================================================

mt_thresholds <- c(
  10,
  15,
  20,
  25,
  30,
  40
)

mt_summary <- data.frame(
  threshold = mt_thresholds,
  cells_remaining = sapply(
    mt_thresholds,
    function(x)
      sum(
        filtered_data$percent.mt <= x
      )
  )
)

mt_summary$percent_remaining <-
  100 *
  mt_summary$cells_remaining /
  ncol(filtered_data)

print(mt_summary)

write.csv(
  mt_summary,
  file.path(
    analysis_dir,
    "mitochondrial_threshold_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STOP POINT — MITOCHONDRIAL QC
# ============================================================
#
# DO NOT choose a threshold automatically.
#
# Inspect:
#
#   Figure2_percent_mt_distribution.pdf
#   Figure2_percent_mt_histogram.pdf
#
# and mt_summary.
#
# We will then select a biologically justified threshold.
#
# ============================================================

cat(
  "\n====================================================\n",
  "MITOCHONDRIAL QC COMPLETE.\n",
  "Review the generated figures and threshold table.\n",
  "No cells have been removed yet.\n",
  "====================================================\n"
)


# ============================================================
# 13. TEMPORARY SAFETY STOP
# ============================================================

stop(
  "STOP HERE AFTER MITOCHONDRIAL QC. Choose the QC threshold before filtering."
)

read.csv(
  "C:/Users/Admin/Documents/results/GSE205506_targeted_analysis/mitochondrial_threshold_comparison.csv"
)

PercentageFeatureSet(filtered_data, pattern = "^MT-")

# ============================================================
# GSE205506
# TARGETED ANALYSIS:
# IL1B+ MONOCYTES <-> CD8+ T CELLS
# ASSOCIATED WITH pCR / non-pCR
#
# Stage 3 is ALREADY COMPLETE.
# No new QC
# No mitochondrial filtering
# No Harmony
# No re-normalization
#
# Main analysis:
# Tumor + post-treatment samples
# pCR vs non-pCR
#
# CellChat:
# IL1B+ Monocytes -> CD8 T
# CD8 T -> IL1B+ Monocytes
# Self-signaling
#
# Then patient-level validation
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "tidyr",
  "Matrix",
  "patchwork",
  "CellChat"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "\nPackage missing: ", pkg,
        "\nPlease install it before continuing."
      )
    )
  }
  
  library(pkg, character.only = TRUE)
}


# ============================================================
# 1. DIRECTORIES
# ============================================================

base_dir <- "C:/Users/Admin/Documents/results"

stage3_dir <- file.path(
  base_dir,
  "GSE205506_STAGE3"
)

target_dir <- file.path(
  base_dir,
  "GSE205506_TARGETED_ANALYSIS"
)

if (!dir.exists(target_dir)) {
  dir.create(
    target_dir,
    recursive = TRUE
  )
}


# ============================================================
# 2. LOAD STAGE 3 OBJECT
# ============================================================

cat("\n============================================\n")
cat("LOADING STAGE 3 OBJECT\n")
cat("============================================\n\n")


# First: use an object already loaded in memory
if (exists("clustered_data") &&
    inherits(clustered_data, "Seurat")) {
  
  filtered_data <- clustered_data
  
  cat("Using object already present in memory: clustered_data\n")
  
} else {
  
  # Second: exact expected Stage 3 checkpoint
  candidate_files <- c(
    
    file.path(
      stage3_dir,
      "GSE205506_postQC_no_harmony_clustered.rds"
    ),
    
    file.path(
      base_dir,
      "GSE205506_postQC_no_harmony_clustered.rds"
    ),
    
    file.path(
      base_dir,
      "GSE205506_stage3_checkpoint_for_targeted_analysis.rds"
    )
  )
  
  existing_files <- candidate_files[
    file.exists(candidate_files)
  ]
  
  if (length(existing_files) == 0) {
    
    stop(
      paste0(
        "\nCould not find the Stage 3 Seurat object.\n",
        "Expected one of:\n",
        paste(candidate_files, collapse = "\n"),
        "\n\nDo NOT continue with another RDS automatically."
      )
    )
  }
  
  stage3_file <- existing_files[1]
  
  cat(
    "Loading:\n",
    stage3_file,
    "\n"
  )
  
  filtered_data <- readRDS(stage3_file)
}


if (!inherits(filtered_data, "Seurat")) {
  
  stop(
    "\nThe loaded object is not a Seurat object."
  )
}


cat(
  "\nCells:",
  ncol(filtered_data),
  "\nGenes:",
  nrow(filtered_data),
  "\n"
)


# ============================================================
# 3. CHECK REQUIRED METADATA
# ============================================================

required_metadata <- c(
  "patient",
  "sample_status",
  "treatment",
  "GSM",
  "seurat_clusters"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(filtered_data@meta.data)
)

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "\nMissing metadata columns:\n",
      paste(missing_metadata, collapse = ", ")
    )
  )
}


# ============================================================
# 4. CREATE TIMEPOINT
# ============================================================

cat("\n============================================\n")
cat("CREATING TIMEPOINT\n")
cat("============================================\n\n")


filtered_data$timepoint <- ifelse(
  filtered_data$treatment == "untreated",
  "pre",
  "post"
)

table(
  filtered_data$timepoint,
  filtered_data$sample_status,
  useNA = "ifany"
)


# ============================================================
# 5. PATIENT-LEVEL RESPONSE MAP
# ============================================================

cat("\n============================================\n")
cat("ADDING pCR / NON-pCR STATUS\n")
cat("============================================\n\n")


# ------------------------------------------------------------
# IMPORTANT
#
# Response is assigned at PATIENT level.
#
# P23 and P33 have only untreated/pre-treatment samples
# in this dataset and therefore are NOT part of the primary
# post-treatment response analysis.
#
# Primary post-treatment cohort:
#
# pCR:
# P11 P14 P15 P17 P19 P21 P24 P25
# P27 P28 P29 P30 P32
#
# non-pCR:
# P12 P18 P26 P31
# ------------------------------------------------------------

response_map <- data.frame(
  
  patient = c(
    "P11",
    "P12",
    "P14",
    "P15",
    "P17",
    "P18",
    "P19",
    "P21",
    "P24",
    "P25",
    "P26",
    "P27",
    "P28",
    "P29",
    "P30",
    "P31",
    "P32"
  ),
  
  response = c(
    "pCR",
    "non-pCR",
    "pCR",
    "pCR",
    "pCR",
    "non-pCR",
    "pCR",
    "pCR",
    "pCR",
    "pCR",
    "non-pCR",
    "pCR",
    "pCR",
    "pCR",
    "pCR",
    "non-pCR",
    "pCR"
  ),
  
  stringsAsFactors = FALSE
)


response_map


# Check no duplicated patients
if (anyDuplicated(response_map$patient) > 0) {
  
  stop(
    "\nDuplicated patient IDs in response_map."
  )
}


# Add response
filtered_data$response <- response_map$response[
  match(
    filtered_data$patient,
    response_map$patient
  )
]


# Make factor
filtered_data$response <- factor(
  filtered_data$response,
  levels = c(
    "pCR",
    "non-pCR"
  )
)


# Check response distribution
cat("\nResponse distribution:\n")

print(
  table(
    filtered_data$response,
    useNA = "ifany"
  )
)


# ============================================================
# 6. CHECK PATIENT / RESPONSE / TREATMENT STRUCTURE
# ============================================================

cat("\n============================================\n")
cat("PATIENT / RESPONSE / TREATMENT CHECK\n")
cat("============================================\n\n")


patient_check <- filtered_data@meta.data %>%
  
  dplyr::filter(
    sample_status == "tumor",
    timepoint == "post"
  ) %>%
  
  dplyr::group_by(
    patient,
    response,
    treatment
  ) %>%
  
  dplyr::summarise(
    cells = dplyr::n(),
    .groups = "drop"
  )

print(patient_check)


write.csv(
  patient_check,
  file.path(
    target_dir,
    "patient_response_treatment_check.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. PRIMARY ANALYSIS COHORT
# ============================================================

cat("\n============================================\n")
cat("PRIMARY COHORT\n")
cat("============================================\n\n")


analysis_data <- subset(
  
  filtered_data,
  
  subset =
    sample_status == "tumor" &
    timepoint == "post" &
    !is.na(response)
)


cat(
  "\nPrimary cohort cells:",
  ncol(analysis_data),
  "\n"
)


cat(
  "\nPatients:\n"
)

print(
  table(
    analysis_data$patient
  )
)


cat(
  "\nResponse:\n"
)

print(
  table(
    analysis_data$response
  )
)


cat(
  "\nTreatment:\n"
)

print(
  table(
    analysis_data$treatment
  )
)


# Save primary cohort
saveRDS(
  analysis_data,
  file.path(
    target_dir,
    "GSE205506_primary_post_treatment_tumor_cohort.rds"
  )
)


# ============================================================
# 8. FIND RNA ASSAY
# ============================================================

DefaultAssay(analysis_data) <- "RNA"

cat(
  "\nDefault assay:",
  DefaultAssay(analysis_data),
  "\n"
)


# ============================================================
# 9. GENE CHECK
# ============================================================

cat("\n============================================\n")
cat("CHECKING TARGET MARKERS\n")
cat("============================================\n\n")


marker_genes <- c(
  
  # T cell
  "CD3D",
  "CD3E",
  "TRBC1",
  "TRBC2",
  
  # CD8
  "CD8A",
  "CD8B",
  
  # NK
  "NKG7",
  "GNLY",
  
  # Myeloid / monocyte
  "LST1",
  "S100A8",
  "S100A9",
  "FCN1",
  "CTSS",
  "LILRB1",
  "TYMP",
  "LGALS3",
  
  # target
  "IL1B"
)

present_markers <- marker_genes[
  marker_genes %in% rownames(analysis_data)
]

missing_markers <- setdiff(
  marker_genes,
  present_markers
)

cat(
  "\nPresent markers:\n",
  paste(
    present_markers,
    collapse = ", "
  ),
  "\n"
)

if (length(missing_markers) > 0) {
  
  cat(
    "\nMissing markers:\n",
    paste(
      missing_markers,
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# 10. TARGETED MODULE SCORES
# ============================================================

cat("\n============================================\n")
cat("CALCULATING TARGETED SCORES\n")
cat("============================================\n\n")


tcell_genes <- intersect(
  c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2"
  ),
  rownames(analysis_data)
)


cd8_genes <- intersect(
  c(
    "CD8A",
    "CD8B"
  ),
  rownames(analysis_data)
)


mono_genes <- intersect(
  c(
    "LST1",
    "S100A8",
    "S100A9",
    "FCN1",
    "CTSS",
    "LILRB1",
    "TYMP",
    "LGALS3"
  ),
  rownames(analysis_data)
)


if (length(tcell_genes) < 2) {
  
  stop(
    "\nNot enough T-cell marker genes available."
  )
}


if (length(cd8_genes) < 1) {
  
  stop(
    "\nNo CD8 marker genes available."
  )
}


if (length(mono_genes) < 3) {
  
  stop(
    "\nNot enough monocyte marker genes available."
  )
}


analysis_data <- AddModuleScore(
  analysis_data,
  features = list(tcell_genes),
  name = "TcellScore",
  assay = "RNA"
)


analysis_data <- AddModuleScore(
  analysis_data,
  features = list(cd8_genes),
  name = "CD8Score",
  assay = "RNA"
)


analysis_data <- AddModuleScore(
  analysis_data,
  features = list(mono_genes),
  name = "MonoScore",
  assay = "RNA"
)


# ============================================================
# 11. EXAMINE SCORE DISTRIBUTIONS
# ============================================================

cat("\n============================================\n")
cat("SCORE DISTRIBUTIONS\n")
cat("============================================\n\n")


print(
  summary(
    analysis_data$TcellScore1
  )
)

print(
  summary(
    analysis_data$CD8Score1
  )
)

print(
  summary(
    analysis_data$MonoScore1
  )
)


# ============================================================
# 12. VISUAL QC OF TARGET POPULATIONS
# ============================================================

cat("\n============================================\n")
cat("GENERATING TARGET MARKER FIGURES\n")
cat("============================================\n\n")


p_target_features <- FeaturePlot(
  analysis_data,
  features = intersect(
    c(
      "CD3D",
      "CD3E",
      "CD8A",
      "CD8B",
      "LST1",
      "FCN1",
      "S100A8",
      "S100A9",
      "IL1B"
    ),
    rownames(analysis_data)
  ),
  reduction = "umap",
  ncol = 3,
  order = TRUE
)

ggsave(
  file.path(
    target_dir,
    "01_target_marker_FeaturePlot.pdf"
  ),
  p_target_features,
  width = 14,
  height = 10
)


p_target_scores <- FeaturePlot(
  analysis_data,
  features = c(
    "TcellScore1",
    "CD8Score1",
    "MonoScore1"
  ),
  reduction = "umap",
  ncol = 3,
  order = TRUE
)

ggsave(
  file.path(
    target_dir,
    "02_target_module_scores.pdf"
  ),
  p_target_scores,
  width = 14,
  height = 5
)


# ============================================================
# 13. DEFINE TARGET CELLS
# ============================================================

cat("\n============================================\n")
cat("DEFINING TARGET CELL POPULATIONS\n")
cat("============================================\n\n")


# ------------------------------------------------------------
# Expression matrix
#
# We use normalized RNA expression.
# A detected gene has expression > 0.
# ------------------------------------------------------------

expr <- GetAssayData(
  analysis_data,
  assay = "RNA",
  slot = "data"
)


# ------------------------------------------------------------
# T-cell detection
# At least one T-cell receptor/CD3 marker detected.
# ------------------------------------------------------------

tcell_detected <- rep(
  FALSE,
  ncol(analysis_data)
)

tcell_detected <- Matrix::colSums(
  expr[
    tcell_genes,
    ,
    drop = FALSE
  ] > 0
) >= 2


# ------------------------------------------------------------
# CD8 detection
# At least one CD8 gene detected.
# ------------------------------------------------------------

cd8_detected <- Matrix::colSums(
  expr[
    cd8_genes,
    ,
    drop = FALSE
  ] > 0
) >= 1


# ------------------------------------------------------------
# Monocyte detection
# At least 3 monocyte markers detected.
# ------------------------------------------------------------

mono_detected <- Matrix::colSums(
  expr[
    mono_genes,
    ,
    drop = FALSE
  ] > 0
) >= 3


# ------------------------------------------------------------
# IL1B detection
# ------------------------------------------------------------

il1b_detected <- rep(
  FALSE,
  ncol(analysis_data)
)

if ("IL1B" %in% rownames(expr)) {
  
  il1b_detected <- as.numeric(
    expr["IL1B", ]
  ) > 0
}


# ============================================================
# 14. ADD CELL LABELS
# ============================================================

analysis_data$target_cell <- "Other"


# CD8 T
analysis_data$target_cell[
  tcell_detected &
    cd8_detected
] <- "CD8_T"


# IL1B+ monocyte
analysis_data$target_cell[
  mono_detected &
    il1b_detected
] <- "IL1B_Mono"


analysis_data$target_cell <- factor(
  analysis_data$target_cell,
  levels = c(
    "IL1B_Mono",
    "CD8_T",
    "Other"
  )
)


# ============================================================
# 15. TARGET CELL COUNTS
# ============================================================

cat("\n============================================\n")
cat("TARGET CELL COUNTS\n")
cat("============================================\n\n")


print(
  table(
    analysis_data$target_cell
  )
)


print(
  table(
    analysis_data$target_cell,
    analysis_data$response
  )
)


# ============================================================
# 16. PATIENT-LEVEL TARGET CELL COUNTS\n
# ============================================================

target_patient_counts <- analysis_data@meta.data %>%
  
  dplyr::filter(
    target_cell != "Other"
  ) %>%
  
  dplyr::count(
    patient,
    response,
    target_cell,
    name = "cells"
  )

print(
  target_patient_counts
)


write.csv(
  target_patient_counts,
  file.path(
    target_dir,
    "03_target_cell_counts_by_patient.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CHECK MINIMUM CELL NUMBERS
# ============================================================

cat("\n============================================\n")
cat("MINIMUM CELL CHECK\n")
cat("============================================\n\n")


cell_count_table <- analysis_data@meta.data %>%
  
  dplyr::filter(
    target_cell != "Other"
  ) %>%
  
  dplyr::count(
    patient,
    response,
    target_cell,
    name = "n"
  ) %>%
  
  tidyr::pivot_wider(
    names_from = target_cell,
    values_from = n,
    values_fill = 0
  )


print(
  cell_count_table
)


# ------------------------------------------------------------
# For CellChat:
# minimum 10 cells per population.
#
# Patients below this threshold are excluded from the
# patient-level CellChat validation, NOT from the global
# analysis.
# ------------------------------------------------------------

eligible_patients <- cell_count_table %>%
  
  dplyr::filter(
    IL1B_Mono >= 10,
    CD8_T >= 10
  ) %>%
  
  dplyr::pull(
    patient
  )


cat(
  "\nPatients eligible for patient-level CellChat:\n"
)

print(
  eligible_patients
)


write.csv(
  cell_count_table,
  file.path(
    target_dir,
    "04_patient_target_cell_eligibility.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. TARGET-ONLY OBJECT
# ============================================================

target_data <- subset(
  
  analysis_data,
  
  subset =
    target_cell %in%
    c(
      "IL1B_Mono",
      "CD8_T"
    )
)


cat(
  "\nTarget-only cells:",
  ncol(target_data),
  "\n"
)


print(
  table(
    target_data$target_cell
  )
)

print(
  table(
    target_data$target_cell,
    target_data$response
  )
)


# Save
saveRDS(
  target_data,
  file.path(
    target_dir,
    "GSE205506_IL1B_Mono_CD8T_target_cells.rds"
  )
)


# ============================================================
# 19. TARGET POPULATION VISUALIZATION
# ============================================================

p_target_umap <- DimPlot(
  analysis_data,
  reduction = "umap",
  group.by = "target_cell",
  cols = NULL
) +
  ggtitle(
    "Target populations: IL1B+ monocytes and CD8+ T cells"
  )

ggsave(
  file.path(
    target_dir,
    "05_target_populations_UMAP.pdf"
  ),
  p_target_umap,
  width = 8,
  height = 6
)


# ============================================================
# 20. RESPONSE VISUALIZATION
# ============================================================

p_response_umap <- DimPlot(
  analysis_data,
  reduction = "umap",
  group.by = "response"
) +
  ggtitle(
    "Post-treatment tumor cells: pCR vs non-pCR"
  )

ggsave(
  file.path(
    target_dir,
    "06_response_UMAP.pdf"
  ),
  p_response_umap,
  width = 8,
  height = 6
)


# ============================================================
# 21. TARGET MARKER DOTPLOT
# ============================================================

dot_genes <- intersect(
  c(
    "IL1B",
    "LST1",
    "FCN1",
    "S100A8",
    "S100A9",
    "CTSS",
    "CD3D",
    "CD3E",
    "TRBC1",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY"
  ),
  rownames(analysis_data)
)


p_dot <- DotPlot(
  analysis_data,
  features = dot_genes,
  group.by = "target_cell"
) +
  RotatedAxis() +
  ggtitle(
    "Marker validation of target populations"
  )

ggsave(
  file.path(
    target_dir,
    "07_target_marker_DotPlot.pdf"
  ),
  p_dot,
  width = 14,
  height = 6
)


# ============================================================
# 22. SAVE CHECKPOINT BEFORE CELLCHAT
# ============================================================

saveRDS(
  analysis_data,
  file.path(
    target_dir,
    "GSE205506_targeted_annotation_checkpoint.rds"
  )
)

saveRDS(
  target_data,
  file.path(
    target_dir,
    "GSE205506_target_only_checkpoint.rds"
  )
)


# ============================================================
# 23. CHECK CELLCHAT VERSION
# ============================================================

cat("\n============================================\n")
cat("CELLCHAT VERSION\n")
cat("============================================\n\n")


cat(
  "CellChat version: ",
  as.character(
    packageVersion("CellChat")
  ),
  "\n"
)


# ============================================================
# 24. CELLCHAT FUNCTION
# ============================================================

run_cellchat_target <- function(
    seu,
    response_label,
    output_prefix
) {
  
  cat(
    "\n--------------------------------------------\n"
  )
  
  cat(
    "Running CellChat:",
    response_label,
    "\n"
  )
  
  cat(
    "--------------------------------------------\n"
  )
  
  
  # ----------------------------------------------------------
  # Subset response
  # ----------------------------------------------------------
  
  obj <- subset(
    seu,
    subset =
      response == response_label
  )
  
  
  # Remove unused factor levels
  obj$target_cell <- droplevels(
    factor(
      obj$target_cell
    )
  )
  
  
  # ----------------------------------------------------------
  # Cell counts
  # ----------------------------------------------------------
  
  print(
    table(
      obj$target_cell
    )
  )
  
  
  # ----------------------------------------------------------
  # Minimum cells
  # ----------------------------------------------------------
  
  cell_counts <- table(
    obj$target_cell
  )
  
  
  if (any(cell_counts < 10)) {
    
    stop(
      paste0(
        "\nCellChat cannot be run safely for ",
        response_label,
        ".\n",
        "At least one target population has <10 cells."
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Seurat v5 layer handling
  #
  # Only done AFTER target selection.
  # This avoids the large memory problem on the full object.
  # ----------------------------------------------------------
  
  if (
    "JoinLayers" %in%
    getNamespaceExports("SeuratObject")
  ) {
    
    obj <- JoinLayers(
      obj,
      assay = "RNA"
    )
  }
  
  
  # ----------------------------------------------------------
  # Create CellChat object
  # ----------------------------------------------------------
  
  cellchat <- createCellChat(
    object = obj,
    group.by = "target_cell"
  )
  
  
  # ----------------------------------------------------------
  # Human database
  # ----------------------------------------------------------
  
  cellchat@DB <- CellChatDB.human
  
  
  # ----------------------------------------------------------
  # Subset database
  # ----------------------------------------------------------
  
  cellchat <- subsetData(
    cellchat
  )
  
  
  # ----------------------------------------------------------
  # Overexpressed genes
  # ----------------------------------------------------------
  
  cellchat <- identifyOverExpressedGenes(
    cellchat
  )
  
  
  # ----------------------------------------------------------
  # Overexpressed ligand-receptor interactions
  # ----------------------------------------------------------
  
  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )
  
  
  # ----------------------------------------------------------
  # Protein interaction projection
  # ----------------------------------------------------------
  
  cellchat <- projectData(
    cellchat,
    PPI.human
  )
  
  
  # ----------------------------------------------------------
  # Communication probability
  # ----------------------------------------------------------
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )
  
  
  # ----------------------------------------------------------
  # Filter low-cell populations
  # ----------------------------------------------------------
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )
  
  
  # ----------------------------------------------------------
  # Pathway-level probabilities
  # ----------------------------------------------------------
  
  cellchat <- computeCommunProbPathway(
    cellchat
  )
  
  
  # ----------------------------------------------------------
  # Aggregate network
  # ----------------------------------------------------------
  
  cellchat <- aggregateNet(
    cellchat
  )
  
  
  # ----------------------------------------------------------
  # Save CellChat object
  # ----------------------------------------------------------
  
  saveRDS(
    cellchat,
    file.path(
      target_dir,
      paste0(
        output_prefix,
        "_CellChat.rds"
      )
    )
  )
  
  
  # ----------------------------------------------------------
  # Extract LR communication
  # ----------------------------------------------------------
  
  communication <- subsetCommunication(
    cellchat
  )
  
  
  if (!is.null(communication) &&
      nrow(communication) > 0) {
    
    write.csv(
      communication,
      file.path(
        target_dir,
        paste0(
          output_prefix,
          "_all_communications.csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Target directions
  # ----------------------------------------------------------
  
  if (!is.null(communication) &&
      nrow(communication) > 0) {
    
    mono_to_cd8 <- communication %>%
      
      dplyr::filter(
        source == "IL1B_Mono",
        target == "CD8_T"
      )
    
    cd8_to_mono <- communication %>%
      
      dplyr::filter(
        source == "CD8_T",
        target == "IL1B_Mono"
      )
    
    mono_self <- communication %>%
      
      dplyr::filter(
        source == "IL1B_Mono",
        target == "IL1B_Mono"
      )
    
    cd8_self <- communication %>%
      
      dplyr::filter(
        source == "CD8_T",
        target == "CD8_T"
      )
    
    
    write.csv(
      mono_to_cd8,
      file.path(
        target_dir,
        paste0(
          output_prefix,
          "_IL1B_Mono_to_CD8T.csv"
        )
      ),
      row.names = FALSE
    )
    
    
    write.csv(
      cd8_to_mono,
      file.path(
        target_dir,
        paste0(
          output_prefix,
          "_CD8T_to_IL1B_Mono.csv"
        )
      ),
      row.names = FALSE
    )
    
    
    write.csv(
      mono_self,
      file.path(
        target_dir,
        paste0(
          output_prefix,
          "_IL1B_Mono_self.csv"
        )
      ),
      row.names = FALSE
    )
    
    
    write.csv(
      cd8_self,
      file.path(
        target_dir,
        paste0(
          output_prefix,
          "_CD8T_self.csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  
  return(cellchat)
}


# ============================================================
# 25. RUN CELLCHAT — pCR
# ============================================================

cat("\n============================================\n")
cat("CELLCHAT pCR\n")
cat("============================================\n\n")


cellchat_pCR <- run_cellchat_target(
  target_data,
  response_label = "pCR",
  output_prefix = "pCR"
)


# ============================================================
# 26. RUN CELLCHAT — NON-pCR
# ============================================================

cat("\n============================================\n")
cat("CELLCHAT NON-pCR\n")
cat("============================================\n\n")


cellchat_nonpCR <- run_cellchat_target(
  target_data,
  response_label = "non-pCR",
  output_prefix = "nonpCR"
)


# ============================================================
# 27. COMPARE pCR vs NON-pCR
# ============================================================

cat("\n============================================\n")
cat("COMPARING pCR VS NON-pCR\n")
cat("============================================\n\n")


comm_pCR <- subsetCommunication(
  cellchat_pCR
)


comm_nonpCR <- subsetCommunication(
  cellchat_nonpCR
)


# ------------------------------------------------------------
# Standardize columns
# ------------------------------------------------------------

required_comm_cols <- c(
  "source",
  "target",
  "ligand",
  "receptor",
  "prob",
  "pval"
)


missing_pcr_cols <- setdiff(
  required_comm_cols,
  colnames(comm_pCR)
)

missing_nonpcr_cols <- setdiff(
  required_comm_cols,
  colnames(comm_nonpCR)
)


if (
  length(missing_pcr_cols) > 0 |
  length(missing_nonpcr_cols) > 0
) {
  
  cat(
    "\nWARNING: CellChat output column names differ.\n"
  )
  
  cat(
    "\npCR missing:\n",
    paste(
      missing_pcr_cols,
      collapse = ", "
    ),
    "\n"
  )
  
  cat(
    "\nnon-pCR missing:\n",
    paste(
      missing_nonpcr_cols,
      collapse = ", "
    ),
    "\n"
  )
  
} else {
  
  
  # ----------------------------------------------------------
  # Keep only the three directions of interest
  # ----------------------------------------------------------
  
  comm_pCR_target <- comm_pCR %>%
    
    dplyr::filter(
      
      (
        source == "IL1B_Mono" &
          target == "CD8_T"
      ) |
        
        (
          source == "CD8_T" &
            target == "IL1B_Mono"
        ) |
        
        (
          source == "IL1B_Mono" &
            target == "IL1B_Mono"
        ) |
        
        (
          source == "CD8_T" &
            target == "CD8_T"
        )
    )
  
  
  comm_nonpCR_target <- comm_nonpCR %>%
    
    dplyr::filter(
      
      (
        source == "IL1B_Mono" &
          target == "CD8_T"
      ) |
        
        (
          source == "CD8_T" &
            target == "IL1B_Mono"
        ) |
        
        (
          source == "IL1B_Mono" &
            target == "IL1B_Mono"
        ) |
        
        (
          source == "CD8_T" &
            target == "CD8_T"
        )
    )
  
  
  # ----------------------------------------------------------
  # Add group
  # ----------------------------------------------------------
  
  comm_pCR_target$response <- "pCR"
  
  comm_nonpCR_target$response <- "non-pCR"
  
  
  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------
  
  all_target_communications <- bind_rows(
    comm_pCR_target,
    comm_nonpCR_target
  )
  
  
  write.csv(
    all_target_communications,
    file.path(
      target_dir,
      "08_target_communications_pCR_vs_nonpCR.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 28. CREATE LR COMPARISON TABLE
  # ==========================================================
  
  lr_pCR <- comm_pCR_target %>%
    
    dplyr::select(
      source,
      target,
      ligand,
      receptor,
      prob,
      pval
    ) %>%
    
    dplyr::rename(
      prob_pCR = prob,
      pval_pCR = pval
    )
  
  
  lr_nonpCR <- comm_nonpCR_target %>%
    
    dplyr::select(
      source,
      target,
      ligand,
      receptor,
      prob,
      pval
    ) %>%
    
    dplyr::rename(
      prob_nonpCR = prob,
      pval_nonpCR = pval
    )
  
  
  lr_comparison <- full_join(
    
    lr_pCR,
    
    lr_nonpCR,
    
    by = c(
      "source",
      "target",
      "ligand",
      "receptor"
    )
  )
  
  
  lr_comparison <- lr_comparison %>%
    
    mutate(
      
      prob_pCR = ifelse(
        is.na(prob_pCR),
        0,
        prob_pCR
      ),
      
      prob_nonpCR = ifelse(
        is.na(prob_nonpCR),
        0,
        prob_nonpCR
      ),
      
      delta_prob =
        prob_pCR -
        prob_nonpCR,
      
      fold_change_pCR_over_nonpCR =
        (prob_pCR + 1e-6) /
        (prob_nonpCR + 1e-6),
      
      log2FC_pCR_over_nonpCR =
        log2(
          fold_change_pCR_over_nonpCR
        )
    )
  
  
  # ----------------------------------------------------------
  # Sort by absolute difference
  #
  # This is NOT inferential patient-level significance.
  # It is a CellChat discovery ranking.
  # ----------------------------------------------------------
  
  lr_comparison <- lr_comparison %>%
    
    arrange(
      desc(
        abs(
          delta_prob
        )
      )
    )
  
  
  write.csv(
    lr_comparison,
    file.path(
      target_dir,
      "09_LR_comparison_pCR_vs_nonpCR.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 29. TOP DIFFERENTIAL INTERACTIONS
  # ==========================================================
  
  top_lr <- lr_comparison %>%
    
    slice_head(
      n = 30
    )
  
  
  write.csv(
    top_lr,
    file.path(
      target_dir,
      "10_top_30_LR_discovery_candidates.csv"
    ),
    row.names = FALSE
  )
  
  
  print(
    top_lr
  )
}


# ============================================================
# 30. CELLCHAT BUBBLE PLOTS
# ============================================================

cat("\n============================================\n")
cat("CELLCHAT BUBBLE PLOTS\n")
cat("============================================\n\n")


pdf(
  file.path(
    target_dir,
    "11_CellChat_bubble_pCR.pdf"
  ),
  width = 12,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_pCR,
    sources.use = "IL1B_Mono",
    targets.use = "CD8_T",
    remove.isolate = FALSE
  )
)

dev.off()


pdf(
  file.path(
    target_dir,
    "12_CellChat_bubble_nonpCR.pdf"
  ),
  width = 12,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_nonpCR,
    sources.use = "IL1B_Mono",
    targets.use = "CD8_T",
    remove.isolate = FALSE
  )
)

dev.off()


# Reverse direction

pdf(
  file.path(
    target_dir,
    "13_CellChat_bubble_CD8T_to_IL1BMono_pCR.pdf"
  ),
  width = 12,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_pCR,
    sources.use = "CD8_T",
    targets.use = "IL1B_Mono",
    remove.isolate = FALSE
  )
)

dev.off()


pdf(
  file.path(
    target_dir,
    "14_CellChat_bubble_CD8T_to_IL1BMono_nonpCR.pdf"
  ),
  width = 12,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_nonpCR,
    sources.use = "CD8_T",
    targets.use = "IL1B_Mono",
    remove.isolate = FALSE
  )
)

dev.off()


# ============================================================
# 31. PATIENT-LEVEL VALIDATION
# ============================================================

cat("\n============================================\n")
cat("PATIENT-LEVEL VALIDATION\n")
cat("============================================\n\n")


# ------------------------------------------------------------
# IMPORTANT:
#
# CellChat pooled by response group is DISCOVERY.
#
# We now run CellChat independently per patient.
#
# This prevents a single patient with many cells from
# determining the apparent response-associated interaction.
# ------------------------------------------------------------


patient_results <- list()

patient_summary <- data.frame()


for (pt in eligible_patients) {
  
  cat(
    "\nPatient:",
    pt,
    "\n"
  )
  
  
  pt_obj <- subset(
    target_data,
    subset =
      patient == pt
  )
  
  
  counts_pt <- table(
    pt_obj$target_cell
  )
  
  
  if (
    !all(
      c(
        "IL1B_Mono",
        "CD8_T"
      ) %in%
      names(counts_pt)
    )
  ) {
    
    cat(
      "Skipping:",
      pt,
      "- missing population\n"
    )
    
    next
  }
  
  
  if (
    any(
      counts_pt[
        c(
          "IL1B_Mono",
          "CD8_T"
        )
      ] < 10
    )
  ) {
    
    cat(
      "Skipping:",
      pt,
      "- fewer than 10 cells in a target population\n"
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # Join layers only on tiny patient object
  # ----------------------------------------------------------
  
  if (
    "JoinLayers" %in%
    getNamespaceExports("SeuratObject")
  ) {
    
    pt_obj <- JoinLayers(
      pt_obj,
      assay = "RNA"
    )
  }
  
  
  # ----------------------------------------------------------
  # CellChat
  # ----------------------------------------------------------
  
  cc_pt <- createCellChat(
    object = pt_obj,
    group.by = "target_cell"
  )
  
  
  cc_pt@DB <- CellChatDB.human
  
  
  cc_pt <- subsetData(
    cc_pt
  )
  
  
  cc_pt <- identifyOverExpressedGenes(
    cc_pt
  )
  
  
  cc_pt <- identifyOverExpressedInteractions(
    cc_pt
  )
  
  
  cc_pt <- projectData(
    cc_pt,
    PPI.human
  )
  
  
  cc_pt <- computeCommunProb(
    cc_pt,
    type = "triMean"
  )
  
  
  cc_pt <- filterCommunication(
    cc_pt,
    min.cells = 10
  )
  
  
  # ----------------------------------------------------------
  # Extract communication
  # ----------------------------------------------------------
  
  comm_pt <- subsetCommunication(
    cc_pt
  )
  
  
  if (
    !is.null(comm_pt) &&
    nrow(comm_pt) > 0
  ) {
    
    # IL1B Mono -> CD8
    m2t <- comm_pt %>%
      
      filter(
        source == "IL1B_Mono",
        target == "CD8_T"
      ) %>%
      
      mutate(
        patient = pt,
        response = unique(
          pt_obj$response
        ),
        direction =
          "IL1B_Mono_to_CD8T"
      )
    
    
    # CD8 -> IL1B Mono
    t2m <- comm_pt %>%
      
      filter(
        source == "CD8_T",
        target == "IL1B_Mono"
      ) %>%
      
      mutate(
        patient = pt,
        response = unique(
          pt_obj$response
        ),
        direction =
          "CD8T_to_IL1B_Mono"
      )
    
    
    # Self mono
    mono_self_pt <- comm_pt %>%
      
      filter(
        source == "IL1B_Mono",
        target == "IL1B_Mono"
      ) %>%
      
      mutate(
        patient = pt,
        response = unique(
          pt_obj$response
        ),
        direction =
          "IL1B_Mono_self"
      )
    
    
    # Self CD8
    cd8_self_pt <- comm_pt %>%
      
      filter(
        source == "CD8_T",
        target == "CD8_T"
      ) %>%
      
      mutate(
        patient = pt,
        response = unique(
          pt_obj$response
        ),
        direction =
          "CD8T_self"
      )
    
    
    patient_results[[pt]] <- bind_rows(
      m2t,
      t2m,
      mono_self_pt,
      cd8_self_pt
    )
    
    
    patient_summary <- bind_rows(
      
      patient_summary,
      
      data.frame(
        patient = pt,
        response =
          as.character(
            unique(
              pt_obj$response
            )
          ),
        IL1B_Mono =
          counts_pt[
            "IL1B_Mono"
          ],
        CD8_T =
          counts_pt[
            "CD8_T"
          ],
        n_interactions =
          nrow(comm_pt)
      )
    )
  }
}


# ============================================================
# 32. COMBINE PATIENT-LEVEL RESULTS
# ============================================================

patient_lr <- bind_rows(
  patient_results
)


write.csv(
  patient_lr,
  file.path(
    target_dir,
    "15_patient_level_CellChat_interactions.csv"
  ),
  row.names = FALSE
)


write.csv(
  patient_summary,
  file.path(
    target_dir,
    "16_patient_level_CellChat_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 33. PATIENT-LEVEL SUMMARY OF TOP LR PAIRS
# ============================================================

if (
  nrow(patient_lr) > 0
) {
  
  patient_lr_summary <- patient_lr %>%
    
    group_by(
      direction,
      ligand,
      receptor,
      response
    ) %>%
    
    summarise(
      
      n_patients =
        n_distinct(
          patient
        ),
      
      mean_prob =
        mean(
          prob,
          na.rm = TRUE
        ),
      
      median_prob =
        median(
          prob,
          na.rm = TRUE
        ),
      
      mean_pval =
        mean(
          pval,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  write.csv(
    patient_lr_summary,
    file.path(
      target_dir,
      "17_patient_level_LR_summary.csv"
    ),
    row.names = FALSE
  )
  
  
  print(
    patient_lr_summary
  )
}


# ============================================================
# 34. PATIENT-LEVEL pCR vs NON-pCR COMPARISON
# ============================================================

if (
  nrow(patient_lr) > 0
) {
  
  patient_lr_test <- patient_lr %>%
    
    group_by(
      direction,
      ligand,
      receptor
    ) %>%
    
    summarise(
      
      n_pCR =
        sum(
          response == "pCR"
        ),
      
      n_nonpCR =
        sum(
          response == "non-pCR"
        ),
      
      median_pCR =
        median(
          prob[
            response == "pCR"
          ],
          na.rm = TRUE
        ),
      
      median_nonpCR =
        median(
          prob[
            response == "non-pCR"
          ],
          na.rm = TRUE
        ),
      
      p_value_wilcox = {
        
        x <- prob[
          response == "pCR"
        ]
        
        y <- prob[
          response == "non-pCR"
        ]
        
        if (
          length(
            unique(
              na.omit(x)
            )
          ) > 1 ||
          length(
            unique(
              na.omit(y)
            )
          ) > 1
        ) {
          
          tryCatch(
            
            wilcox.test(
              x,
              y,
              exact = FALSE
            )$p.value,
            
            error = function(e) {
              NA_real_
            }
          )
          
        } else {
          NA_real_
        }
      },
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      FDR =
        p.adjust(
          p_value_wilcox,
          method = "BH"
        ),
      
      delta_median =
        median_pCR -
        median_nonpCR
    ) %>%
    
    arrange(
      FDR
    )
  
  
  write.csv(
    patient_lr_test,
    file.path(
      target_dir,
      "18_patient_level_LR_pCR_vs_nonpCR_statistics.csv"
    ),
    row.names = FALSE
  )
  
  
  print(
    patient_lr_test
  )
}


# ============================================================
# 35. FINAL CHECKPOINT
# ============================================================

saveRDS(
  analysis_data,
  file.path(
    target_dir,
    "GSE205506_FINAL_targeted_analysis_object.rds"
  )
)


saveRDS(
  target_data,
  file.path(
    target_dir,
    "GSE205506_FINAL_target_cells.rds"
  )
)


# ============================================================
# 36. FINAL MESSAGE
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  target_dir,
  "\n\n"
)

cat(
  "Main objects:\n"
)

cat(
  "- GSE205506_FINAL_targeted_analysis_object.rds\n"
)

cat(
  "- GSE205506_FINAL_target_cells.rds\n"
)

cat(
  "\nMain CellChat files:\n"
)

cat(
  "- pCR_CellChat.rds\n"
)

cat(
  "- nonpCR_CellChat.rds\n"
)

cat(
  "\nMain LR comparison:\n"
)

cat(
  "- 09_LR_comparison_pCR_vs_nonpCR.csv\n"
)

cat(
  "\nPatient-level validation:\n"
)

cat(
  "- 18_patient_level_LR_pCR_vs_nonpCR_statistics.csv\n"
)

cat("\n============================================================\n")

FeaturePlot(..., reduction = "umap")

FeaturePlot(..., reduction = "umap")

Reductions(analysis_data)


p_target_features <- FeaturePlot(
  analysis_data,
  features = intersect(
    c(
      "CD3D",
      "CD3E",
      "CD8A",
      "CD8B",
      "LST1",
      "FCN1",
      "S100A8",
      "S100A9",
      "IL1B"
    ),
    rownames(analysis_data)
  ),
  reduction = "umap_before",
  ncol = 3,
  order = TRUE
)

ggsave(
  file.path(
    target_dir,
    "01_target_marker_FeaturePlot.pdf"
  ),
  p_target_features,
  width = 14,
  height = 10
)

# ============================================================
# 17. TARGET UMAP
# ============================================================

# Detect the available UMAP reduction

cat("target_cell existe :", "target_cell" %in% colnames(analysis_data@meta.data), "\n\n"); cat("Colonnes pertinentes :\n"); print(colnames(analysis_data@meta.data)[grepl("target|CD8|Mono|IL1B|score", colnames(analysis_data@meta.data), ignore.case = TRUE)])
analysis_data$target_cell <- ifelse(
  analysis_data$TcellScore1 > 0 & analysis_data$CD8Score1 > 0,
  "CD8_T",
  ifelse(
    analysis_data$MonoScore1 > 0 &
      GetAssayData(analysis_data, assay = "RNA", layer = "data")["IL1B", ] > 0,
    "IL1B_Mono",
    "Other"
  )
); analysis_data$target_cell <- factor(analysis_data$target_cell, levels = c("IL1B_Mono", "CD8_T", "Other")); print(table(analysis_data$target_cell))

analysis_data <- JoinLayers(analysis_data, assay = "RNA"); expr_IL1B <- GetAssayData(analysis_data, assay = "RNA", layer = "data")["IL1B", ]; analysis_data$target_cell <- ifelse(analysis_data$TcellScore1 > 0 & analysis_data$CD8Score1 > 0, "CD8_T", ifelse(analysis_data$MonoScore1 > 0 & expr_IL1B > 0, "IL1B_Mono", "Other")); analysis_data$target_cell <- factor(analysis_data$target_cell, levels = c("IL1B_Mono", "CD8_T", "Other")); print(table(analysis_data$target_cell))

target_reduction <- if ("umap_before" %in% Reductions(analysis_data)) "umap_before" else if ("umap" %in% Reductions(analysis_data)) "umap" else stop("No UMAP reduction found."); p_target_umap <- DimPlot(analysis_data, reduction = target_reduction, group.by = "target_cell") + ggtitle("IL1B+ monocytes and CD8+ T cells"); ggsave(file.path(target_dir, "03_target_populations_UMAP.pdf"), p_target_umap, width = 8, height = 6); print(p_target_umap)

p_target_umap <- DimPlot(
  analysis_data,
  reduction = target_reduction,
  group.by = "target_cell",
  cols = NULL
) +
  ggtitle(
    "Target populations: IL1B+ monocytes and CD8+ T cells"
  )

ggsave(
  file.path(
    target_dir,
    "03_target_populations_UMAP.pdf"
  ),
  p_target_umap,
  width = 8,
  height = 6
)

print(p_target_umap)
# ============================================================
# 20. RESPONSE UMAP
# ============================================================

p_response_umap <- DimPlot(
  analysis_data,
  reduction = target_reduction,
  group.by = "response"
) +
  ggtitle(
    "Post-treatment tumor cells: pCR vs non-pCR"
  )

ggsave(
  file.path(
    target_dir,
    "04_response_UMAP.pdf"
  ),
  p_response_umap,
  width = 8,
  height = 6
)

print(p_response_umap)


# ============================================================
# 21. TARGET MARKER VISUALIZATION
# ============================================================

target_markers <- c(
  "IL1B",
  "LST1",
  "FCN1",
  "S100A8",
  "S100A9",
  "CD3D",
  "CD3E",
  "TRBC1",
  "TRBC2",
  "CD8A",
  "CD8B"
)

target_markers <- target_markers[
  target_markers %in% rownames(analysis_data)
]

p_target_dotplot <- DotPlot(
  analysis_data,
  features = target_markers,
  group.by = "target_cell"
) +
  RotatedAxis() +
  ggtitle(
    "Target population marker expression"
  )

ggsave(
  file.path(
    target_dir,
    "05_target_population_DotPlot.pdf"
  ),
  p_target_dotplot,
  width = 11,
  height = 6
)

print(p_target_dotplot)


# ============================================================
# 22. PATIENT-LEVEL TARGET COUNTS
# ============================================================

target_patient_counts <- analysis_data@meta.data %>%
  dplyr::filter(
    target_cell != "Other"
  ) %>%
  dplyr::count(
    patient,
    target_cell,
    name = "n_cells"
  )

print(target_patient_counts)

write.csv(
  target_patient_counts,
  file.path(
    target_dir,
    "06_target_population_counts_by_patient.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. ELIGIBLE PATIENTS
# ============================================================

eligible_patients <- target_patient_counts %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(
    n_target_cells = sum(n_cells),
    n_target_types = dplyr::n_distinct(target_cell),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_target_cells >= 10
  )

print(eligible_patients)

write.csv(
  eligible_patients,
  file.path(
    target_dir,
    "07_eligible_patients.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. TARGET-ONLY OBJECT
# ============================================================

target_data <- subset(
  analysis_data,
  subset = target_cell != "Other"
)

target_data$target_cell <- droplevels(
  target_data$target_cell
)

print(
  table(
    target_data$target_cell,
    useNA = "ifany"
  )
)

print(
  table(
    target_data$response,
    target_data$target_cell,
    useNA = "ifany"
  )
)


# ============================================================
# 25. CELLCHAT FUNCTION
# ============================================================

run_cellchat <- function(
    seu,
    response_name,
    prefix
) {
  
  message(
    "\n============================================================"
  )
  
  message(
    "Running CellChat: ",
    response_name
  )
  
  message(
    "============================================================"
  )
  
  obj <- subset(
    seu,
    subset = response == response_name
  )
  
  obj$target_cell <- droplevels(
    factor(obj$target_cell)
  )
  
  print(
    table(
      obj$target_cell,
      useNA = "ifany"
    )
  )
  
  if (ncol(obj) == 0) {
    stop(
      "No cells found for response group: ",
      response_name
    )
  }
  
  cell_counts <- table(
    obj$target_cell
  )
  
  if (any(cell_counts < 10)) {
    stop(
      "At least one target population has fewer than 10 cells in ",
      response_name,
      "."
    )
  }
  
  # ----------------------------------------------------------
  # Seurat v5: join RNA layers before CellChat
  # ----------------------------------------------------------
  
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    
    obj <- JoinLayers(
      obj,
      assay = "RNA"
    )
    
  }
  
  # ----------------------------------------------------------
  # Create CellChat object
  # ----------------------------------------------------------
  
  cellchat <- createCellChat(
    object = obj,
    group.by = "target_cell"
  )
  
  # ----------------------------------------------------------
  # Human database
  # ----------------------------------------------------------
  
  cellchat@DB <- CellChatDB.human
  
  # ----------------------------------------------------------
  # Preprocessing
  # ----------------------------------------------------------
  
  cellchat <- subsetData(
    cellchat
  )
  
  cellchat <- identifyOverExpressedGenes(
    cellchat
  )
  
  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )
  
  cellchat <- projectData(
    cellchat,
    PPI.human
  )
  
  # ----------------------------------------------------------
  # Communication probability
  # ----------------------------------------------------------
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )
  
  cellchat <- computeCommunProbPathway(
    cellchat
  )
  
  cellchat <- aggregateNet(
    cellchat
  )
  
  # ----------------------------------------------------------
  # Save CellChat object
  # ----------------------------------------------------------
  
  saveRDS(
    cellchat,
    file.path(
      target_dir,
      paste0(
        prefix,
        "_CellChat.rds"
      )
    )
  )
  
  message(
    "CellChat object saved."
  )
  
  return(
    cellchat
  )
}


# ============================================================
# 26. CELLCHAT - pCR
# ============================================================

cellchat_pCR <- run_cellchat(
  target_data,
  response_name = "pCR",
  prefix = "08_pCR"
)


# ============================================================
# 27. CELLCHAT - NON-pCR
# ============================================================

cellchat_nonpCR <- run_cellchat(
  target_data,
  response_name = "non-pCR",
  prefix = "09_nonpCR"
)


# ============================================================
# 28. EXTRACT COMMUNICATION
# ============================================================

comm_pCR <- subsetCommunication(
  cellchat_pCR
)

comm_nonpCR <- subsetCommunication(
  cellchat_nonpCR
)

write.csv(
  comm_pCR,
  file.path(
    target_dir,
    "10_CellChat_communication_pCR.csv"
  ),
  row.names = FALSE
)

write.csv(
  comm_nonpCR,
  file.path(
    target_dir,
    "11_CellChat_communication_nonpCR.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 29. TARGETED LIGAND-RECEPTOR DIRECTIONS
# ============================================================

target_sources <- c(
  "IL1B_Mono",
  "CD8_T"
)

target_targets <- c(
  "IL1B_Mono",
  "CD8_T"
)

comm_pCR_target <- comm_pCR %>%
  dplyr::filter(
    source %in% target_sources,
    target %in% target_targets
  )

comm_nonpCR_target <- comm_nonpCR %>%
  dplyr::filter(
    source %in% target_sources,
    target %in% target_targets
  )

write.csv(
  comm_pCR_target,
  file.path(
    target_dir,
    "12_targeted_LR_pCR.csv"
  ),
  row.names = FALSE
)

write.csv(
  comm_nonpCR_target,
  file.path(
    target_dir,
    "13_targeted_LR_nonpCR.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 30. pCR VS NON-pCR COMPARISON
# ============================================================

comm_pCR_target$response_group <- "pCR"

comm_nonpCR_target$response_group <- "non-pCR"

LR_comparison <- dplyr::bind_rows(
  comm_pCR_target,
  comm_nonpCR_target
)

write.csv(
  LR_comparison,
  file.path(
    target_dir,
    "14_LR_comparison_pCR_vs_nonpCR.csv"
  ),
  row.names = FALSE
)

print(
  head(
    LR_comparison
  )
)


# ============================================================
# 31. TOP CANDIDATE LIGAND-RECEPTOR INTERACTIONS
# ============================================================

top_pCR <- comm_pCR_target %>%
  dplyr::arrange(
    dplyr::desc(prob)
  ) %>%
  dplyr::slice_head(
    n = 20
  )

top_nonpCR <- comm_nonpCR_target %>%
  dplyr::arrange(
    dplyr::desc(prob)
  ) %>%
  dplyr::slice_head(
    n = 20
  )

write.csv(
  top_pCR,
  file.path(
    target_dir,
    "15_top_LR_pCR.csv"
  ),
  row.names = FALSE
)

write.csv(
  top_nonpCR,
  file.path(
    target_dir,
    "16_top_LR_nonpCR.csv"
  ),
  row.names = FALSE
)

print(
  top_pCR
)

print(
  top_nonpCR
)


# ============================================================
# 32. CELLCHAT BUBBLE PLOT - pCR
# ============================================================

pdf(
  file.path(
    target_dir,
    "17_CellChat_bubble_pCR.pdf"
  ),
  width = 10,
  height = 7
)

print(
  netVisual_bubble(
    cellchat_pCR,
    sources.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    targets.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    remove.isolate = FALSE
  )
)

dev.off()


# ============================================================
# 33. CELLCHAT BUBBLE PLOT - NON-pCR
# ============================================================

pdf(
  file.path(
    target_dir,
    "18_CellChat_bubble_nonpCR.pdf"
  ),
  width = 10,
  height = 7
)

print(
  netVisual_bubble(
    cellchat_nonpCR,
    sources.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    targets.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    remove.isolate = FALSE
  )
)

dev.off()


# ============================================================
# 34. PATIENT-LEVEL TARGET ABUNDANCE
# ============================================================

patient_target <- analysis_data@meta.data %>%
  dplyr::filter(
    response %in% c(
      "pCR",
      "non-pCR"
    ),
    target_cell != "Other"
  ) %>%
  dplyr::count(
    patient,
    response,
    target_cell,
    name = "n_cells"
  )

write.csv(
  patient_target,
  file.path(
    target_dir,
    "19_patient_level_target_abundance.csv"
  ),
  row.names = FALSE
)

print(
  patient_target
)


# ============================================================
# 35. PATIENT-LEVEL TARGET PERCENTAGES
# ============================================================

patient_total <- analysis_data@meta.data %>%
  dplyr::filter(
    response %in% c(
      "pCR",
      "non-pCR"
    )
  ) %>%
  dplyr::count(
    patient,
    response,
    name = "total_cells"
  )

patient_target_percent <- patient_target %>%
  dplyr::left_join(
    patient_total,
    by = c(
      "patient",
      "response"
    )
  ) %>%
  dplyr::mutate(
    percent = 100 * n_cells / total_cells
  )

write.csv(
  patient_target_percent,
  file.path(
    target_dir,
    "20_patient_level_target_percentages.csv"
  ),
  row.names = FALSE
)

print(
  patient_target_percent
)


# ============================================================
# 36. PATIENT-LEVEL STATISTICS
# ============================================================

patient_stats <- patient_target_percent %>%
  dplyr::group_by(
    target_cell
  ) %>%
  dplyr::summarise(
    n_pCR = sum(
      response == "pCR"
    ),
    n_nonpCR = sum(
      response == "non-pCR"
    ),
    median_pCR = median(
      percent[
        response == "pCR"
      ],
      na.rm = TRUE
    ),
    median_nonpCR = median(
      percent[
        response == "non-pCR"
      ],
      na.rm = TRUE
    ),
    p_value = tryCatch(
      wilcox.test(
        percent ~ response
      )$p.value,
      error = function(e) {
        NA_real_
      }
    ),
    .groups = "drop"
  )

write.csv(
  patient_stats,
  file.path(
    target_dir,
    "21_patient_level_statistics.csv"
  ),
  row.names = FALSE
)

print(
  patient_stats
)


# ============================================================
# 37. FINAL SAVE
# ============================================================

saveRDS(
  analysis_data,
  file.path(
    target_dir,
    "GSE205506_targeted_analysis_final.rds"
  )
)

saveRDS(
  target_data,
  file.path(
    target_dir,
    "GSE205506_target_cells_final.rds"
  )
)

write.csv(
  analysis_data@meta.data,
  file.path(
    target_dir,
    "22_final_metadata.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 38. FINAL CHECK
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "TARGETED ANALYSIS COMPLETED\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Final cells:",
  ncol(analysis_data),
  "\n\n"
)

cat(
  "Target populations:\n"
)

print(
  table(
    analysis_data$target_cell,
    useNA = "ifany"
  )
)

cat(
  "\nResponse:\n"
)

print(
  table(
    analysis_data$response,
    useNA = "ifany"
  )
)

cat(
  "\nTarget populations x response:\n"
)

print(
  table(
    analysis_data$target_cell,
    analysis_data$response,
    useNA = "ifany"
  )
)

cat(
  "\nFinal object saved at:\n",
  file.path(
    target_dir,
    "GSE205506_targeted_analysis_final.rds"
  ),
  "\n"
)

cat(
  "\nAnalysis finished.\n"
)


cellchat_pCR <- run_cellchat(
  analysis_data,
  "pCR",
  file.path(target_dir, "CellChat_pCR.rds")
)

cellchat_nonpCR <- run_cellchat(
  analysis_data,
  "non-pCR",
  file.path(target_dir, "CellChat_nonpCR.rds")
)

















# ============================================================
# 25. CELLCHAT + TARGETED LR ANALYSIS + PATIENT-LEVEL VALIDATION
# ============================================================

message("============================================================")
message("STARTING CELLCHAT / TARGETED LR ANALYSIS")
message("============================================================")


# ------------------------------------------------------------
# 25.1 Define CellChat function
# ------------------------------------------------------------

run_cellchat <- function(seurat_obj, response_group, output_rds) {
  
  message("========================================")
  message("Running CellChat: ", response_group)
  message("========================================")
  
  # ----------------------------------------------------------
  # Subset response group
  # ----------------------------------------------------------
  
  obj <- subset(
    seurat_obj,
    subset = response == response_group
  )
  
  # ----------------------------------------------------------
  # Keep only target populations
  # ----------------------------------------------------------
  
  obj <- subset(
    obj,
    subset = target_cell %in% c("IL1B_Mono", "CD8_T")
  )
  
  obj$target_cell <- droplevels(obj$target_cell)
  
  message(
    "Target cells in ",
    response_group,
    ": ",
    ncol(obj)
  )
  
  print(table(obj$target_cell))
  
  # ----------------------------------------------------------
  # Explicit sample information
  # ----------------------------------------------------------
  
  if ("GSM" %in% colnames(obj@meta.data)) {
    
    obj$samples <- as.character(obj$GSM)
    
  } else if ("orig.ident" %in% colnames(obj@meta.data)) {
    
    obj$samples <- as.character(obj$orig.ident)
    
  } else {
    
    obj$samples <- "sample1"
    
  }
  
  # ----------------------------------------------------------
  # Join Seurat v5 RNA layers
  # ----------------------------------------------------------
  
  obj <- JoinLayers(
    obj,
    assay = "RNA"
  )
  
  # ----------------------------------------------------------
  # Create CellChat object
  # ----------------------------------------------------------
  
  cellchat <- createCellChat(
    object = obj,
    group.by = "target_cell",
    assay = "RNA"
  )
  
  # Explicitly restore sample information
  cellchat@meta$samples <- obj$samples
  
  # ----------------------------------------------------------
  # Human CellChat database
  # ----------------------------------------------------------
  
  cellchat@DB <- CellChatDB.human
  
  # ----------------------------------------------------------
  # Subset data
  # ----------------------------------------------------------
  
  cellchat <- subsetData(cellchat)
  
  # ----------------------------------------------------------
  # Identify overexpressed genes
  #
  # presto is installed, therefore use fast implementation
  # ----------------------------------------------------------
  
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    do.fast = TRUE
  )
  
  # ----------------------------------------------------------
  # Identify overexpressed ligand-receptor interactions
  # ----------------------------------------------------------
  
  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )
  
  # ----------------------------------------------------------
  # IMPORTANT:
  # projectData() removed because it is unavailable
  # in the installed CellChat version.
  # It is not required for this LR analysis.
  # ----------------------------------------------------------
  
  # ----------------------------------------------------------
  # Compute communication probability
  # ----------------------------------------------------------
  
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )
  
  # ----------------------------------------------------------
  # Remove interactions supported by too few cells
  # ----------------------------------------------------------
  
  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )
  
  # ----------------------------------------------------------
  # Compute communication pathways
  # ----------------------------------------------------------
  
  cellchat <- computeCommunProbPathway(
    cellchat
  )
  
  # ----------------------------------------------------------
  # Aggregate network
  # ----------------------------------------------------------
  
  cellchat <- aggregateNet(
    cellchat
  )
  
  # ----------------------------------------------------------
  # Save CellChat object
  # ----------------------------------------------------------
  
  saveRDS(
    cellchat,
    output_rds
  )
  
  message(
    "CellChat completed successfully: ",
    response_group
  )
  
  return(cellchat)
}


# ============================================================
# 26. CELLCHAT — pCR
# ============================================================

message("============================================================")
message("26. CELLCHAT: pCR")
message("============================================================")

cellchat_pCR <- run_cellchat(
  analysis_data,
  "pCR",
  file.path(
    target_dir,
    "CellChat_pCR.rds"
  )
)


# ============================================================
# 27. CELLCHAT — non-pCR
# ============================================================

message("============================================================")
message("27. CELLCHAT: non-pCR")
message("============================================================")

cellchat_nonpCR <- run_cellchat(
  analysis_data,
  "non-pCR",
  file.path(
    target_dir,
    "CellChat_nonpCR.rds"
  )
)


# ============================================================
# 28. EXTRACT COMMUNICATIONS
# ============================================================

message("============================================================")
message("28. EXTRACTING COMMUNICATIONS")
message("============================================================")

comm_pCR <- subsetCommunication(
  cellchat_pCR
)

comm_nonpCR <- subsetCommunication(
  cellchat_nonpCR
)

write.csv(
  comm_pCR,
  file.path(
    target_dir,
    "05_CellChat_all_interactions_pCR.csv"
  ),
  row.names = FALSE
)

write.csv(
  comm_nonpCR,
  file.path(
    target_dir,
    "06_CellChat_all_interactions_nonpCR.csv"
  ),
  row.names = FALSE
)

message(
  "pCR interactions: ",
  nrow(comm_pCR)
)

message(
  "non-pCR interactions: ",
  nrow(comm_nonpCR)
)


# ============================================================
# 29. TARGETED IL1B_Mono <-> CD8_T INTERACTIONS
# ============================================================

message("============================================================")
message("29. TARGETED IL1B_Mono <-> CD8_T INTERACTIONS")
message("============================================================")

target_sources <- c(
  "IL1B_Mono",
  "CD8_T"
)

target_targets <- c(
  "IL1B_Mono",
  "CD8_T"
)

comm_pCR_target <- comm_pCR[
  comm_pCR$source %in% target_sources &
    comm_pCR$target %in% target_targets,
  ,
  drop = FALSE
]

comm_nonpCR_target <- comm_nonpCR[
  comm_nonpCR$source %in% target_sources &
    comm_nonpCR$target %in% target_targets,
  ,
  drop = FALSE
]

write.csv(
  comm_pCR_target,
  file.path(
    target_dir,
    "07_targeted_LR_pCR.csv"
  ),
  row.names = FALSE
)

write.csv(
  comm_nonpCR_target,
  file.path(
    target_dir,
    "08_targeted_LR_nonpCR.csv"
  ),
  row.names = FALSE
)

message(
  "Targeted pCR interactions: ",
  nrow(comm_pCR_target)
)

message(
  "Targeted non-pCR interactions: ",
  nrow(comm_nonpCR_target)
)


# ============================================================
# 30. COMBINE pCR AND NON-pCR
# ============================================================

message("============================================================")
message("30. COMPARING pCR VS NON-pCR")
message("============================================================")

comm_pCR_target$response_group <- "pCR"

comm_nonpCR_target$response_group <- "non-pCR"

comm_target_all <- rbind(
  comm_pCR_target,
  comm_nonpCR_target
)

write.csv(
  comm_target_all,
  file.path(
    target_dir,
    "09_targeted_LR_pCR_vs_nonpCR.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 31. TOP TARGETED INTERACTIONS
# ============================================================

message("============================================================")
message("31. TOP TARGETED INTERACTIONS")
message("============================================================")

if (nrow(comm_target_all) > 0 &&
    "prob" %in% colnames(comm_target_all)) {
  
  top_targeted <- comm_target_all[
    order(
      comm_target_all$prob,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  top_targeted <- head(
    top_targeted,
    30
  )
  
  write.csv(
    top_targeted,
    file.path(
      target_dir,
      "10_top_targeted_interactions.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  top_targeted <- comm_target_all
  
  message(
    "No targeted interactions available for ranking."
  )
}


# ============================================================
# 32. BUBBLE PLOT — pCR
# ============================================================

message("============================================================")
message("32. BUBBLE PLOT — pCR")
message("============================================================")

if (nrow(comm_pCR_target) > 0) {
  
  p_bubble_pCR <- netVisual_bubble(
    cellchat_pCR,
    sources.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    targets.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    remove.isolate = FALSE
  )
  
  ggsave(
    file.path(
      target_dir,
      "11_CellChat_bubble_pCR.pdf"
    ),
    p_bubble_pCR,
    width = 10,
    height = 8
  )
  
  print(p_bubble_pCR)
  
} else {
  
  message(
    "No targeted pCR interactions available for bubble plot."
  )
}


# ============================================================
# 33. BUBBLE PLOT — non-pCR
# ============================================================

message("============================================================")
message("33. BUBBLE PLOT — non-pCR")
message("============================================================")

if (nrow(comm_nonpCR_target) > 0) {
  
  p_bubble_nonpCR <- netVisual_bubble(
    cellchat_nonpCR,
    sources.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    targets.use = c(
      "IL1B_Mono",
      "CD8_T"
    ),
    remove.isolate = FALSE
  )
  
  ggsave(
    file.path(
      target_dir,
      "12_CellChat_bubble_nonpCR.pdf"
    ),
    p_bubble_nonpCR,
    width = 10,
    height = 8
  )
  
  print(p_bubble_nonpCR)
  
} else {
  
  message(
    "No targeted non-pCR interactions available for bubble plot."
  )
}


# ============================================================
# 34. PATIENT-LEVEL TARGET CELL ABUNDANCE
# ============================================================

message("============================================================")
message("34. PATIENT-LEVEL TARGET CELL ABUNDANCE")
message("============================================================")

target_meta <- analysis_data@meta.data

target_meta <- target_meta[
  !is.na(target_meta$response) &
    target_meta$response %in% c(
      "pCR",
      "non-pCR"
    ) &
    target_meta$target_cell %in% c(
      "IL1B_Mono",
      "CD8_T"
    ),
  ,
  drop = FALSE
]

patient_target_counts <- as.data.frame(
  table(
    target_meta$patient,
    target_meta$target_cell
  )
)

colnames(patient_target_counts) <- c(
  "patient",
  "cell_type",
  "n_cells"
)

patient_response <- unique(
  target_meta[
    ,
    c(
      "patient",
      "response"
    ),
    drop = FALSE
  ]
)

patient_target_counts <- merge(
  patient_target_counts,
  patient_response,
  by = "patient",
  all.x = TRUE
)

write.csv(
  patient_target_counts,
  file.path(
    target_dir,
    "13_patient_level_target_cell_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 35. PATIENT-LEVEL TARGET CELL PERCENTAGES
# ============================================================

message("============================================================")
message("35. PATIENT-LEVEL TARGET CELL PERCENTAGES")
message("============================================================")

patient_total_cells <- as.data.frame(
  table(
    target_meta$patient
  )
)

colnames(patient_total_cells) <- c(
  "patient",
  "total_target_cells"
)

patient_target_percent <- merge(
  patient_target_counts,
  patient_total_cells,
  by = "patient",
  all.x = TRUE
)

patient_target_percent$percent_target <- (
  patient_target_percent$n_cells /
    patient_target_percent$total_target_cells
) * 100

write.csv(
  patient_target_percent,
  file.path(
    target_dir,
    "14_patient_level_target_percentages.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 36. PATIENT-LEVEL SUMMARY STATISTICS
# ============================================================

message("============================================================")
message("36. PATIENT-LEVEL SUMMARY STATISTICS")
message("============================================================")

if (nrow(patient_target_percent) > 0) {
  
  patient_stats <- do.call(
    rbind,
    lapply(
      unique(patient_target_percent$cell_type),
      function(ct) {
        
        dat <- patient_target_percent[
          patient_target_percent$cell_type == ct,
          ,
          drop = FALSE
        ]
        
        pcr_values <- dat$percent_target[
          dat$response == "pCR"
        ]
        
        nonpcr_values <- dat$percent_target[
          dat$response == "non-pCR"
        ]
        
        p_value <- NA_real_
        
        if (
          length(pcr_values) >= 2 &&
          length(nonpcr_values) >= 2
        ) {
          
          p_value <- tryCatch(
            wilcox.test(
              pcr_values,
              nonpcr_values,
              exact = FALSE
            )$p.value,
            error = function(e) NA_real_
          )
        }
        
        data.frame(
          cell_type = ct,
          n_pCR = length(pcr_values),
          n_nonpCR = length(nonpcr_values),
          median_pCR = ifelse(
            length(pcr_values) > 0,
            median(
              pcr_values,
              na.rm = TRUE
            ),
            NA
          ),
          median_nonpCR = ifelse(
            length(nonpcr_values) > 0,
            median(
              nonpcr_values,
              na.rm = TRUE
            ),
            NA
          ),
          p_value = p_value
        )
      }
    )
  )
  
  write.csv(
    patient_stats,
    file.path(
      target_dir,
      "15_patient_level_statistics.csv"
    ),
    row.names = FALSE
  )
  
  print(patient_stats)
  
} else {
  
  patient_stats <- data.frame()
  
  message(
    "No patient-level target data available."
  )
}


# ============================================================
# 37. SAVE FINAL ANALYSIS OBJECTS
# ============================================================

message("============================================================")
message("37. SAVING FINAL OBJECTS")
message("============================================================")

saveRDS(
  analysis_data,
  file.path(
    target_dir,
    "GSE205506_targeted_analysis_final.rds"
  )
)

saveRDS(
  cellchat_pCR,
  file.path(
    target_dir,
    "CellChat_pCR_final.rds"
  )
)

saveRDS(
  cellchat_nonpCR,
  file.path(
    target_dir,
    "CellChat_nonpCR_final.rds"
  )
)

saveRDS(
  comm_target_all,
  file.path(
    target_dir,
    "targeted_LR_interactions_final.rds"
  )
)

saveRDS(
  patient_target_percent,
  file.path(
    target_dir,
    "patient_level_target_percentages_final.rds"
  )
)


# ============================================================
# 38. FINAL CHECK
# ============================================================

message("============================================================")
message("38. FINAL CHECK")
message("============================================================")

cat("\n")
cat("============================================================\n")
cat("GSE205506 TARGETED ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "Total cells in analysis_data: ",
  ncol(analysis_data),
  "\n"
)

cat(
  "IL1B_Mono cells: ",
  sum(
    analysis_data$target_cell == "IL1B_Mono",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "CD8_T cells: ",
  sum(
    analysis_data$target_cell == "CD8_T",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Targeted pCR interactions: ",
  nrow(comm_pCR_target),
  "\n"
)

cat(
  "Targeted non-pCR interactions: ",
  nrow(comm_nonpCR_target),
  "\n"
)

cat(
  "Patient-level rows: ",
  nrow(patient_target_percent),
  "\n"
)

cat(
  "Output directory: ",
  target_dir,
  "\n"
)

cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")






# ============================================================
# SAVE EVERYTHING — FINAL RESULTS
# ============================================================

# Vérifier le dossier actuel
cat("TARGET DIRECTORY:\n")
print(target_dir)

# Créer un dossier de résultats final
final_results_dir <- file.path(
  target_dir,
  "FINAL_RESULTS"
)

if (!dir.exists(final_results_dir)) {
  dir.create(
    final_results_dir,
    recursive = TRUE
  )
}

cat("\nFINAL RESULTS DIRECTORY:\n")
print(final_results_dir)


# ------------------------------------------------------------
# 1. MAIN SEURAT OBJECT
# ------------------------------------------------------------

saveRDS(
  analysis_data,
  file.path(
    final_results_dir,
    "GSE205506_targeted_analysis_FINAL.rds"
  )
)


# ------------------------------------------------------------
# 2. CELLCHAT OBJECTS
# ------------------------------------------------------------

saveRDS(
  cellchat_pCR,
  file.path(
    final_results_dir,
    "CellChat_pCR_FINAL.rds"
  )
)

saveRDS(
  cellchat_nonpCR,
  file.path(
    final_results_dir,
    "CellChat_nonpCR_FINAL.rds"
  )
)


# ------------------------------------------------------------
# 3. TARGETED LR RESULTS
# ------------------------------------------------------------

saveRDS(
  comm_pCR_target,
  file.path(
    final_results_dir,
    "LR_targeted_pCR_FINAL.rds"
  )
)

saveRDS(
  comm_nonpCR_target,
  file.path(
    final_results_dir,
    "LR_targeted_nonpCR_FINAL.rds"
  )
)

saveRDS(
  comm_target_all,
  file.path(
    final_results_dir,
    "LR_targeted_pCR_vs_nonpCR_FINAL.rds"
  )
  
  
  # ------------------------------------------------------------
  # 4. CSV FILES
  # ------------------------------------------------------------
  
  write.csv(
    comm_pCR_target,
    file.path(
      final_results_dir,
      "LR_targeted_pCR_FINAL.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    comm_nonpCR_target,
    file.path(
      final_results_dir,
      "LR_targeted_nonpCR_FINAL.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    comm_target_all,
    file.path(
      final_results_dir,
      "LR_targeted_pCR_vs_nonpCR_FINAL.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    patient_target_counts,
    file.path(
      final_results_dir,
      "Patient_target_cell_counts.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    patient_target_percent,
    file.path(
      final_results_dir,
      "Patient_target_cell_percentages.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    patient_stats,
    file.path(
      final_results_dir,
      "Patient_level_statistics.csv"
    ),
    row.names = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 5. SAVE COMPLETE R SESSION
  # ------------------------------------------------------------
  
  save.image(
    file.path(
      final_results_dir,
      "GSE205506_COMPLETE_FINAL_SESSION.RData"
    )
  )
  
  
  # ------------------------------------------------------------
  # 6. LIST EVERYTHING THAT WAS SAVED
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("FILES SAVED SUCCESSFULLY\n")
  cat("============================================================\n")
  
  print(
    list.files(
      final_results_dir,
      full.names = TRUE
    )
  )
  
  cat("\n")
  cat("============================================================\n")
  cat("RESULTS LOCATION:\n")
  cat(final_results_dir)
  cat("\n============================================================\n")
  
  # ============================================================
  # SAVE FINAL RESULTS — SAFE VERSION
  # ============================================================
  
  final_report_dir <- file.path(
    target_dir,
    "FINAL_REPORT"
  )
  
  if (!dir.exists(final_report_dir)) {
    dir.create(
      final_report_dir,
      recursive = TRUE
    )
  }
  
  # 1. Save IL1B-CD8 interactions
  write.csv(
    comm_pCR_target,
    file.path(
      final_report_dir,
      "IL1B_CD8_interactions_pCR.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    comm_nonpCR_target,
    file.path(
      final_report_dir,
      "IL1B_CD8_interactions_nonpCR.csv"
    ),
    row.names = FALSE
  )
  
  # 2. Save patient-level results if they exist
  if (exists("patient_target_counts")) {
    write.csv(
      patient_target_counts,
      file.path(
        final_report_dir,
        "Patient_target_cell_counts.csv"
      ),
      row.names = FALSE
    )
  }
  
  if (exists("patient_target_percent")) {
    write.csv(
      patient_target_percent,
      file.path(
        final_report_dir,
        "Patient_target_cell_percentages.csv"
      ),
      row.names = FALSE
    )
  }
  
  if (exists("patient_stats")) {
    write.csv(
      patient_stats,
      file.path(
        final_report_dir,
        "Patient_level_statistics.csv"
      ),
      row.names = FALSE
    )
  }
  
  # 3. Save CellChat objects
  saveRDS(
    cellchat_pCR,
    file.path(
      final_report_dir,
      "CellChat_pCR_FINAL.rds"
    )
  )
  
  saveRDS(
    cellchat_nonpCR,
    file.path(
      final_report_dir,
      "CellChat_nonpCR_FINAL.rds"
    )
  )
  
  # 4. Save Seurat object
  saveRDS(
    analysis_data,
    file.path(
      final_report_dir,
      "GSE205506_analysis_FINAL.rds"
    )
  )
  
  # 5. Save complete R session
  save.image(
    file.path(
      final_report_dir,
      "GSE205506_COMPLETE_SESSION.RData"
    )
  )
  
  # 6. Show exactly where everything was saved
  cat("\n============================================\n")
  cat("SAUVEGARDE TERMINEE\n")
  cat("============================================\n")
  cat("\nDossier :\n")
  cat(normalizePath(final_report_dir))
  cat("\n\nFichiers :\n")
  
  print(
    list.files(
      final_report_dir
    )
  )
  
  cat("\n============================================\n")
  
  # ============================================================
  # SAVE FINAL RESULTS — SAFE VERSION
  # ============================================================
  
  final_report_dir <- file.path(
    target_dir,
    "FINAL_REPORT"
  )
  
  if (!dir.exists(final_report_dir)) {
    dir.create(
      final_report_dir,
      recursive = TRUE
    )
  }
  
  # 1. Save IL1B-CD8 interactions
  write.csv(
    comm_pCR_target,
    file.path(
      final_report_dir,
      "IL1B_CD8_interactions_pCR.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    comm_nonpCR_target,
    file.path(
      final_report_dir,
      "IL1B_CD8_interactions_nonpCR.csv"
    ),
    row.names = FALSE
  )
  
  # 2. Save patient-level results if they exist
  if (exists("patient_target_counts")) {
    write.csv(
      patient_target_counts,
      file.path(
        final_report_dir,
        "Patient_target_cell_counts.csv"
      ),
      row.names = FALSE
    )
  }
  
  if (exists("patient_target_percent")) {
    write.csv(
      patient_target_percent,
      file.path(
        final_report_dir,
        "Patient_target_cell_percentages.csv"
      ),
      row.names = FALSE
    )
  }
  
  if (exists("patient_stats")) {
    write.csv(
      patient_stats,
      file.path(
        final_report_dir,
        "Patient_level_statistics.csv"
      ),
      row.names = FALSE
    )
  }
  
  # 3. Save CellChat objects
  saveRDS(
    cellchat_pCR,
    file.path(
      final_report_dir,
      "CellChat_pCR_FINAL.rds"
    )
  )
  
  saveRDS(
    cellchat_nonpCR,
    file.path(
      final_report_dir,
      "CellChat_nonpCR_FINAL.rds"
    )
  )
  
  # 4. Save Seurat object
  saveRDS(
    analysis_data,
    file.path(
      final_report_dir,
      "GSE205506_analysis_FINAL.rds"
    )
  )
  
  # 5. Save complete R session
  save.image(
    file.path(
      final_report_dir,
      "GSE205506_COMPLETE_SESSION.RData"
    )
  )
  
  # 6. Show exactly where everything was saved
  cat("\n============================================\n")
  cat("SAUVEGARDE TERMINEE\n")
  cat("============================================\n")
  cat("\nDossier :\n")
  cat(normalizePath(final_report_dir))
  cat("\n\nFichiers :\n")
  
  print(
    list.files(
      final_report_dir
    )
  )
  
  cat("\n============================================\n")
  
  
  # ============================================================
  # FINAL REPORT PACKAGE
  # GSE205506
  #
  # QUESTION:
  # Are there ligand-receptor interactions between
  # IL1B+ monocytes and CD8+ T cells, and do they differ
  # between pCR and non-pCR?
  #
  # THIS SCRIPT:
  #   1. Extracts IL1B_Mono <-> CD8_T interactions
  #   2. Creates the main result tables
  #   3. Creates clear report figures
  #   4. Saves everything in one FINAL_REPORT folder
  #
  # NO CELLCHAT RECALCULATION
  # ============================================================
  
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  
  # ============================================================
  # 1. CREATE FINAL REPORT DIRECTORY
  # ============================================================
  
  final_report_dir <- file.path(
    target_dir,
    "FINAL_REPORT"
  )
  
  if (!dir.exists(final_report_dir)) {
    dir.create(
      final_report_dir,
      recursive = TRUE
    )
  }
  
  cat("\n")
  cat("============================================================\n")
  cat("FINAL REPORT DIRECTORY\n")
  cat(final_report_dir)
  cat("\n============================================================\n")
  
  
  # ============================================================
  # 2. CHECK MAIN OBJECTS
  # ============================================================
  
  if (!exists("analysis_data")) {
    stop("analysis_data is not available.")
  }
  
  if (!exists("comm_pCR_target")) {
    stop("comm_pCR_target is not available.")
  }
  
  if (!exists("comm_nonpCR_target")) {
    stop("comm_nonpCR_target is not available.")
  }
  
  
  # ============================================================
  # 3. EXTRACT IL1B_Mono <-> CD8_T INTERACTIONS
  # ============================================================
  
  cat("\nExtracting IL1B_Mono <-> CD8_T interactions...\n")
  
  
  pCR_cross <- comm_pCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (
        source == "IL1B_Mono" &
          target == "CD8_T"
      ) |
        (
          source == "CD8_T" &
            target == "IL1B_Mono"
        )
    ) %>%
    mutate(
      response = "pCR"
    )
  
  
  nonpCR_cross <- comm_nonpCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (
        source == "IL1B_Mono" &
          target == "CD8_T"
      ) |
        (
          source == "CD8_T" &
            target == "IL1B_Mono"
        )
    ) %>%
    mutate(
      response = "non-pCR"
    )
  
  
  # ============================================================
  # 4. CREATE ROBUST LR NAME
  # ============================================================
  
  make_LR_name <- function(x) {
    
    if (
      "interaction_name" %in%
      colnames(x)
    ) {
      
      name <- as.character(
        x$interaction_name
      )
      
    } else if (
      all(
        c(
          "ligand",
          "receptor"
        ) %in%
        colnames(x)
      )
    ) {
      
      name <- paste(
        x$ligand,
        x$receptor,
        sep = " - "
      )
      
    } else {
      
      name <- paste(
        x$source,
        x$target,
        sep = " -> "
      )
    }
    
    name
  }
  
  
  pCR_cross$LR <- make_LR_name(
    pCR_cross
  )
  
  nonpCR_cross$LR <- make_LR_name(
    nonpCR_cross
  )
  
  
  # ============================================================
  # 5. CREATE MASTER TABLE
  # ============================================================
  
  LR_results <- bind_rows(
    pCR_cross,
    nonpCR_cross
  )
  
  
  # Keep useful columns
  
  preferred_columns <- c(
    "response",
    "source",
    "target",
    "LR",
    "ligand",
    "receptor",
    "interaction_name",
    "interaction_name_2",
    "pathway_name",
    "prob",
    "pval"
  )
  
  available_columns <- intersect(
    preferred_columns,
    colnames(LR_results)
  )
  
  LR_results <- LR_results[
    ,
    available_columns,
    drop = FALSE
  ]
  
  
  # Order
  
  if ("prob" %in% colnames(LR_results)) {
    
    LR_results <- LR_results %>%
      arrange(
        LR,
        response,
        desc(prob)
      )
  }
  
  
  # ============================================================
  # 6. SAVE COMPLETE IL1B-CD8 TABLE
  # ============================================================
  
  write.csv(
    LR_results,
    file.path(
      final_report_dir,
      "TABLE_1_IL1B_CD8_ALL_INTERACTIONS.csv"
    ),
    row.names = FALSE
  )
  
  
  # ============================================================
  # 7. CREATE WIDE COMPARISON TABLE
  # ============================================================
  
  if ("prob" %in% colnames(LR_results)) {
    
    LR_probability <- LR_results %>%
      select(
        response,
        source,
        target,
        LR,
        prob
      ) %>%
      group_by(
        response,
        source,
        target,
        LR
      ) %>%
      summarise(
        probability = max(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = response,
        values_from = probability,
        values_fill = 0
      )
    
    # Make sure both columns exist
    
    if (!"pCR" %in% colnames(LR_probability)) {
      LR_probability$pCR <- 0
    }
    
    if (!"non-pCR" %in% colnames(LR_probability)) {
      LR_probability$`non-pCR` <- 0
    }
    
    LR_probability <- LR_probability %>%
      mutate(
        difference_pCR_minus_nonpCR =
          pCR - `non-pCR`
      ) %>%
      arrange(
        desc(
          abs(
            difference_pCR_minus_nonpCR
          )
        )
      )
    
    
    write.csv(
      LR_probability,
      file.path(
        final_report_dir,
        "TABLE_2_IL1B_CD8_pCR_vs_nonpCR.csv"
      ),
      row.names = FALSE
    )
    
  } else {
    
    LR_probability <- NULL
  }
  
  
  # ============================================================
  # 8. PRINT THE ACTUAL ANSWER TO THE BIOLOGICAL QUESTION
  # ============================================================
  
  cat("\n\n")
  cat("============================================================\n")
  cat("BIOLOGICAL QUESTION: IL1B+ MONOCYTES <-> CD8+ T CELLS\n")
  cat("============================================================\n\n")
  
  
  if (nrow(LR_results) == 0) {
    
    cat(
      "NO IL1B_Mono <-> CD8_T interactions were detected.\n"
    )
    
  } else {
    
    cat(
      "INTERACTIONS DETECTED:\n\n"
    )
    
    print(
      LR_results,
      row.names = FALSE
    )
  }
  
  
  # ============================================================
  # FIGURE 1
  # TARGET POPULATION IDENTIFICATION
  # ============================================================
  
  cat("\nCreating Figure 1...\n")
  
  
  target_reduction <- if (
    "umap_before" %in%
    Reductions(analysis_data)
  ) {
    
    "umap_before"
    
  } else if (
    "umap" %in%
    Reductions(analysis_data)
  ) {
    
    "umap"
    
  } else {
    
    stop(
      "No UMAP reduction found."
    )
  }
  
  
  p_fig1 <- DimPlot(
    analysis_data,
    reduction = target_reduction,
    group.by = "target_cell"
  ) +
    ggtitle(
      "Identification of IL1B+ monocytes and CD8+ T cells"
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      legend.title = element_blank()
    )
  
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_1_Target_populations_UMAP.pdf"
    ),
    p_fig1,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_1_Target_populations_UMAP.png"
    ),
    p_fig1,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # FIGURE 1B
  # MARKER VALIDATION
  # ============================================================
  
  marker_genes <- c(
    "IL1B",
    "LST1",
    "S100A8",
    "FCN1",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD8B"
  )
  
  marker_genes <- marker_genes[
    marker_genes %in%
      rownames(analysis_data)
  ]
  
  
  p_fig1b <- DotPlot(
    analysis_data,
    features = marker_genes,
    group.by = "target_cell"
  ) +
    RotatedAxis() +
    ggtitle(
      "Marker expression supporting target-cell annotation"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13
      )
    )
  
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_1B_Target_marker_validation.pdf"
    ),
    p_fig1b,
    width = 10,
    height = 6
  )
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_1B_Target_marker_validation.png"
    ),
    p_fig1b,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # FIGURE 2
  # PATIENT-LEVEL TARGET ABUNDANCE
  # ============================================================
  
  cat("\nCreating Figure 2...\n")
  
  
  meta <- analysis_data@meta.data
  
  
  # Keep response-labeled cells
  
  meta_response <- meta[
    !is.na(meta$response) &
      meta$response %in% c(
        "pCR",
        "non-pCR"
      ),
    ,
    drop = FALSE
  ]
  
  
  # If sample_status exists, keep tumor cells
  
  if ("sample_status" %in% colnames(meta_response)) {
    
    meta_response <- meta_response[
      meta_response$sample_status == "tumor",
      ,
      drop = FALSE
    ]
  }
  
  
  # Total cells per patient
  
  total_cells <- meta_response %>%
    group_by(patient, response) %>%
    summarise(
      total_cells = n(),
      .groups = "drop"
    )
  
  
  # Target cells
  
  target_meta <- meta_response[
    meta_response$target_cell %in%
      c(
        "IL1B_Mono",
        "CD8_T"
      ),
    ,
    drop = FALSE
  ]
  
  
  target_counts <- target_meta %>%
    group_by(
      patient,
      response,
      target_cell
    ) %>%
    summarise(
      target_cells = n(),
      .groups = "drop"
    )
  
  
  # Complete patient x population grid
  
  patients <- unique(
    meta_response$patient
  )
  
  responses <- unique(
    meta_response$response
  )
  
  populations <- c(
    "IL1B_Mono",
    "CD8_T"
  )
  
  complete_grid <- expand.grid(
    patient = patients,
    target_cell = populations,
    stringsAsFactors = FALSE
  )
  
  
  patient_response <- meta_response %>%
    select(
      patient,
      response
    ) %>%
    distinct()
  
  
  abundance_table <- complete_grid %>%
    left_join(
      patient_response,
      by = "patient"
    ) %>%
    left_join(
      total_cells,
      by = c(
        "patient",
        "response"
      )
    ) %>%
    left_join(
      target_counts,
      by = c(
        "patient",
        "response",
        "target_cell"
      )
    ) %>%
    mutate(
      target_cells = ifelse(
        is.na(target_cells),
        0,
        target_cells
      ),
      percent_target =
        100 *
        target_cells /
        total_cells
    )
  
  
  abundance_table$cell_type <- recode(
    abundance_table$target_cell,
    "IL1B_Mono" = "IL1B+ monocytes",
    "CD8_T" = "CD8+ T cells"
  )
  
  
  write.csv(
    abundance_table,
    file.path(
      final_report_dir,
      "TABLE_3_Patient_level_target_abundance.csv"
    ),
    row.names = FALSE
  )
  
  
  p_fig2 <- ggplot(
    abundance_table,
    aes(
      x = response,
      y = percent_target
    )
  ) +
    
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.10,
      size = 2.5,
      alpha = 0.85
    ) +
    
    facet_wrap(
      ~ cell_type,
      scales = "free_y"
    ) +
    
    labs(
      title = "Target-cell abundance according to clinical response",
      x = "Response",
      y = "Target cells (% of response-labeled tumor cells)"
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
  
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_2_Patient_level_target_abundance.pdf"
    ),
    p_fig2,
    width = 9,
    height = 6
  )
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_2_Patient_level_target_abundance.png"
    ),
    p_fig2,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # FIGURE 3
  # IL1B <-> CD8 LIGAND-RECEPTOR INTERACTIONS
  # ============================================================
  
  cat("\nCreating Figure 3...\n")
  
  
  if (
    nrow(LR_results) > 0 &&
    "prob" %in% colnames(LR_results)
  ) {
    
    # Select up to 20 strongest interactions
    # based on maximum communication probability
    
    top_LR <- LR_results %>%
      group_by(LR) %>%
      summarise(
        max_prob = max(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(
        desc(max_prob)
      ) %>%
      slice_head(
        n = 20
      ) %>%
      pull(LR)
    
    
    lr_plot <- LR_results %>%
      filter(
        LR %in% top_LR
      )
    
    
    lr_plot$LR <- factor(
      lr_plot$LR,
      levels = rev(top_LR)
    )
    
    
    p_fig3 <- ggplot(
      lr_plot,
      aes(
        x = response,
        y = LR,
        size = prob
      )
    ) +
      
      geom_point(
        alpha = 0.85
      ) +
      
      facet_grid(
        source ~ target,
        scales = "free_y",
        space = "free_y"
      ) +
      
      labs(
        title = "Ligand–receptor interactions between IL1B+ monocytes and CD8+ T cells",
        x = "Clinical response",
        y = "Ligand–receptor pair",
        size = "CellChat communication probability"
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        axis.text.y = element_text(
          size = 8
        ),
        strip.text = element_text(
          face = "bold"
        )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_3_IL1B_CD8_LR_interactions.pdf"
      ),
      p_fig3,
      width = 11,
      height = 9
    )
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_3_IL1B_CD8_LR_interactions.png"
      ),
      p_fig3,
      width = 11,
      height = 9,
      dpi = 300
    )
    
  } else {
    
    cat(
      "No IL1B-CD8 interactions available for Figure 3.\n"
    )
  }
  
  
  # ============================================================
  # FIGURE 4
  # DIFFERENCE IN COMMUNICATION PROBABILITY
  # ============================================================
  
  cat("\nCreating Figure 4...\n")
  
  
  if (
    !is.null(LR_probability) &&
    nrow(LR_probability) > 0
  ) {
    
    # Select strongest differences
    
    plot_diff <- LR_probability %>%
      arrange(
        desc(
          abs(
            difference_pCR_minus_nonpCR
          )
        )
      ) %>%
      slice_head(
        n = 15
      )
    
    
    plot_diff$LR <- factor(
      plot_diff$LR,
      levels = rev(
        plot_diff$LR
      )
    )
    
    
    p_fig4 <- ggplot(
      plot_diff,
      aes(
        x = difference_pCR_minus_nonpCR,
        y = LR
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      
      geom_col(
        aes(
          fill = source
        )
      ) +
      
      labs(
        title = "Difference in modeled communication probability",
        subtitle = "pCR minus non-pCR",
        x = "Communication probability difference",
        y = "Ligand–receptor pair",
        fill = "Source cell"
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        axis.text.y = element_text(
          size = 8
        )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_4_LR_probability_difference_pCR_vs_nonpCR.pdf"
      ),
      p_fig4,
      width = 10,
      height = 8
    )
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_4_LR_probability_difference_pCR_vs_nonpCR.png"
      ),
      p_fig4,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  
  # ============================================================
  # 9. SAVE ALL IMPORTANT R OBJECTS
  # ============================================================
  
  saveRDS(
    LR_results,
    file.path(
      final_report_dir,
      "RESULTS_IL1B_CD8_ALL.rds"
    )
  )
  
  if (!is.null(LR_probability)) {
    
    saveRDS(
      LR_probability,
      file.path(
        final_report_dir,
        "RESULTS_IL1B_CD8_pCR_vs_nonpCR.rds"
      )
    )
  }
  
  
  saveRDS(
    abundance_table,
    file.path(
      final_report_dir,
      "RESULTS_patient_level_abundance.rds"
    )
  )
  
  
  # ============================================================
  # 10. SAVE COMPLETE SESSION
  # ============================================================
  
  save.image(
    file.path(
      final_report_dir,
      "GSE205506_FINAL_COMPLETE_SESSION.RData"
    )
  )
  
  
  # ============================================================
  # 11. FINAL FILE LIST
  # ============================================================
  
  cat("\n\n")
  cat("============================================================\n")
  cat("                 FINAL REPORT COMPLETED\n")
  cat("============================================================\n\n")
  
  cat(
    "Results folder:\n",
    final_report_dir,
    "\n\n"
  )
  
  cat("Files created:\n\n")
  
  print(
    list.files(
      final_report_dir,
      full.names = FALSE
    )
  )
  
  cat("\n============================================================\n")
  cat("IMPORTANT BIOLOGICAL RESULT TABLE\n")
  cat("============================================================\n\n")
  
  if (nrow(LR_results) == 0) {
    
    cat(
      "No IL1B_Mono <-> CD8_T ligand-receptor interactions were detected.\n"
    )
    
  } else {
    
    print(
      LR_results,
      row.names = FALSE
    )
  }
  
  cat("\n============================================================\n")
  cat("DONE — YOU CAN NOW CLOSE R\n")
  cat("============================================================\n")
  
  # ============================================================
  # GSE205506 — FINAL BIOLOGICAL REPORT
  # ============================================================
  #
  # Biological question:
  # Which ligand–receptor interactions between IL1B+ monocytes
  # and CD8+ T cells are associated with response to PD-1
  # blockade in dMMR/MSI-H colorectal cancer?
  #
  # This block:
  #   1. Identifies target populations
  #   2. Quantifies target populations at patient level
  #   3. Extracts IL1B_Mono <-> CD8_T CellChat interactions
  #   4. Compares modeled communication between pCR and non-pCR
  #   5. Creates final tables and publication-quality figures
  #   6. Saves all final outputs
  #
  # NO QC / clustering / Harmony / CellChat recomputation
  # ============================================================
  
  
  # ============================================================
  # 0. PACKAGES
  # ============================================================
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  
  # ============================================================
  # 1. FINAL OUTPUT DIRECTORY
  # ============================================================
  
  final_report_dir <- file.path(
    target_dir,
    "FINAL_REPORT"
  )
  
  if (!dir.exists(final_report_dir)) {
    dir.create(
      final_report_dir,
      recursive = TRUE
    )
  }
  
  cat("\n============================================\n")
  cat("FINAL REPORT DIRECTORY\n")
  cat("============================================\n")
  
  print(
    normalizePath(final_report_dir)
  )
  
  
  # ============================================================
  # 2. CHECK REQUIRED OBJECTS
  # ============================================================
  
  required_objects <- c(
    "analysis_data",
    "cellchat_pCR",
    "cellchat_nonpCR",
    "comm_pCR_target",
    "comm_nonpCR_target"
  )
  
  missing_objects <- required_objects[
    !sapply(
      required_objects,
      exists
    )
  ]
  
  if (length(missing_objects) > 0) {
    
    stop(
      paste(
        "Missing required objects:",
        paste(
          missing_objects,
          collapse = ", "
        )
      )
    )
  }
  
  
  # ============================================================
  # FIGURE 1
  # TARGET POPULATION IDENTIFICATION
  # ============================================================
  
  target_reduction <- if (
    "umap_before" %in% Reductions(analysis_data)
  ) {
    
    "umap_before"
    
  } else if (
    "umap" %in% Reductions(analysis_data)
  ) {
    
    "umap"
    
  } else {
    
    stop(
      "No UMAP reduction found in analysis_data."
    )
  }
  
  
  p_target_umap <- DimPlot(
    analysis_data,
    reduction = target_reduction,
    group.by = "target_cell"
  ) +
    ggtitle(
      "Target populations: IL1B+ monocytes and CD8+ T cells"
    ) +
    labs(
      subtitle =
        "GSE205506 dMMR/MSI-H colorectal cancer"
    ) +
    theme_classic() +
    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 14
        ),
      plot.subtitle =
        element_text(
          size = 10
        )
    )
  
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_1_TARGET_POPULATIONS_UMAP.pdf"
    ),
    p_target_umap,
    width = 9,
    height = 7
  )
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_1_TARGET_POPULATIONS_UMAP.png"
    ),
    p_target_umap,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  
  # ============================================================
  # FIGURE 2 + TABLE 1
  # PATIENT-LEVEL TARGET POPULATION ABUNDANCE
  # ============================================================
  
  # Keep only post-treatment samples with known response
  post_treatment <- analysis_data@meta.data %>%
    mutate(
      patient = as.character(patient),
      response = as.character(response),
      target_cell = as.character(target_cell)
    ) %>%
    filter(
      response %in% c(
        "pCR",
        "non-pCR"
      )
    )
  
  
  # Number of cells per patient
  patient_total <- post_treatment %>%
    group_by(
      patient,
      response
    ) %>%
    summarise(
      total_cells = n(),
      .groups = "drop"
    )
  
  
  # Target population counts
  patient_target_counts <- post_treatment %>%
    filter(
      target_cell %in% c(
        "IL1B_Mono",
        "CD8_T"
      )
    ) %>%
    count(
      patient,
      response,
      target_cell,
      name = "target_cells"
    )
  
  
  # Add patients with zero target cells
  patient_target_counts <- patient_total %>%
    tidyr::crossing(
      target_cell = c(
        "IL1B_Mono",
        "CD8_T"
      )
    ) %>%
    left_join(
      patient_target_counts,
      by = c(
        "patient",
        "response",
        "target_cell"
      )
    ) %>%
    mutate(
      target_cells =
        ifelse(
          is.na(target_cells),
          0,
          target_cells
        )
    )
  
  
  # Percentage among ALL response-labelled post-treatment cells
  patient_target_percentages <- patient_target_counts %>%
    mutate(
      percentage =
        100 *
        target_cells /
        total_cells
    )
  
  
  # Save patient-level tables
  
  write.csv(
    patient_target_counts,
    file.path(
      final_report_dir,
      "TABLE_1_PATIENT_TARGET_CELL_COUNTS.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    patient_target_percentages,
    file.path(
      final_report_dir,
      "TABLE_2_PATIENT_TARGET_CELL_PERCENTAGES.csv"
    ),
    row.names = FALSE
  )
  
  
  # Statistical summary
  patient_level_statistics <- patient_target_percentages %>%
    group_by(
      target_cell,
      response
    ) %>%
    summarise(
      n_patients = n(),
      mean_percentage =
        mean(
          percentage,
          na.rm = TRUE
        ),
      median_percentage =
        median(
          percentage,
          na.rm = TRUE
        ),
      SD =
        sd(
          percentage,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
  
  
  write.csv(
    patient_level_statistics,
    file.path(
      final_report_dir,
      "TABLE_3_PATIENT_LEVEL_SUMMARY.csv"
    ),
    row.names = FALSE
  )
  
  
  # Figure 2
  
  p_abundance <- ggplot(
    patient_target_percentages,
    aes(
      x = response,
      y = percentage
    )
  ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.12,
      size = 2.5,
      alpha = 0.8
    ) +
    facet_wrap(
      ~ target_cell,
      scales = "free_y"
    ) +
    labs(
      title =
        "Patient-level abundance of target populations",
      subtitle =
        "Post-treatment tumor cells",
      x = "Clinical response",
      y =
        "Target population (% of post-treatment cells)"
    ) +
    theme_classic() +
    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 14
        ),
      strip.text =
        element_text(
          face = "bold"
        )
    )
  
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_2_PATIENT_LEVEL_TARGET_ABUNDANCE.pdf"
    ),
    p_abundance,
    width = 10,
    height = 6
  )
  
  ggsave(
    file.path(
      final_report_dir,
      "FIGURE_2_PATIENT_LEVEL_TARGET_ABUNDANCE.png"
    ),
    p_abundance,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # CELLCHAT — EXTRACT ONLY IL1B_Mono <-> CD8_T
  # ============================================================
  
  pCR_IL1B_CD8 <- comm_pCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (source == "IL1B_Mono" &
         target == "CD8_T") |
        (source == "CD8_T" &
           target == "IL1B_Mono")
    ) %>%
    mutate(
      response = "pCR"
    )
  
  
  nonpCR_IL1B_CD8 <- comm_nonpCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (source == "IL1B_Mono" &
         target == "CD8_T") |
        (source == "CD8_T" &
           target == "IL1B_Mono")
    ) %>%
    mutate(
      response = "non-pCR"
    )
  
  
  # ============================================================
  # CREATE LR LABEL
  # ============================================================
  
  if (
    "interaction_name" %in%
    colnames(pCR_IL1B_CD8)
  ) {
    
    pCR_IL1B_CD8$LR <-
      as.character(
        pCR_IL1B_CD8$interaction_name
      )
    
    nonpCR_IL1B_CD8$LR <-
      as.character(
        nonpCR_IL1B_CD8$interaction_name
      )
    
  } else if (
    all(
      c(
        "ligand",
        "receptor"
      ) %in%
      colnames(pCR_IL1B_CD8)
    )
  ) {
    
    pCR_IL1B_CD8$LR <-
      paste(
        pCR_IL1B_CD8$ligand,
        pCR_IL1B_CD8$receptor,
        sep = " - "
      )
    
    nonpCR_IL1B_CD8$LR <-
      paste(
        nonpCR_IL1B_CD8$ligand,
        nonpCR_IL1B_CD8$receptor,
        sep = " - "
      )
    
  } else {
    
    pCR_IL1B_CD8$LR <-
      paste(
        pCR_IL1B_CD8$source,
        pCR_IL1B_CD8$target,
        sep = " -> "
      )
    
    nonpCR_IL1B_CD8$LR <-
      paste(
        nonpCR_IL1B_CD8$source,
        nonpCR_IL1B_CD8$target,
        sep = " -> "
      )
  }
  
  
  # ============================================================
  # TABLE 4
  # ALL IL1B/CD8 INTERACTIONS
  # ============================================================
  
  LR_results <- bind_rows(
    pCR_IL1B_CD8,
    nonpCR_IL1B_CD8
  )
  
  
  write.csv(
    LR_results,
    file.path(
      final_report_dir,
      "TABLE_4_IL1B_CD8_ALL_INTERACTIONS.csv"
    ),
    row.names = FALSE
  )
  
  
  # ============================================================
  # FIGURE 3 + TABLE 5
  # MODELED LR COMMUNICATION pCR vs non-pCR
  # ============================================================
  
  if (
    nrow(LR_results) > 0 &&
    "prob" %in%
    colnames(LR_results)
  ) {
    
    # ------------------------------------------
    # Comparison table
    # ------------------------------------------
    
    LR_comparison <- LR_results %>%
      group_by(
        response,
        source,
        target,
        LR
      ) %>%
      summarise(
        probability =
          max(
            prob,
            na.rm = TRUE
          ),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = response,
        values_from = probability,
        values_fill = 0
      )
    
    
    if (
      !"pCR" %in%
      colnames(LR_comparison)
    ) {
      LR_comparison$pCR <- 0
    }
    
    
    if (
      !"non-pCR" %in%
      colnames(LR_comparison)
    ) {
      LR_comparison$`non-pCR` <- 0
    }
    
    
    LR_comparison <-
      LR_comparison %>%
      mutate(
        difference_pCR_minus_nonpCR =
          pCR - `non-pCR`
      ) %>%
      arrange(
        desc(
          abs(
            difference_pCR_minus_nonpCR
          )
        )
      )
    
    
    write.csv(
      LR_comparison,
      file.path(
        final_report_dir,
        "TABLE_5_IL1B_CD8_pCR_vs_nonpCR.csv"
      ),
      row.names = FALSE
    )
    
    
    # ------------------------------------------
    # Figure 3
    # ------------------------------------------
    
    plot_data <-
      LR_results %>%
      mutate(
        interaction_id =
          paste(
            source,
            target,
            LR,
            sep = " | "
          )
      )
    
    
    top_ids <-
      plot_data %>%
      group_by(
        interaction_id
      ) %>%
      summarise(
        max_prob =
          max(
            prob,
            na.rm = TRUE
          ),
        .groups = "drop"
      ) %>%
      arrange(
        desc(max_prob)
      ) %>%
      slice_head(
        n = 20
      ) %>%
      pull(
        interaction_id
      )
    
    
    plot_data <-
      plot_data %>%
      filter(
        interaction_id %in%
          top_ids
      )
    
    
    label_table <-
      plot_data %>%
      distinct(
        interaction_id,
        LR
      )
    
    
    plot_data <-
      plot_data %>%
      left_join(
        label_table,
        by = "interaction_id",
        suffix =
          c(
            "",
            "_label"
          )
      )
    
    
    # Unique factor levels
    plot_data$interaction_id <-
      factor(
        plot_data$interaction_id,
        levels =
          rev(
            unique(
              plot_data$interaction_id
            )
          )
      )
    
    
    p_LR <-
      ggplot(
        plot_data,
        aes(
          x = response,
          y = interaction_id,
          size = prob
        )
      ) +
      geom_point(
        alpha = 0.85
      ) +
      facet_grid(
        source ~ target,
        scales = "free_y",
        space = "free_y"
      ) +
      scale_y_discrete(
        labels =
          function(x) {
            label_table$LR[
              match(
                x,
                label_table$interaction_id
              )
            ]
          }
      ) +
      labs(
        title =
          "IL1B+ monocyte–CD8+ T-cell ligand–receptor interactions",
        subtitle =
          "CellChat-modeled communication probability",
        x =
          "Clinical response",
        y =
          "Ligand–receptor pair",
        size =
          "Communication probability"
      ) +
      theme_classic() +
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),
        axis.text.y =
          element_text(
            size = 8
          ),
        strip.text =
          element_text(
            face = "bold"
          )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_3_IL1B_CD8_LR_INTERACTIONS.pdf"
      ),
      p_LR,
      width = 11,
      height = 9
    )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_3_IL1B_CD8_LR_INTERACTIONS.png"
      ),
      p_LR,
      width = 11,
      height = 9,
      dpi = 300
    )
    
    
    # ========================================================
    # FIGURE 4
    # pCR MINUS non-pCR
    # ========================================================
    
    difference_plot <-
      LR_comparison %>%
      arrange(
        desc(
          abs(
            difference_pCR_minus_nonpCR
          )
        )
      ) %>%
      slice_head(
        n = 15
      ) %>%
      mutate(
        plot_id =
          paste(
            source,
            target,
            LR,
            sep = " | "
          )
      )
    
    
    difference_plot$plot_id <-
      factor(
        difference_plot$plot_id,
        levels =
          rev(
            unique(
              difference_plot$plot_id
            )
          )
      )
    
    
    p_difference <-
      ggplot(
        difference_plot,
        aes(
          x =
            difference_pCR_minus_nonpCR,
          y = plot_id
        )
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      geom_col() +
      scale_y_discrete(
        labels =
          function(x) {
            difference_plot$LR[
              match(
                x,
                difference_plot$plot_id
              )
            ]
          }
      ) +
      labs(
        title =
          "Difference in modeled IL1B–CD8 communication",
        subtitle =
          "pCR minus non-pCR",
        x =
          "Communication probability difference",
        y =
          "Ligand–receptor pair"
      ) +
      theme_classic() +
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),
        axis.text.y =
          element_text(
            size = 8
          )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_4_IL1B_CD8_pCR_vs_nonpCR.pdf"
      ),
      p_difference,
      width = 10,
      height = 8
    )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_4_IL1B_CD8_pCR_vs_nonpCR.png"
      ),
      p_difference,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  
  # ============================================================
  # SAVE FINAL OBJECTS
  # ============================================================
  
  saveRDS(
    cellchat_pCR,
    file.path(
      final_report_dir,
      "CellChat_pCR_FINAL.rds"
    )
  )
  
  saveRDS(
    cellchat_nonpCR,
    file.path(
      final_report_dir,
      "CellChat_nonpCR_FINAL.rds"
    )
  )
  
  saveRDS(
    LR_results,
    file.path(
      final_report_dir,
      "IL1B_CD8_LR_results_FINAL.rds"
    )
  )
  
  
  # ============================================================
  # SAVE ANALYSIS OBJECT
  # ============================================================
  
  if (exists("analysis_data")) {
    
    saveRDS(
      analysis_data,
      file.path(
        final_report_dir,
        "GSE205506_analysis_FINAL.rds"
      )
    )
  }
  
  
  # ============================================================
  # FINAL INVENTORY
  # ============================================================
  
  cat("\n\n")
  cat("====================================================\n")
  cat("        GSE205506 FINAL REPORT COMPLETED\n")
  cat("====================================================\n\n")
  
  cat("Output folder:\n")
  print(
    normalizePath(
      final_report_dir
    )
  )
  
  cat("\n\nFinal files:\n")
  
  print(
    list.files(
      final_report_dir
    )
  )
  
  cat("\n\nDONE.\n")
  
  shell.exec(normalizePath(final_report_dir))
  
  # ============================================================
  # FINAL IL1B-CD8 RESULTS
  # START DIRECTLY FROM EXISTING CELLCHAT RESULTS
  # ============================================================
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # ------------------------------------------------------------
  # 1. CREATE OUTPUT FOLDER
  # ------------------------------------------------------------
  
  final_report_dir <- file.path(
    target_dir,
    "FINAL_REPORT"
  )
  
  if (!dir.exists(final_report_dir)) {
    dir.create(
      final_report_dir,
      recursive = TRUE
    )
  }
  
  cat("\nSaving to:\n")
  cat(normalizePath(final_report_dir))
  cat("\n")
  
  
  # ------------------------------------------------------------
  # 2. EXTRACT ONLY IL1B_Mono <-> CD8_T
  # ------------------------------------------------------------
  
  pCR_IL1B_CD8 <- comm_pCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (source == "IL1B_Mono" & target == "CD8_T") |
        (source == "CD8_T" & target == "IL1B_Mono")
    ) %>%
    mutate(
      response = "pCR"
    )
  
  
  nonpCR_IL1B_CD8 <- comm_nonpCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (source == "IL1B_Mono" & target == "CD8_T") |
        (source == "CD8_T" & target == "IL1B_Mono")
    ) %>%
    mutate(
      response = "non-pCR"
    )
  
  
  # ------------------------------------------------------------
  # 3. CREATE LR NAME
  # ------------------------------------------------------------
  
  if ("interaction_name" %in% colnames(pCR_IL1B_CD8)) {
    
    pCR_IL1B_CD8$LR <- as.character(
      pCR_IL1B_CD8$interaction_name
    )
    
    nonpCR_IL1B_CD8$LR <- as.character(
      nonpCR_IL1B_CD8$interaction_name
    )
    
  } else {
    
    pCR_IL1B_CD8$LR <- paste(
      pCR_IL1B_CD8$ligand,
      pCR_IL1B_CD8$receptor,
      sep = " - "
    )
    
    nonpCR_IL1B_CD8$LR <- paste(
      nonpCR_IL1B_CD8$ligand,
      nonpCR_IL1B_CD8$receptor,
      sep = " - "
    )
  }
  
  
  # ------------------------------------------------------------
  # 4. COMBINE pCR + non-pCR
  # ------------------------------------------------------------
  
  LR_results <- bind_rows(
    pCR_IL1B_CD8,
    nonpCR_IL1B_CD8
  )
  
  
  # ------------------------------------------------------------
  # 5. PRINT THE ACTUAL RESULT
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("IL1B+ MONOCYTE <-> CD8+ T CELL INTERACTIONS\n")
  cat("============================================================\n\n")
  
  if (nrow(LR_results) == 0) {
    
    cat(
      "NO IL1B_Mono <-> CD8_T interactions detected.\n"
    )
    
  } else {
    
    print(
      LR_results,
      row.names = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # 6. SAVE COMPLETE TABLE
  # ------------------------------------------------------------
  
  write.csv(
    LR_results,
    file.path(
      final_report_dir,
      "TABLE_IL1B_CD8_ALL_INTERACTIONS.csv"
    ),
    row.names = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 7. CREATE SIMPLE pCR vs non-pCR TABLE
  # ------------------------------------------------------------
  
  if (
    nrow(LR_results) > 0 &&
    "prob" %in% colnames(LR_results)
  ) {
    
    LR_comparison <- LR_results %>%
      select(
        response,
        source,
        target,
        LR,
        prob
      ) %>%
      group_by(
        response,
        source,
        target,
        LR
      ) %>%
      summarise(
        probability = max(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = response,
        values_from = probability,
        values_fill = 0
      )
    
    if (!"pCR" %in% colnames(LR_comparison)) {
      LR_comparison$pCR <- 0
    }
    
    if (!"non-pCR" %in% colnames(LR_comparison)) {
      LR_comparison$`non-pCR` <- 0
    }
    
    LR_comparison <- LR_comparison %>%
      mutate(
        difference =
          pCR - `non-pCR`
      ) %>%
      arrange(
        desc(
          abs(difference)
        )
      )
    
    write.csv(
      LR_comparison,
      file.path(
        final_report_dir,
        "TABLE_IL1B_CD8_pCR_vs_nonpCR.csv"
      ),
      row.names = FALSE
    )
    
    cat("\n")
    cat("============================================================\n")
    cat("pCR vs non-pCR COMPARISON\n")
    cat("============================================================\n\n")
    
    print(
      LR_comparison,
      row.names = FALSE
    )
    
    
    # ----------------------------------------------------------
    # 8. FIGURE: LR INTERACTIONS
    # ----------------------------------------------------------
    
    top_LR <- LR_results %>%
      group_by(LR) %>%
      summarise(
        max_prob = max(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(
        desc(max_prob)
      ) %>%
      slice_head(n = 20) %>%
      pull(LR)
    
    
    plot_data <- LR_results %>%
      filter(
        LR %in% top_LR
      )
    
    
    plot_data$LR <- factor(
      plot_data$LR,
      levels = rev(top_LR)
    )
    
    
    p_LR <- ggplot(
      plot_data,
      aes(
        x = response,
        y = LR,
        size = prob
      )
    ) +
      
      geom_point(
        alpha = 0.85
      ) +
      
      facet_grid(
        source ~ target,
        scales = "free_y",
        space = "free_y"
      ) +
      
      labs(
        title =
          "IL1B+ monocyte–CD8+ T-cell ligand–receptor interactions",
        x = "Clinical response",
        y = "Ligand–receptor pair",
        size = "Communication probability"
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        axis.text.y = element_text(
          size = 8
        ),
        strip.text = element_text(
          face = "bold"
        )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_IL1B_CD8_LR_INTERACTIONS.pdf"
      ),
      p_LR,
      width = 11,
      height = 9
    )
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_IL1B_CD8_LR_INTERACTIONS.png"
      ),
      p_LR,
      width = 11,
      height = 9,
      dpi = 300
    )
    
    
    # ----------------------------------------------------------
    # 9. FIGURE: DIFFERENCE pCR - non-pCR
    # ----------------------------------------------------------
    
    p_difference <- LR_comparison %>%
      slice_head(n = 15)
    
    p_difference$LR <- factor(
      p_difference$LR,
      levels = rev(p_difference$LR)
    )
    
    
    p_diff <- ggplot(
      p_difference,
      aes(
        x = difference,
        y = LR
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      
      geom_col() +
      
      labs(
        title =
          "Difference in CellChat communication probability",
        subtitle =
          "pCR minus non-pCR",
        x =
          "Communication probability difference",
        y =
          "Ligand–receptor pair"
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        axis.text.y = element_text(
          size = 8
        )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_IL1B_CD8_pCR_minus_nonpCR.pdf"
      ),
      p_diff,
      width = 10,
      height = 8
    )
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_IL1B_CD8_pCR_minus_nonpCR.png"
      ),
      p_diff,
      width = 10,
      height = 8,
      dpi = 300
    )
    
  }
  
  
  # ------------------------------------------------------------
  # 10. SAVE THE EXISTING CELLCHAT OBJECTS
  # ------------------------------------------------------------
  
  saveRDS(
    cellchat_pCR,
    file.path(
      final_report_dir,
      "CellChat_pCR_FINAL.rds"
    )
  )
  
  saveRDS(
    cellchat_nonpCR,
    file.path(
      final_report_dir,
      "CellChat_nonpCR_FINAL.rds"
    )
  )
  
  
  # ------------------------------------------------------------
  # 11. SAVE COMPLETE SESSION
  # ------------------------------------------------------------
  
  save.image(
    file.path(
      final_report_dir,
      "GSE205506_FINAL_SESSION.RData"
    )
  )
  
  
  # ------------------------------------------------------------
  # 12. FINAL MESSAGE
  # ------------------------------------------------------------
  
  cat("\n\n")
  cat("============================================================\n")
  cat("                 FINISHED\n")
  cat("============================================================\n")
  
  cat("\nAll results are here:\n")
  cat(normalizePath(final_report_dir))
  
  cat("\n\nFiles:\n")
  
  print(
    list.files(
      final_report_dir
    )
  )
  
  cat("\n============================================================\n")
  cat("END\n")
  cat("============================================================\n")
  
  
  
  
  
  # ============================================================
  # FINAL IL1B-CD8 RESULTS — ROBUST VERSION
  # ============================================================
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # ------------------------------------------------------------
  # 1. OUTPUT DIRECTORY
  # ------------------------------------------------------------
  
  final_report_dir <- file.path(
    target_dir,
    "FINAL_REPORT"
  )
  
  if (!dir.exists(final_report_dir)) {
    dir.create(
      final_report_dir,
      recursive = TRUE
    )
  }
  
  
  # ------------------------------------------------------------
  # 2. EXTRACT IL1B_Mono <-> CD8_T
  # ------------------------------------------------------------
  
  pCR_IL1B_CD8 <- comm_pCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (source == "IL1B_Mono" & target == "CD8_T") |
        (source == "CD8_T" & target == "IL1B_Mono")
    ) %>%
    mutate(
      response = "pCR"
    )
  
  nonpCR_IL1B_CD8 <- comm_nonpCR_target %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    ) %>%
    filter(
      (source == "IL1B_Mono" & target == "CD8_T") |
        (source == "CD8_T" & target == "IL1B_Mono")
    ) %>%
    mutate(
      response = "non-pCR"
    )
  
  
  # ------------------------------------------------------------
  # 3. CREATE LR NAME
  # ------------------------------------------------------------
  
  if ("interaction_name" %in% colnames(pCR_IL1B_CD8)) {
    
    pCR_IL1B_CD8$LR <- as.character(
      pCR_IL1B_CD8$interaction_name
    )
    
    nonpCR_IL1B_CD8$LR <- as.character(
      nonpCR_IL1B_CD8$interaction_name
    )
    
  } else if (
    all(
      c("ligand", "receptor") %in%
      colnames(pCR_IL1B_CD8)
    )
  ) {
    
    pCR_IL1B_CD8$LR <- paste(
      pCR_IL1B_CD8$ligand,
      pCR_IL1B_CD8$receptor,
      sep = " - "
    )
    
    nonpCR_IL1B_CD8$LR <- paste(
      nonpCR_IL1B_CD8$ligand,
      nonpCR_IL1B_CD8$receptor,
      sep = " - "
    )
    
  } else {
    
    pCR_IL1B_CD8$LR <- paste(
      pCR_IL1B_CD8$source,
      pCR_IL1B_CD8$target,
      sep = " -> "
    )
    
    nonpCR_IL1B_CD8$LR <- paste(
      nonpCR_IL1B_CD8$source,
      nonpCR_IL1B_CD8$target,
      sep = " -> "
    )
  }
  
  
  # ------------------------------------------------------------
  # 4. COMBINE
  # ------------------------------------------------------------
  
  LR_results <- bind_rows(
    pCR_IL1B_CD8,
    nonpCR_IL1B_CD8
  )
  
  
  # ------------------------------------------------------------
  # 5. PRINT ACTUAL BIOLOGICAL RESULT
  # ------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat("IL1B+ MONOCYTE <-> CD8+ T CELL INTERACTIONS\n")
  cat("============================================================\n\n")
  
  if (nrow(LR_results) == 0) {
    
    cat(
      "NO IL1B_Mono <-> CD8_T interactions detected.\n"
    )
    
  } else {
    
    print(
      LR_results,
      row.names = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # 6. SAVE COMPLETE TABLE
  # ------------------------------------------------------------
  
  write.csv(
    LR_results,
    file.path(
      final_report_dir,
      "TABLE_1_IL1B_CD8_ALL_INTERACTIONS.csv"
    ),
    row.names = FALSE
  )
  
  
  # ============================================================
  # 7. CREATE COMPARISON TABLE
  # ============================================================
  
  if (
    nrow(LR_results) > 0 &&
    "prob" %in% colnames(LR_results)
  ) {
    
    LR_comparison <- LR_results %>%
      group_by(
        response,
        source,
        target,
        LR
      ) %>%
      summarise(
        probability = max(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = response,
        values_from = probability,
        values_fill = 0
      )
    
    
    if (!"pCR" %in% colnames(LR_comparison)) {
      LR_comparison$pCR <- 0
    }
    
    if (!"non-pCR" %in% colnames(LR_comparison)) {
      LR_comparison$`non-pCR` <- 0
    }
    
    
    LR_comparison <- LR_comparison %>%
      mutate(
        difference_pCR_minus_nonpCR =
          pCR - `non-pCR`
      ) %>%
      arrange(
        desc(
          abs(
            difference_pCR_minus_nonpCR
          )
        )
      )
    
    
    write.csv(
      LR_comparison,
      file.path(
        final_report_dir,
        "TABLE_2_IL1B_CD8_pCR_vs_nonpCR.csv"
      ),
      row.names = FALSE
    )
    
    
    cat("\n")
    cat("============================================================\n")
    cat("pCR vs non-pCR COMPARISON\n")
    cat("============================================================\n\n")
    
    print(
      LR_comparison,
      row.names = FALSE
    )
    
    
    # ==========================================================
    # 8. PREPARE FIGURE DATA
    # ==========================================================
    
    # Use unique ID for each direction + LR
    # This prevents duplicated factor levels.
    
    plot_data <- LR_results %>%
      mutate(
        interaction_id = paste(
          source,
          target,
          LR,
          sep = " | "
        )
      )
    
    
    # Select strongest interactions
    
    top_ids <- plot_data %>%
      group_by(interaction_id) %>%
      summarise(
        max_prob = max(
          prob,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(
        desc(max_prob)
      ) %>%
      slice_head(
        n = 20
      ) %>%
      pull(interaction_id)
    
    
    plot_data <- plot_data %>%
      filter(
        interaction_id %in% top_ids
      )
    
    
    # Create UNIQUE labels for plot
    
    label_table <- plot_data %>%
      distinct(
        interaction_id,
        LR
      )
    
    
    plot_data <- plot_data %>%
      left_join(
        label_table,
        by = "interaction_id",
        suffix = c("", "_label")
      )
    
    
    # Use interaction_id as factor — guaranteed unique
    
    plot_data$interaction_id <- factor(
      plot_data$interaction_id,
      levels = rev(
        unique(
          plot_data$interaction_id
        )
      )
    )
    
    
    # ==========================================================
    # 9. FIGURE 3 — LR INTERACTIONS
    # ==========================================================
    
    p_LR <- ggplot(
      plot_data,
      aes(
        x = response,
        y = interaction_id,
        size = prob
      )
    ) +
      
      geom_point(
        alpha = 0.85
      ) +
      
      facet_grid(
        source ~ target,
        scales = "free_y",
        space = "free_y"
      ) +
      
      scale_y_discrete(
        labels = function(x) {
          label_table$LR[
            match(
              x,
              label_table$interaction_id
            )
          ]
        }
      ) +
      
      labs(
        title =
          "IL1B+ monocyte–CD8+ T-cell ligand–receptor interactions",
        x = "Clinical response",
        y = "Ligand–receptor pair",
        size = "Communication probability"
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        axis.text.y = element_text(
          size = 8
        ),
        strip.text = element_text(
          face = "bold"
        )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_1_IL1B_CD8_LR_INTERACTIONS.pdf"
      ),
      p_LR,
      width = 11,
      height = 9
    )
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_1_IL1B_CD8_LR_INTERACTIONS.png"
      ),
      p_LR,
      width = 11,
      height = 9,
      dpi = 300
    )
    
    
    # ==========================================================
    # 10. FIGURE 2 — DIFFERENCE pCR vs non-pCR
    # ==========================================================
    
    difference_plot <- LR_comparison %>%
      arrange(
        desc(
          abs(
            difference_pCR_minus_nonpCR
          )
        )
      ) %>%
      slice_head(
        n = 15
      ) %>%
      mutate(
        plot_id = paste(
          source,
          target,
          LR,
          sep = " | "
        )
      )
    
    
    difference_plot$plot_id <- factor(
      difference_plot$plot_id,
      levels = rev(
        unique(
          difference_plot$plot_id
        )
      )
    )
    
    
    p_difference <- ggplot(
      difference_plot,
      aes(
        x = difference_pCR_minus_nonpCR,
        y = plot_id
      )
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      
      geom_col() +
      
      scale_y_discrete(
        labels = function(x) {
          difference_plot$LR[
            match(
              x,
              difference_plot$plot_id
            )
          ]
        }
      ) +
      
      labs(
        title =
          "Difference in modeled IL1B–CD8 communication",
        subtitle =
          "pCR minus non-pCR",
        x =
          "Communication probability difference",
        y =
          "Ligand–receptor pair"
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 13
        ),
        axis.text.y = element_text(
          size = 8
        )
      )
    
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_2_IL1B_CD8_pCR_vs_nonpCR.pdf"
      ),
      p_difference,
      width = 10,
      height = 8
    )
    
    ggsave(
      file.path(
        final_report_dir,
        "FIGURE_2_IL1B_CD8_pCR_vs_nonpCR.png"
      ),
      p_difference,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  
  # ============================================================
  # 11. SAVE CELLCHAT OBJECTS
  # ============================================================
  
  saveRDS(
    cellchat_pCR,
    file.path(
      final_report_dir,
      "CellChat_pCR_FINAL.rds"
    )
  )
  
  saveRDS(
    cellchat_nonpCR,
    file.path(
      final_report_dir,
      "CellChat_nonpCR_FINAL.rds"
    )
  )
  
  
  # ============================================================
  # 12. SAVE COMPLETE SESSION
  # ============================================================
  
  save.image(
    file.path(
      final_report_dir,
      "GSE205506_FINAL_SESSION.RData"
    )
  )
  
  
  # ============================================================
  # 13. FINAL CHECK
  # ============================================================
  
  cat("\n\n")
  cat("============================================================\n")
  cat("                 ANALYSIS SAVED\n")
  cat("============================================================\n")
  
  cat("\nFolder:\n")
  cat(normalizePath(final_report_dir))
  
  cat("\n\nFiles:\n")
  
  print(
    list.files(
      final_report_dir
    )
  )
  
  cat("\n============================================================\n")
  cat("IMPORTANT: READ THE TABLE PRINTED ABOVE\n")
  cat("It contains the actual IL1B_Mono <-> CD8_T interactions.\n")
  cat("============================================================\n")
  
  # ============================================================
  # SAVE ONLY THE IMPORTANT FINAL OBJECTS
  # ============================================================
  
  saveRDS(
    cellchat_pCR,
    file.path(
      final_report_dir,
      "CellChat_pCR_FINAL.rds"
    )
  )
  
  saveRDS(
    cellchat_nonpCR,
    file.path(
      final_report_dir,
      "CellChat_nonpCR_FINAL.rds"
    )
  )
  
  saveRDS(
    LR_results,
    file.path(
      final_report_dir,
      "IL1B_CD8_LR_results_FINAL.rds"
    )
  )
  
  if (exists("LR_comparison")) {
    write.csv(
      LR_comparison,
      file.path(
        final_report_dir,
        "TABLE_2_IL1B_CD8_pCR_vs_nonpCR.csv"
      ),
      row.names = FALSE
    )
  }
  
  cat("\n============================================\n")
  cat("FINAL RESULTS SAVED\n")
  cat("============================================\n\n")
  
  cat("Folder:\n")
  print(normalizePath(final_report_dir))
  
  cat("\nFiles:\n")
  print(list.files(final_report_dir))
  
  normalizePath(final_report_dir)
  
  getwd()
  list.files()
  list.files("R")
  
  list.files("results", recursive = TRUE)
  
  
  list.files(
    path = "C:/Users/Admin/Documents",
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  list.files(
    path = "C:/Users/Admin/Documents/results",
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )


