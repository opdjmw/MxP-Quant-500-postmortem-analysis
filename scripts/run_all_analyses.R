# =============================================================================
# run_all_analyses.R
# Master script to run complete metabolomics analysis pipeline
# =============================================================================

# Clear workspace
rm(list = ls())

# Record start time
start_time <- Sys.time()

cat("\n")
cat("============================================================\n")
cat("  CARDIAC METABOLOMICS ANALYSIS PIPELINE\n")
cat("============================================================\n\n")

# =============================================================================
# STEP 0: SETUP AND PACKAGE INSTALLATION
# =============================================================================

cat("[Step 0/7] Running setup and package installation...\n")
source("scripts/00_setup.R")
cat("\n✓ Setup complete\n\n")

# =============================================================================
# STEP 1: DATA LOADING
# =============================================================================

cat("[Step 1/7] Loading and preprocessing data...\n")
source("scripts/01_data_loading.R")

# Load data
data <- load_metabolomics_data(file_path = "data/Dataset1_µM.xlsx")

# Extract objects for use in subsequent scripts
rawdata_clean <- data$rawdata
long_data <- data$long_data
metabolite_info <- data$metabolite_info
metabolite_names <- data$metabolite_names
lod_combined <- data$lod_values

cat("\n✓ Data loading complete\n\n")

# =============================================================================
# STEP 2: DATA INSPECTION
# =============================================================================

cat("[Step 2/7] Running data inspection and quality control...\n")
tryCatch({
  source("scripts/02_inspection.R")
  cat("\n✓ Data inspection complete\n\n")
}, error = function(e) {
  cat("\n✗ Error in data inspection:", e$message, "\n\n")
})

# =============================================================================
# STEP 3: STATISTICAL TESTS
# =============================================================================

cat("[Step 3/7] Running statistical tests (Kruskal-Wallis, PERMANOVA)...\n")
tryCatch({
  source("scripts/03_statistical_tests.R")
  cat("\n✓ Statistical tests complete\n\n")
}, error = function(e) {
  cat("\n✗ Error in statistical tests:", e$message, "\n\n")
})

# =============================================================================
# STEP 4: VOLCANO PLOTS
# =============================================================================

cat("[Step 4/7] Generating volcano plots...\n")
tryCatch({
  source("scripts/04_volcano_plots.R")
  cat("\n✓ Volcano plots complete\n\n")
}, error = function(e) {
  cat("\n✗ Error in volcano plots:", e$message, "\n\n")
})

# =============================================================================
# STEP 5: PCA ANALYSIS
# =============================================================================

cat("[Step 5/7] Running Principal Component Analysis...\n")
tryCatch({
  source("scripts/05_pca_analysis.R")
  cat("\n✓ PCA analysis complete\n\n")
}, error = function(e) {
  cat("\n✗ Error in PCA:", e$message, "\n\n")
})

# =============================================================================
# STEP 6: PLS-DA CLASSIFICATION
# =============================================================================

cat("[Step 6/7] Running PLS-DA classification...\n")
tryCatch({
  source("scripts/06_plsda_analysis.R")
  cat("\n✓ PLS-DA analysis complete\n\n")
}, error = function(e) {
  cat("\n✗ Error in PLS-DA:", e$message, "\n\n")
})

# =============================================================================
# STEP 7: SUPPLEMENTARY MATERIALS
# =============================================================================

cat("[Step 7/8] Generating supplementary tables...\n")
tryCatch({
  source("scripts/07_supplementary.R")
  cat("\n✓ Supplementary materials complete\n\n")
}, error = function(e) {
  cat("\n✗ Error in supplementary materials:", e$message, "\n\n")
})

# =============================================================================
# STEP 8: DATA QUALITY ASSESSMENT (OPTIONAL - REQUIRES DATASET2)
# =============================================================================

cat("[Step 8/8] Running data quality assessment...\n")
if (file.exists("data/Dataset2_µM.xlsx")) {
  tryCatch({
    source("scripts/08_data_quality.R")
    cat("\n✓ Data quality assessment complete\n\n")
  }, error = function(e) {
    cat("\n✗ Error in data quality assessment:", e$message, "\n\n")
  })
} else {
  cat("  - Dataset2_µM.xlsx not found - skipping (optional)\n\n")
}

# =============================================================================
# COMPLETION AND SESSION INFO
# =============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("============================================================\n")
cat("  ANALYSIS PIPELINE COMPLETE\n")
cat("============================================================\n\n")
cat("Total runtime:", round(duration, 2), "minutes\n\n")
cat("Output locations:\n")
cat("  - Figures: output/figures/\n")
cat("  - Tables: output/tables/\n")
cat("  - Data: output/data/\n\n")

# Save session info
session_file <- "output/session_info.txt"
writeLines(capture.output(sessionInfo()), session_file)
cat("Session info saved to:", session_file, "\n\n")

# List output files
cat("Generated files:\n")

# Figures
figures <- list.files("output/figures", pattern = "\\.(png|pdf)$", recursive = TRUE)
if (length(figures) > 0) {
  cat("\n  Figures (", length(figures), "):\n")
  for (f in head(figures, 10)) {
    cat("    -", f, "\n")
  }
  if (length(figures) > 10) {
    cat("    ... and", length(figures) - 10, "more\n")
  }
}

# Tables
tables <- list.files("output/tables", pattern = "\\.(csv|docx)$", recursive = TRUE)
if (length(tables) > 0) {
  cat("\n  Tables (", length(tables), "):\n")
  for (t in head(tables, 10)) {
    cat("    -", t, "\n")
  }
  if (length(tables) > 10) {
    cat("    ... and", length(tables) - 10, "more\n")
  }
}

cat("\n")
cat("============================================================\n")
cat("  Analysis completed successfully!\n")
cat("============================================================\n\n")
