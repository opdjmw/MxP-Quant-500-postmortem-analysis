# =============================================================================
# 01_data_loading.R
# Central data loading and preprocessing function for metabolomics analysis
#
# DATASET INFORMATION:
# This script loads Dataset1_µM.xlsx, which contains:
#   - 463 retained metabolites (≥75% sample presence)
#   - Log2-transformed concentrations: log2(µM)
#   - Pre-imputed values: 46 < LOD, 18 > ULOQ
#   - Target normalization: Median QC Level 2
#
# Dataset2_µM.xlsx (optional, used only by script 08):
#   - 628 total metabolites (original scale, µM)
#   - Contains quality flags: "< LOD", "< LLOQ", "> ULOQ"
# =============================================================================

library(readxl)
library(tidyverse)

#' Load and preprocess metabolomics data
#'
#' @param file_path Path to the Excel file (Dataset1_µM.xlsx or Dataset2_µM.xlsx)
#' @param impute_method Method for imputing missing values ("LOD" or "min")
#' @return List containing cleaned data, long format data, and metadata
#' @export
load_metabolomics_data <- function(file_path = "data/Dataset1_µM.xlsx", 
                                   impute_method = "LOD") {
  
  message("Loading data from: ", file_path)
  
  # ---------------------- Step 1: Read Raw Data ----------------------
  rawdata <- read_excel(file_path, skip = 1)[, 1:485]
  
  # Extract metabolite names and class labels
  metabolite_names <- colnames(rawdata[, 23:485])
  metabolite_classes <- as.character(unlist(rawdata[1, 23:485]))
  metabolite_classes[is.na(metabolite_classes) | metabolite_classes == ""] <- "Other Information"
  
  # Create class info table
  metabolite_info <- data.frame(
    Metabolite = metabolite_names,
    Class = metabolite_classes,
    stringsAsFactors = FALSE
  )
  
  message("  - Found ", length(metabolite_names), " metabolites")
  message("  - Found ", length(unique(metabolite_classes)), " metabolite classes")
  
  # ---------------------- Step 2: Extract LOD Rows ----------------------
  # Row 22 = MS LOD, Row 25 = FIA LOD – use whichever is not -Inf
  lod_row_ms  <- suppressWarnings(as.numeric(unlist(rawdata[22, 23:485])))
  lod_row_fia <- suppressWarnings(as.numeric(unlist(rawdata[25, 23:485])))
  
  # Use MS LOD if available, else FIA LOD
  lod_combined <- ifelse(is.finite(lod_row_ms), lod_row_ms, lod_row_fia)
  names(lod_combined) <- metabolite_names
  
  message("  - Extracted LOD values")
  
  # ---------------------- Step 3: Prepare Data ----------------------
  # Remove non-data rows (e.g., labels, LODs, etc.)
  rawdata_clean <- rawdata[-(1:30), ]
  
  # Ensure numeric conversion for metabolite matrix
  rawdata_clean[, 23:485] <- lapply(rawdata_clean[, 23:485], as.numeric)
  
  message("  - Removed metadata rows (1-30)")
  message("  - Sample size: ", nrow(rawdata_clean), " samples")
  
  # ---------------------- Step 4: Impute Missing Values ----------------------
  metab_matrix <- as.data.frame(rawdata_clean[, 23:485])
  
  n_missing_before <- sum(is.na(metab_matrix))
  
  for (metab in metabolite_names) {
    vec <- metab_matrix[[metab]]
    
    if (impute_method == "LOD") {
      # Use LOD/2 if finite, else min(vec)/2
      lod_val <- lod_combined[[metab]]
      if (is.finite(lod_val)) {
        impute_val <- lod_val / 2
      } else {
        min_val <- min(vec, na.rm = TRUE)
        impute_val <- min_val / 2
      }
    } else if (impute_method == "min") {
      # Use minimum/2
      min_val <- min(vec, na.rm = TRUE)
      impute_val <- min_val / 2
    }
    
    # Replace NAs with the imputed value
    vec[is.na(vec)] <- impute_val
    metab_matrix[[metab]] <- vec
  }
  
  n_missing_after <- sum(is.na(metab_matrix))
  message("  - Imputed ", n_missing_before - n_missing_after, " missing values using ", impute_method, " method")
  
  # Replace the cleaned matrix in main data
  rawdata_clean[, 23:485] <- metab_matrix
  
  # ---------------------- Step 5: Long Format for Analysis ----------------------
  long_data <- rawdata_clean %>%
    pivot_longer(cols = 23:485, names_to = "Metabolite", values_to = "Value") %>%
    left_join(metabolite_info, by = "Metabolite")
  
  # Standardize group names if needed: DM2 -> T2D
  if ("Group" %in% colnames(long_data)) {
    long_data$Group <- recode(long_data$Group, DM2 = "T2D", Dm2 = "T2D", dm2 = "T2D")
    rawdata_clean$Group <- recode(rawdata_clean$Group, DM2 = "T2D", Dm2 = "T2D", dm2 = "T2D")
    message("  - Standardized group labels (DM2 -> T2D)")
  }
  
  message("  - Created long format dataset")
  
  # ---------------------- Return Results ----------------------
  result <- list(
    rawdata = rawdata_clean,
    long_data = long_data,
    metabolite_info = metabolite_info,
    metabolite_names = metabolite_names,
    lod_values = lod_combined
  )
  
  message("=== Data loading complete ===\n")
  
  return(result)

}
