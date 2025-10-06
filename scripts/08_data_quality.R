# =============================================================================
# 08_data_quality.R
# Data quality assessment: LOD/LLOQ/ULOQ statistics and retention analysis
#
# REQUIRES BOTH DATASETS:
#   Dataset1_µM.xlsx (463 metabolites, log2-transformed, cleaned)
#     - Used for concentration statistics of retained metabolites
#   
#   Dataset2_µM.xlsx (628 metabolites, original µM scale with TEXT FLAGS)
#     - WebIDQ configured to replace threshold values with text flags
#     - Flags: "< LOD", "< threshold", "< LLOQ", "> ULOQ"
#     - Used for detection limit statistics
# =============================================================================

library(tidyverse)
library(readxl)
library(flextable)
library(officer)

# Avoid MASS::select conflict
select <- dplyr::select

# Check if Dataset2 exists
if (!file.exists("data/Dataset2_µM.xlsx")) {
  message("\n=== WARNING ===")
  message("Dataset2_µM.xlsx not found in data/ directory")
  message("This script requires Dataset2 for LOD/LLOQ/ULOQ statistics")
  message("Skipping data quality analysis...")
  message("===============\n")
  quit(save = "no", status = 0)
}

message("\n=== Starting Data Quality Analysis ===")

# Load Dataset1 (retained metabolites)
source("scripts/01_data_loading.R")
data1 <- load_metabolomics_data("data/Dataset1_µM.xlsx")

# =============================================================================
# LOAD DATASET2 (FULL DATASET WITH TEXT FLAGS)
# =============================================================================

message("\n=== Loading Dataset2 (full dataset with quality flags) ===")

# Load Dataset1 first to get structure
dataset1_raw <- read_excel("data/Dataset1_µM.xlsx", skip = 1)[, 1:485]

# Load Dataset2 - keep as much as possible
dataset2_raw <- read_excel("data/Dataset2_µM.xlsx", skip = 1)[, 1:650]

# Extract all metabolite info from Dataset2
all_metabolite_names <- colnames(dataset2_raw[, 23:650])
all_metabolite_classes <- as.character(unlist(dataset2_raw[1, 23:650]))
all_metabolite_classes[is.na(all_metabolite_classes) | all_metabolite_classes == ""] <- "Other Information"

all_metabolite_info <- data.frame(
  Metabolite = all_metabolite_names,
  Class = all_metabolite_classes,
  stringsAsFactors = FALSE
)

message("  - Dataset2 has ", length(all_metabolite_names), " total metabolites")

# =============================================================================
# EXTRACT LOD/LLOQ/ULOQ STATISTICS FROM DATASET2 TEXT FLAGS
# =============================================================================

message("\n=== Extracting LOD/LLOQ/ULOQ statistics from text flags ===")

# Get sample rows AS TEXT (rows 32-71, do NOT convert to numeric)
dataset2_samples_text <- dataset2_raw[32:71, 23:650]
n_samples <- nrow(dataset2_samples_text)

# Define threshold text patterns (from WebIDQ flags)
lod_patterns <- c("< LOD", "< threshold")
lloq_patterns <- c("< LLOQ")
uloq_patterns <- c("> ULOQ")

# Count flags for all metabolites
lod_lloq_uloq <- data.frame()

message("  - Counting quality flags across ", length(all_metabolite_names), " metabolites...")

for (met in all_metabolite_names) {
  if (met %in% colnames(dataset2_samples_text)) {
    # Keep as character to preserve text flags
    values <- as.character(dataset2_samples_text[[met]])
    
    # Count each flag type
    lod_count <- sum(values %in% lod_patterns, na.rm = TRUE)
    lloq_count <- sum(values %in% lloq_patterns, na.rm = TRUE)
    uloq_count <- sum(values %in% uloq_patterns, na.rm = TRUE)
    
    lod_lloq_uloq <- rbind(
      lod_lloq_uloq,
      data.frame(
        Metabolite = met,
        Values_Below_LOD = lod_count,
        Values_Below_LLOQ = lloq_count,
        Values_Above_ULOQ = uloq_count,
        Percent_LOD = round(100 * lod_count / n_samples, 1),
        Percent_LLOQ = round(100 * lloq_count / n_samples, 1),
        Percent_ULOQ = round(100 * uloq_count / n_samples, 1),
        stringsAsFactors = FALSE
      )
    )
  }
}

# Add class information
lod_lloq_uloq <- lod_lloq_uloq %>%
  left_join(all_metabolite_info, by = "Metabolite")

message("  - Found quality flags in ", 
        sum(lod_lloq_uloq$Values_Below_LOD > 0), " metabolites (< LOD)")
message("  - Found quality flags in ", 
        sum(lod_lloq_uloq$Values_Below_LLOQ > 0), " metabolites (< LLOQ)")
message("  - Found quality flags in ", 
        sum(lod_lloq_uloq$Values_Above_ULOQ > 0), " metabolites (> ULOQ)")

# =============================================================================
# COMPARE RETAINED VS EXCLUDED METABOLITES
# =============================================================================

message("\n=== Comparing retained vs excluded metabolites ===")

# Identify retained and excluded metabolites
retained_metabolites <- data1$metabolite_names
excluded_metabolites <- setdiff(all_metabolite_names, retained_metabolites)

cat("  - Retained:", length(retained_metabolites), "\n")
cat("  - Excluded:", length(excluded_metabolites), "\n")

# Count by class
retention_by_class <- all_metabolite_info %>%
  mutate(Status = ifelse(Metabolite %in% retained_metabolites, "Retained", "Excluded")) %>%
  group_by(Class, Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Status, values_from = Count, values_fill = 0) %>%
  mutate(
    Total = Retained + Excluded,
    Retention_Percent = round(100 * Retained / Total, 1)
  ) %>%
  arrange(desc(Total))

# =============================================================================
# GET CONCENTRATION STATISTICS FOR RETAINED METABOLITES
# =============================================================================

message("\n=== Calculating concentration statistics for retained metabolites ===")

# Back-transform from log2 to µM
concentration_stats <- data1$long_data %>%
  mutate(Concentration_uM = 2^Value) %>%
  group_by(Metabolite, Class) %>%
  summarise(
    Mean_uM = mean(Concentration_uM, na.rm = TRUE),
    SD_uM = sd(Concentration_uM, na.rm = TRUE),
    Median_uM = median(Concentration_uM, na.rm = TRUE),
    Min_uM = min(Concentration_uM, na.rm = TRUE),
    Max_uM = max(Concentration_uM, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# CREATE COMPREHENSIVE SUMMARY TABLE
# =============================================================================

message("\n=== Creating comprehensive summary table ===")

# Filter to classes with ≥8 total metabolites
classes_to_include <- retention_by_class %>%
  filter(Total >= 8) %>%
  pull(Class)

# Aggregate by class
summary_table <- retention_by_class %>%
  filter(Class %in% classes_to_include) %>%
  select(Class, Total, Retained, Retention_Percent)

# Add concentration stats (only for retained metabolites in these classes)
conc_by_class <- concentration_stats %>%
  filter(Class %in% classes_to_include) %>%
  group_by(Class) %>%
  summarise(
    Mean_Conc = round(mean(Mean_uM, na.rm = TRUE), 1),
    SD_Conc = round(sd(Mean_uM, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_SD = paste0(Mean_Conc, " (", SD_Conc, ")")
  )

# Add LOD/LLOQ/ULOQ stats (for all metabolites in these classes)
lod_by_class <- lod_lloq_uloq %>%
  filter(Class %in% classes_to_include) %>%
  group_by(Class) %>%
  summarise(
    Mean_LOD = round(mean(Percent_LOD, na.rm = TRUE), 1),
    Mean_LLOQ = round(mean(Percent_LLOQ, na.rm = TRUE), 1),
    Mean_ULOQ = round(mean(Percent_ULOQ, na.rm = TRUE), 1),
    .groups = "drop"
  )

# Combine everything
summary_table <- summary_table %>%
  left_join(conc_by_class[, c("Class", "Mean_SD")], by = "Class") %>%
  left_join(lod_by_class, by = "Class") %>%
  arrange(desc(Retained))

# Add example metabolites
examples_by_class <- concentration_stats %>%
  filter(Class %in% classes_to_include) %>%
  group_by(Class) %>%
  slice_head(n = 2) %>%
  summarise(Examples = paste(Metabolite, collapse = ", "), .groups = "drop")

summary_table <- summary_table %>%
  left_join(examples_by_class, by = "Class")

# Reorder columns
summary_table <- summary_table %>%
  select(
    Class,
    Total,
    Retained,
    `Retention (%)` = Retention_Percent,
    Examples,
    `Mean Conc (SD) µM` = Mean_SD,
    `% < LOD` = Mean_LOD,
    `% < LLOQ` = Mean_LLOQ,
    `% > ULOQ` = Mean_ULOQ
  )

# =============================================================================
# EXPORT RESULTS
# =============================================================================

message("\n=== Exporting results ===")

# Save comprehensive summary
write.csv(summary_table, "output/tables/data_quality_summary.csv", row.names = FALSE)
message("Saved: output/tables/data_quality_summary.csv")

# Create formatted table
ft_summary <- flextable(summary_table) %>%
  set_caption("Data Quality Summary: Metabolite Retention and Detection Limits") %>%
  set_header_labels(
    Class = "Metabolite Class",
    Total = "Total Metabolites",
    Retained = "Retained",
    `Retention (%)` = "Retention (%)",
    Examples = "Example Metabolites",
    `Mean Conc (SD) µM` = "Mean Conc (SD) µM",
    `% < LOD` = "% < LOD",
    `% < LLOQ` = "% < LLOQ",
    `% > ULOQ` = "% > ULOQ"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(j = 2:9, align = "right", part = "all")

# Export to Word
doc <- read_docx() %>%
  body_add_par("Data Quality Assessment", style = "heading 1") %>%
  body_add_par("Summary of metabolite retention, concentration statistics, and detection limit information for classes with ≥8 total metabolites. Detection limit percentages based on WebIDQ text flags in Dataset2.", 
               style = "Normal") %>%
  body_add_flextable(ft_summary)

print(doc, target = "output/tables/data_quality_summary.docx")
message("Saved: output/tables/data_quality_summary.docx")

# Save detailed per-metabolite statistics (for retained metabolites only)
detailed_stats <- concentration_stats %>%
  left_join(lod_lloq_uloq %>% 
              select(Metabolite, Values_Below_LOD, Values_Below_LLOQ, Values_Above_ULOQ), 
            by = "Metabolite") %>%
  arrange(Class, Metabolite) %>%
  mutate(across(c(Mean_uM, SD_uM, Median_uM, Min_uM, Max_uM), ~ round(., 3)))

write.csv(detailed_stats, "output/tables/metabolite_detailed_statistics.csv", row.names = FALSE)
message("Saved: output/tables/metabolite_detailed_statistics.csv")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Data quality analysis complete ===")
cat("\nSummary:\n")
cat("  - Total metabolites in Dataset2:", length(all_metabolite_names), "\n")
cat("  - Retained metabolites:", length(retained_metabolites), "\n")
cat("  - Excluded metabolites:", length(excluded_metabolites), "\n")
cat("  - Classes with ≥8 metabolites:", nrow(summary_table), "\n")
cat("\nTables created:\n")
cat("  - data_quality_summary (CSV + Word)\n")
cat("  - metabolite_detailed_statistics (CSV)\n")
