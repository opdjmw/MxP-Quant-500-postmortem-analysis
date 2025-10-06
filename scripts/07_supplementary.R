# =============================================================================
# 07_supplementary.R
# Generate supplementary tables for publication
# =============================================================================

library(tidyverse)
library(readxl)
library(flextable)
library(officer)

# Load the data
source("scripts/01_data_loading.R")
data <- load_metabolomics_data()

rawdata_clean <- data$rawdata
long_data <- data$long_data
metabolite_info <- data$metabolite_info
metabolite_names <- data$metabolite_names

# =============================================================================
# SUPPLEMENTARY TABLE 1: METABOLITE CONCENTRATION SUMMARY
# =============================================================================

message("\n=== Creating Supplementary Table 1: Metabolite Concentrations ===")

# Back-transform from log2 to original µM scale
supplementary_table_1 <- long_data %>%
  mutate(Concentration_uM = 2^Value) %>%
  group_by(Class, Metabolite) %>%
  summarise(
    Min_uM = round(min(Concentration_uM, na.rm = TRUE), 3),
    Median_uM = round(median(Concentration_uM, na.rm = TRUE), 3),
    Max_uM = round(max(Concentration_uM, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(Class, Metabolite)

# Create flextable
ft_s1 <- flextable(supplementary_table_1) %>%
  set_header_labels(
    Class = "Metabolite Class",
    Metabolite = "Metabolite",
    Min_uM = "Min (µM)",
    Median_uM = "Median (µM)",
    Max_uM = "Max (µM)"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("Supplementary Table 1. Concentration ranges for all detected metabolites")

# Export to Word
doc_s1 <- read_docx() %>%
  body_add_par("Supplementary Table 1", style = "heading 1") %>%
  body_add_par("Concentration Ranges for All Detected Metabolites", style = "heading 2") %>%
  body_add_par("Values represent minimum, median, and maximum concentrations (µM) across all samples.", 
               style = "Normal") %>%
  body_add_flextable(ft_s1)

print(doc_s1, target = "output/tables/supplementary_table_1_concentrations.docx")
write.csv(supplementary_table_1, "output/tables/supplementary_table_1_concentrations.csv", row.names = FALSE)
message("Saved: output/tables/supplementary_table_1_concentrations.docx")
message("Saved: output/tables/supplementary_table_1_concentrations.csv")

# =============================================================================
# SUPPLEMENTARY TABLE 2: METABOLITE CLASS SUMMARY
# =============================================================================

message("\n=== Creating Supplementary Table 2: Metabolite Class Summary ===")

# Summary by metabolite class
supplementary_table_2 <- long_data %>%
  mutate(Concentration_uM = 2^Value) %>%
  group_by(Class) %>%
  summarise(
    N_Metabolites = n_distinct(Metabolite),
    Example_Metabolites = paste(head(unique(Metabolite), 3), collapse = ", "),
    Mean_Concentration = round(mean(Concentration_uM, na.rm = TRUE), 2),
    SD_Concentration = round(sd(Concentration_uM, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_SD = paste0(Mean_Concentration, " (", SD_Concentration, ")")
  ) %>%
  select(Class, N_Metabolites, Example_Metabolites, Mean_SD) %>%
  arrange(desc(N_Metabolites))

# Create flextable
ft_s2 <- flextable(supplementary_table_2) %>%
  set_header_labels(
    Class = "Metabolite Class",
    N_Metabolites = "Number of Metabolites",
    Example_Metabolites = "Example Metabolites",
    Mean_SD = "Mean Concentration (SD) µM"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("Supplementary Table 2. Summary of metabolite classes")

# Export to Word
doc_s2 <- read_docx() %>%
  body_add_par("Supplementary Table 2", style = "heading 1") %>%
  body_add_par("Summary of Metabolite Classes", style = "heading 2") %>%
  body_add_par("Number of metabolites, examples, and mean concentration (± SD) for each metabolite class.", 
               style = "Normal") %>%
  body_add_flextable(ft_s2)

print(doc_s2, target = "output/tables/supplementary_table_2_class_summary.docx")
write.csv(supplementary_table_2, "output/tables/supplementary_table_2_class_summary.csv", row.names = FALSE)
message("Saved: output/tables/supplementary_table_2_class_summary.docx")
message("Saved: output/tables/supplementary_table_2_class_summary.csv")

# =============================================================================
# SUPPLEMENTARY TABLE 3: LOD/LLOQ/ULOQ STATISTICS
# =============================================================================

message("\n=== Creating Supplementary Table 3: LOD/LLOQ/ULOQ Statistics ===")

# Load full dataset to check for LOD markers
if (file.exists("data/Dataset2_µM.xlsx")) {
  dataset2 <- read_excel("data/Dataset2_µM.xlsx", skip = 1)[, 1:650]
  
  # Get text values from samples (rows 32-71)
  dataset2_samples <- dataset2[32:71, ]
  
  # Count LOD/LLOQ/ULOQ occurrences for retained metabolites
  lod_like <- c("< LOD", "< threshold")
  lloq_like <- c("< LLOQ")
  uloq_like <- c("> ULOQ")
  
  lod_summary <- data.frame()
  
  for (met in metabolite_names) {
    if (met %in% colnames(dataset2_samples)) {
      values <- dataset2_samples[[met]]
      
      lod_count <- sum(values %in% lod_like, na.rm = TRUE)
      lloq_count <- sum(values %in% lloq_like, na.rm = TRUE)
      uloq_count <- sum(values %in% uloq_like, na.rm = TRUE)
      
      lod_summary <- rbind(
        lod_summary,
        data.frame(
          Metabolite = met,
          Below_LOD = lod_count,
          Below_LLOQ = lloq_count,
          Above_ULOQ = uloq_count
        )
      )
    }
  }
  
  # Add class information
  lod_summary <- lod_summary %>%
    left_join(metabolite_info, by = "Metabolite") %>%
    select(Class, Metabolite, Below_LOD, Below_LLOQ, Above_ULOQ) %>%
    arrange(Class, Metabolite)
  
  # Create flextable
  ft_s3 <- flextable(lod_summary) %>%
    set_header_labels(
      Class = "Metabolite Class",
      Metabolite = "Metabolite",
      Below_LOD = "Values < LOD",
      Below_LLOQ = "Values < LLOQ",
      Above_ULOQ = "Values > ULOQ"
    ) %>%
    theme_booktabs() %>%
    autofit() %>%
    set_caption("Supplementary Table 3. Detection limit statistics for metabolites")
  
  # Export to Word
  doc_s3 <- read_docx() %>%
    body_add_par("Supplementary Table 3", style = "heading 1") %>%
    body_add_par("Detection Limit Statistics", style = "heading 2") %>%
    body_add_par("Number of samples with values below limit of detection (LOD), below lower limit of quantification (LLOQ), or above upper limit of quantification (ULOQ).", 
                 style = "Normal") %>%
    body_add_flextable(ft_s3)
  
  print(doc_s3, target = "output/tables/supplementary_table_3_detection_limits.docx")
  write.csv(lod_summary, "output/tables/supplementary_table_3_detection_limits.csv", row.names = FALSE)
  message("Saved: output/tables/supplementary_table_3_detection_limits.docx")
  message("Saved: output/tables/supplementary_table_3_detection_limits.csv")
  
} else {
  message("Note: Dataset2_µM.xlsx not found. Skipping LOD/LLOQ/ULOQ table.")
}

# =============================================================================
# SUPPLEMENTARY TABLE 4: GROUP CHARACTERISTICS
# =============================================================================

message("\n=== Creating Supplementary Table 4: Group Characteristics ===")

# Calculate group-level summary statistics
group_summary <- long_data %>%
  mutate(Concentration_uM = 2^Value) %>%
  group_by(Group) %>%
  summarise(
    N_Samples = n_distinct(`Sample identification`),
    N_Metabolites = n_distinct(Metabolite),
    Mean_Concentration = round(mean(Concentration_uM, na.rm = TRUE), 2),
    SD_Concentration = round(sd(Concentration_uM, na.rm = TRUE), 2),
    Median_Concentration = round(median(Concentration_uM, na.rm = TRUE), 2),
    .groups = "drop"
  )

# Create flextable
ft_s4 <- flextable(group_summary) %>%
  set_header_labels(
    Group = "Group",
    N_Samples = "Number of Samples",
    N_Metabolites = "Metabolites Detected",
    Mean_Concentration = "Mean Conc. (µM)",
    SD_Concentration = "SD (µM)",
    Median_Concentration = "Median Conc. (µM)"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("Supplementary Table 4. Group characteristics")

# Export to Word
doc_s4 <- read_docx() %>%
  body_add_par("Supplementary Table 4", style = "heading 1") %>%
  body_add_par("Group Characteristics", style = "heading 2") %>%
  body_add_par("Summary statistics for each experimental group.", style = "Normal") %>%
  body_add_flextable(ft_s4)

print(doc_s4, target = "output/tables/supplementary_table_4_group_characteristics.docx")
write.csv(group_summary, "output/tables/supplementary_table_4_group_characteristics.csv", row.names = FALSE)
message("Saved: output/tables/supplementary_table_4_group_characteristics.docx")
message("Saved: output/tables/supplementary_table_4_group_characteristics.csv")

# =============================================================================
# COMBINED SUPPLEMENTARY DOCUMENT
# =============================================================================

message("\n=== Creating combined supplementary document ===")

# Create one Word document with all supplementary tables
doc_combined <- read_docx()

# Add each table with proper headings
doc_combined <- doc_combined %>%
  body_add_par("Supplementary Materials", style = "Title") %>%
  body_add_par("Cardiac Metabolomics Analysis", style = "heading 1") %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_par("Supplementary Table 1: Metabolite Concentrations", style = "heading 2") %>%
  body_add_par("Concentration ranges (minimum, median, maximum in µM) for all detected metabolites.", 
               style = "Normal") %>%
  body_add_flextable(ft_s1) %>%
  body_add_break() %>%
  
  body_add_par("Supplementary Table 2: Metabolite Class Summary", style = "heading 2") %>%
  body_add_par("Number of metabolites per class with mean concentration and standard deviation.", 
               style = "Normal") %>%
  body_add_flextable(ft_s2) %>%
  body_add_break()

# Add LOD table if it exists
if (exists("ft_s3")) {
  doc_combined <- doc_combined %>%
    body_add_par("Supplementary Table 3: Detection Limit Statistics", style = "heading 2") %>%
    body_add_par("Number of measurements below/above detection limits for each metabolite.", 
                 style = "Normal") %>%
    body_add_flextable(ft_s3) %>%
    body_add_break()
}

doc_combined <- doc_combined %>%
  body_add_par("Supplementary Table 4: Group Characteristics", style = "heading 2") %>%
  body_add_par("Summary statistics for each experimental group.", style = "Normal") %>%
  body_add_flextable(ft_s4)

print(doc_combined, target = "output/tables/supplementary_materials_combined.docx")
message("Saved: output/tables/supplementary_materials_combined.docx")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Supplementary tables complete ===")
cat("\nGenerated supplementary materials:\n")
cat("  - Table S1: Metabolite concentrations (", nrow(supplementary_table_1), "metabolites)\n")
cat("  - Table S2: Metabolite class summary (", nrow(supplementary_table_2), "classes)\n")
if (exists("lod_summary")) {
  cat("  - Table S3: Detection limit statistics (", nrow(lod_summary), "metabolites)\n")
}
cat("  - Table S4: Group characteristics (", nrow(group_summary), "groups)\n")
cat("  - Combined supplementary document\n")
cat("\nAll files saved to: output/tables/\n")