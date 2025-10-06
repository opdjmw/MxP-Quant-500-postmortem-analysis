# =============================================================================
# 05_pca_analysis.R
# Principal Component Analysis of metabolomics data
# =============================================================================

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(flextable)
library(officer)

# Load the data
source("scripts/01_data_loading.R")
data <- load_metabolomics_data()

rawdata_clean <- data$rawdata
metabolite_names <- data$metabolite_names

# Load plot settings
load("output/data/plot_settings.RData")

# =============================================================================
# PERFORM PCA
# =============================================================================

message("\n=== Performing PCA ===")

# Extract metabolite matrix
metabolite_data <- as.data.frame(rawdata_clean[, metabolite_names])

# Perform PCA with scaling
pca_result <- PCA(metabolite_data, scale.unit = TRUE, ncp = 10, graph = FALSE)

message("  - PCA complete with 10 components")

# Extract PCA scores
pca_scores <- as.data.frame(pca_result$ind$coord)
pca_scores$Group <- rawdata_clean$Group
pca_scores$Sample_ID <- rownames(rawdata_clean)

# Extract variance explained
explained_variance <- pca_result$eig[, 2]

cat("\nVariance explained by first 5 PCs:\n")
for (i in 1:5) {
  cat(sprintf("  PC%d: %.1f%%\n", i, explained_variance[i]))
}

# =============================================================================
# PCA PLOTTING FUNCTION
# =============================================================================

#' Create PCA score plot for two components
#'
#' @param pc_x PC number for x-axis
#' @param pc_y PC number for y-axis
#' @return ggplot object
plot_pca_pair <- function(pc_x, pc_y) {
  var_x <- round(explained_variance[pc_x], 1)
  var_y <- round(explained_variance[pc_y], 1)
  
  pc_x_col <- paste0("Dim.", pc_x)
  pc_y_col <- paste0("Dim.", pc_y)
  
  # Calculate convex hulls for each group
  hulls <- pca_scores %>%
    group_by(Group) %>%
    slice(chull(.data[[pc_x_col]], .data[[pc_y_col]])) %>%
    ungroup()
  
  p <- ggplot(pca_scores, aes_string(x = pc_x_col, y = pc_y_col, 
                                     color = "Group", shape = "Group")) +
    geom_polygon(
      data = hulls,
      aes_string(x = pc_x_col, y = pc_y_col, fill = "Group", group = "Group"),
      alpha = 0.2, color = NA, show.legend = FALSE
    ) +
    geom_point(size = 3, stroke = 1, alpha = 1) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    scale_fill_manual(values = group_colors) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste0("PCA: PC", pc_x, " vs PC", pc_y),
      x = paste0("PC", pc_x, " (", var_x, "%)"),
      y = paste0("PC", pc_y, " (", var_y, "%)")
    )
  
  return(p)
}

# =============================================================================
# GENERATE PCA PLOTS
# =============================================================================

message("\n=== Creating PCA plots ===")

# PC1 vs PC2
p_12 <- plot_pca_pair(1, 2)
ggsave("output/figures/pca_pc1_pc2.png", p_12, 
       width = 8, height = 6, dpi = 600, bg = "white")
message("Saved: output/figures/pca_pc1_pc2.png")

# PC1 vs PC3
p_13 <- plot_pca_pair(1, 3)
ggsave("output/figures/pca_pc1_pc3.png", p_13, 
       width = 8, height = 6, dpi = 600, bg = "white")
message("Saved: output/figures/pca_pc1_pc3.png")

# PC2 vs PC3
p_23 <- plot_pca_pair(2, 3)
ggsave("output/figures/pca_pc2_pc3.png", p_23, 
       width = 8, height = 6, dpi = 600, bg = "white")
message("Saved: output/figures/pca_pc2_pc3.png")

# =============================================================================
# SUPPLEMENTARY TABLES
# =============================================================================

message("\n=== Creating supplementary tables ===")

# Table 1: Eigenvalues and variance explained (PC1-10)
eig_df <- as.data.frame(pca_result$eig)
colnames(eig_df) <- c("Eigenvalue", "Variance_Percent", "Cumulative_Percent")
eig_df$PC <- paste0("PC", 1:nrow(eig_df))
eig_df <- eig_df[1:10, c("PC", "Eigenvalue", "Variance_Percent", "Cumulative_Percent")]
eig_df[, 2:4] <- round(eig_df[, 2:4], 2)

ft_eigenvalues <- flextable(eig_df) %>%
  set_header_labels(
    PC = "Principal Component",
    Eigenvalue = "Eigenvalue",
    Variance_Percent = "Variance (%)",
    Cumulative_Percent = "Cumulative (%)"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("PCA Eigenvalues and Explained Variance (PC1-10)")

# Table 2: Top contributing metabolite per PC
contrib_df <- as.data.frame(pca_result$var$contrib)
top_contrib <- data.frame(
  PC = character(),
  Top_Metabolite = character(),
  Contribution_Percent = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:10) {
  colname <- paste0("Dim.", i)
  if (colname %in% colnames(contrib_df)) {
    max_idx <- which.max(contrib_df[[colname]])
    top_contrib <- rbind(
      top_contrib,
      data.frame(
        PC = paste0("PC", i),
        Top_Metabolite = rownames(contrib_df)[max_idx],
        Contribution_Percent = round(contrib_df[max_idx, colname], 2)
      )
    )
  }
}

ft_contrib <- flextable(top_contrib) %>%
  set_header_labels(
    PC = "Principal Component",
    Top_Metabolite = "Top Contributing Metabolite",
    Contribution_Percent = "Contribution (%)"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("Top Contributing Metabolite per Principal Component (PC1-10)")

# Export both tables to Word
doc <- read_docx() %>%
  body_add_par("Supplementary Table: PCA Results", style = "heading 1") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table S1. PCA Eigenvalues and Explained Variance", style = "heading 2") %>%
  body_add_flextable(ft_eigenvalues) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table S2. Top Contributing Metabolites", style = "heading 2") %>%
  body_add_flextable(ft_contrib)

print(doc, target = "output/tables/pca_supplementary_tables.docx")
message("Saved: output/tables/pca_supplementary_tables.docx")

# Save PCA scores
write.csv(pca_scores, "output/tables/pca_scores.csv", row.names = FALSE)
message("Saved: output/tables/pca_scores.csv")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== PCA analysis complete ===")
cat("\nSummary:\n")
cat("  - Total variance explained by PC1-3:", 
    round(sum(explained_variance[1:3]), 1), "%\n")
cat("  - Plots saved: PC1-2, PC1-3, PC2-3\n")
cat("  - Supplementary tables created\n")