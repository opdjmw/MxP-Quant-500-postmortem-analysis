# =============================================================================
# 04_volcano_plots.R
# Generate volcano plots for pairwise group comparisons
# =============================================================================

library(tidyverse)
library(ggrepel)
library(flextable)
library(officer)

# Load the data
source("scripts/01_data_loading.R")
data <- load_metabolomics_data()

long_data <- data$long_data

# Load plot settings
load("output/data/plot_settings.RData")

# =============================================================================
# VOLCANO PLOT PARAMETERS
# =============================================================================

# Significance thresholds
p_threshold <- 0.05
fc_threshold <- 0.58  # log2(1.5)

# Color palette
volcano_colors <- c(
  "Upregulated" = "#D7191C",      # red
  "Downregulated" = "#2C7BB6",    # blue
  "Not Significant" = "grey70"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Calculate volcano plot data for two groups
#'
#' @param data Long format data
#' @param group1 First group name
#' @param group2 Second group name (reference)
#' @return Dataframe with log2FC and p-values
compute_volcano_data <- function(data, group1, group2) {
  subset_data <- data %>% filter(Group %in% c(group1, group2))
  
  results <- subset_data %>%
    group_by(Metabolite) %>%
    summarise(
      mean_g1 = mean(Value[Group == group1], na.rm = TRUE),
      mean_g2 = mean(Value[Group == group2], na.rm = TRUE),
      p_value = wilcox.test(
        Value[Group == group1], 
        Value[Group == group2], 
        exact = FALSE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      log2FC = mean_g1 - mean_g2,
      neg_log10_pval = -log10(p_value)
    ) %>%
    select(Metabolite, log2FC, p_value, neg_log10_pval)
  
  return(results)
}

# =============================================================================
# CALCULATE VOLCANO DATA FOR ALL COMPARISONS
# =============================================================================

message("\n=== Calculating volcano plot data ===")

volcano_ihd_ctrl <- compute_volcano_data(long_data, "IHD", "Control") %>%
  mutate(Comparison = "IHD vs Control")

volcano_t2d_ctrl <- compute_volcano_data(long_data, "T2D", "Control") %>%
  mutate(Comparison = "T2D vs Control")

volcano_t2d_ihd <- compute_volcano_data(long_data, "T2D", "IHD") %>%
  mutate(Comparison = "T2D vs IHD")

# Combine all comparisons
volcano_combined <- bind_rows(
  volcano_ihd_ctrl, 
  volcano_t2d_ctrl, 
  volcano_t2d_ihd
) %>%
  mutate(
    significance = case_when(
      p_value < p_threshold & log2FC > fc_threshold ~ "Upregulated",
      p_value < p_threshold & log2FC < -fc_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  group_by(Comparison) %>%
  mutate(rank_p = rank(p_value)) %>%
  ungroup()

# Identify top metabolites for labeling (top 10 by p-value per comparison)
top_labeled <- volcano_combined %>%
  group_by(Comparison) %>%
  filter(rank_p <= 10) %>%
  ungroup()

# Get all significant metabolites for labeling
sig_labeled <- volcano_combined %>%
  filter(significance != "Not Significant")

message("Calculated volcano data for 3 comparisons")

# =============================================================================
# CREATE FACETED VOLCANO PLOT
# =============================================================================

message("\n=== Creating volcano plots ===")

volcano_plot <- ggplot(
  volcano_combined, 
  aes(x = log2FC, y = neg_log10_pval, color = significance)
) +
  geom_point(size = 2.2, alpha = 0.7) +
  scale_color_manual(values = volcano_colors) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = sig_labeled,
    aes(label = Metabolite),
    size = 3,
    color = "black",
    max.overlaps = 15,
    box.padding = 0.5
  ) +
  facet_wrap(~Comparison, nrow = 1, scales = "free_x") +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("output/figures/volcano_plot_comparisons.png", 
       volcano_plot, width = 14, height = 6, dpi = 600, bg = "white")
message("Saved: output/figures/volcano_plot_comparisons.png")

# =============================================================================
# EXPORT SUPPLEMENTARY TABLE
# =============================================================================

message("\n=== Creating supplementary table ===")

# Create supplementary table with significant metabolites
si_table <- volcano_combined %>%
  filter(significance != "Not Significant") %>%
  select(Metabolite, log2FC, p_value, Comparison, significance) %>%
  arrange(Comparison, p_value) %>%
  mutate(
    log2FC = round(log2FC, 2),
    p_value = signif(p_value, 3)
  ) %>%
  rename(
    `Metabolite` = Metabolite,
    `Log2 Fold Change` = log2FC,
    `P-value` = p_value,
    `Comparison` = Comparison,
    `Direction` = significance
  )

# Save CSV
write.csv(si_table, "output/tables/volcano_plot_significant.csv", row.names = FALSE)
message("Saved: output/tables/volcano_plot_significant.csv")

# Create formatted flextable
ft <- flextable(si_table) %>%
  autofit() %>%
  color(i = ~ Direction == "Upregulated", j = "Direction", color = "#D7191C") %>%
  color(i = ~ Direction == "Downregulated", j = "Direction", color = "#2C7BB6") %>%
  bold(i = ~ `P-value` < 0.01, bold = TRUE) %>%
  theme_booktabs() %>%
  set_caption("Significantly altered metabolites in pairwise comparisons")

# Export to Word
doc <- read_docx() %>%
  body_add_par("Supplementary Table: Volcano Plot Results", style = "heading 1") %>%
  body_add_par("Metabolites with significant differences (p < 0.05, |log2FC| > 0.58) between groups.", 
               style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "output/tables/volcano_plot_results.docx")
message("Saved: output/tables/volcano_plot_results.docx")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

summary_volcano <- volcano_combined %>%
  group_by(Comparison, significance) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = significance, values_from = Count, values_fill = 0)

cat("\nVolcano Plot Summary:\n")
print(summary_volcano)

write.csv(summary_volcano, "output/tables/volcano_summary.csv", row.names = FALSE)

message("\n=== Volcano plot analysis complete ===")
message("All outputs saved to output/ directory")