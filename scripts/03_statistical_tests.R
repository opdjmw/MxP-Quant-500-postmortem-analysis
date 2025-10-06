# =============================================================================
# 03_statistical_tests.R
# Statistical analysis: Kruskal-Wallis, Dunn's test, PERMANOVA
# =============================================================================

library(tidyverse)
library(FSA)
library(vegan)
library(flextable)
library(officer)

# Load the data
source("scripts/01_data_loading.R")
data <- load_metabolomics_data()

rawdata_clean <- data$rawdata
long_data <- data$long_data
metabolite_names <- data$metabolite_names

# Load plot settings
load("output/data/plot_settings.RData")

# =============================================================================
# 1. KRUSKAL-WALLIS TESTS FOR ALL METABOLITES
# =============================================================================

message("\n=== Running Kruskal-Wallis tests ===")

# Perform Kruskal-Wallis test for each metabolite
kruskal_results <- long_data %>%
  group_by(Metabolite) %>%
  summarise(
    H_stat = kruskal.test(Value ~ Group)$statistic,
    p_value = kruskal.test(Value ~ Group)$p.value,
    n = n(),
    k = n_distinct(Group),
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate epsilon squared effect size
    Epsilon_sq = (H_stat - k + 1) / (n - k),
    # FDR correction
    p_adjusted = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(p_value)

# Print summary
cat("\nKruskal-Wallis Test Summary:\n")
cat("  Total metabolites tested:", nrow(kruskal_results), "\n")
cat("  Significant (raw p < 0.05):", sum(kruskal_results$p_value < 0.05), "\n")
cat("  Significant (FDR < 0.05):", sum(kruskal_results$p_adjusted < 0.05), "\n")

# Save full results
write.csv(kruskal_results, "output/tables/kruskal_wallis_results.csv", row.names = FALSE)
message("Saved: output/tables/kruskal_wallis_results.csv")

# Create formatted table for top 20 metabolites
top_20 <- kruskal_results %>%
  head(20) %>%
  select(Metabolite, p_value, p_adjusted, Epsilon_sq) %>%
  mutate(
    p_value = sprintf("%.4f", p_value),
    p_adjusted = sprintf("%.4f", p_adjusted),
    Epsilon_sq = sprintf("%.3f", Epsilon_sq)
  )

ft_kruskal <- flextable(top_20) %>%
  set_caption("Top 20 Metabolites - Kruskal-Wallis Test Results") %>%
  set_header_labels(
    Metabolite = "Metabolite",
    p_value = "P-value",
    p_adjusted = "FDR-adjusted P",
    Epsilon_sq = "Effect Size (ε²)"
  ) %>%
  theme_vanilla() %>%
  autofit()

# Export to Word
doc <- read_docx() %>%
  body_add_par("Kruskal-Wallis Test Results", style = "heading 1") %>%
  body_add_flextable(ft_kruskal)

print(doc, target = "output/tables/kruskal_results_top20.docx")
message("Saved: output/tables/kruskal_results_top20.docx")

# =============================================================================
# 2. DUNN'S POST-HOC TEST FOR SIGNIFICANT METABOLITES
# =============================================================================

message("\n=== Running Dunn's post-hoc tests ===")

# Filter metabolites with significant raw p-value
significant_metabolites <- kruskal_results %>%
  filter(p_value < 0.05) %>%
  pull(Metabolite)

cat("  Running Dunn's test for", length(significant_metabolites), "significant metabolites\n")

# Run Dunn's test for each significant metabolite
dunn_results_list <- list()

for (met in significant_metabolites) {
  met_data <- long_data %>% filter(Metabolite == met)
  
  # Perform Dunn's test with Bonferroni correction
  dunn_test <- dunnTest(Value ~ Group, data = met_data, method = "bonferroni")
  
  # Extract results and add metabolite name
  dunn_df <- dunn_test$res
  dunn_df$Metabolite <- met
  dunn_results_list[[met]] <- dunn_df
}

# Combine all results
dunn_results_combined <- bind_rows(dunn_results_list) %>%
  mutate(across(c(Z, P.unadj, P.adj), ~ round(., 3))) %>%
  select(Metabolite, Comparison, Z, P.unadj, P.adj)

# Save results
write.csv(dunn_results_combined, "output/tables/dunn_posthoc_results.csv", row.names = FALSE)
message("Saved: output/tables/dunn_posthoc_results.csv")

# Create formatted table
ft_dunn <- flextable(dunn_results_combined) %>%
  set_caption("Dunn's Post-Hoc Test Results for Significant Metabolites") %>%
  set_header_labels(
    Metabolite = "Metabolite",
    Comparison = "Comparison",
    Z = "Z-statistic",
    P.unadj = "P-value",
    P.adj = "Adjusted P"
  ) %>%
  theme_vanilla() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Dunn's Post-Hoc Test Results", style = "heading 1") %>%
  body_add_flextable(ft_dunn)

print(doc, target = "output/tables/dunn_posthoc_results.docx")
message("Saved: output/tables/dunn_posthoc_results.docx")

# =============================================================================
# 3. BOXPLOTS FOR SIGNIFICANT METABOLITES
# =============================================================================

message("\n=== Creating boxplots for significant metabolites ===")

# Boxplot for all significant metabolites
boxplot_sig <- ggplot(
  long_data %>% filter(Metabolite %in% significant_metabolites),
  aes(x = Group, y = Value, fill = Group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  facet_wrap(~ Metabolite, scales = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Log2 Concentration (µM)")

ggsave("output/figures/boxplot_significant_metabolites.png", 
       boxplot_sig, width = 12, height = 8, dpi = 600, bg = "white")
message("Saved: output/figures/boxplot_significant_metabolites.png")

# =============================================================================
# 4. PERMANOVA
# =============================================================================

message("\n=== Running PERMANOVA ===")

# Extract metabolite matrix (samples × features)
metabolite_matrix <- as.matrix(rawdata_clean[, metabolite_names])

# Ensure all numeric
metabolite_matrix <- apply(metabolite_matrix, 2, as.numeric)

# Create grouping vector
group_vector <- rawdata_clean$Group

# Calculate distance matrix using Bray-Curtis
distance_matrix <- vegdist(metabolite_matrix, method = "bray")

# Run PERMANOVA
set.seed(123)
permanova_result <- adonis2(distance_matrix ~ group_vector, permutations = 999)

cat("\nPERMANOVA Results:\n")
print(permanova_result)

# Save PERMANOVA results
permanova_df <- as.data.frame(permanova_result)
permanova_df$Term <- rownames(permanova_df)
permanova_df <- permanova_df %>%
  select(Term, everything())

write.csv(permanova_df, "output/tables/permanova_results.csv", row.names = FALSE)
message("Saved: output/tables/permanova_results.csv")

# Create formatted table
ft_permanova <- flextable(permanova_df) %>%
  set_caption("PERMANOVA Results - Group Effect on Metabolite Profile") %>%
  theme_vanilla() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("PERMANOVA Results", style = "heading 1") %>%
  body_add_par("Multivariate analysis testing for differences in metabolite profiles between groups.", 
               style = "Normal") %>%
  body_add_flextable(ft_permanova)

print(doc, target = "output/tables/permanova_results.docx")
message("Saved: output/tables/permanova_results.docx")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Statistical Testing Complete ===")
message("Outputs saved:")
message("  - Kruskal-Wallis results (CSV and Word)")
message("  - Dunn's post-hoc results (CSV and Word)")
message("  - Boxplots of significant metabolites (PNG)")
message("  - PERMANOVA results (CSV and Word)")

# Create summary report
summary_stats <- data.frame(
  Test = c("Kruskal-Wallis", "Dunn's Post-Hoc", "PERMANOVA"),
  Description = c(
    paste(nrow(kruskal_results), "metabolites tested;", 
          sum(kruskal_results$p_adjusted < 0.05), "significant (FDR < 0.05)"),
    paste(nrow(dunn_results_combined), "pairwise comparisons for", 
          length(significant_metabolites), "metabolites"),
    paste("Overall group effect: p =", 
          round(permanova_result$`Pr(>F)`[1], 4))
  )
)

write.csv(summary_stats, "output/tables/statistical_tests_summary.csv", row.names = FALSE)
message("\nSummary saved to: output/tables/statistical_tests_summary.csv")