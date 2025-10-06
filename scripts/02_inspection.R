# =============================================================================
# 02_inspection.R
# Data inspection and quality control for metabolomics dataset
# =============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(viridis)

# Load the data using central function
source("scripts/01_data_loading.R")
data <- load_metabolomics_data()

# Extract data objects
rawdata_clean <- data$rawdata
long_data <- data$long_data
metabolite_info <- data$metabolite_info
metabolite_names <- data$metabolite_names

# Load plot settings
load("output/data/plot_settings.RData")

# =============================================================================
# 1. Q-Q PLOTS FOR NORMALITY ASSESSMENT
# =============================================================================

#' Generate Q-Q plots for a specific metabolite class
#'
#' @param data Long format data
#' @param class_name Name of metabolite class to plot
plot_qq_by_class <- function(data, class_name) {
  data %>%
    filter(Class == class_name) %>%
    ggplot(aes(sample = Value)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    facet_wrap(~Metabolite, scales = "free") +
    labs(title = paste("Q-Q Plots:", class_name)) +
    theme_minimal()
}

# Generate Q-Q plots for all classes
qq_all <- long_data %>%
  ggplot(aes(sample = Value)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~Class, scales = "free") +
  labs(title = "Q-Q Plots for Metabolite Classes") +
  theme_minimal()

ggsave("output/figures/qq_plots_all_classes.png", qq_all, width = 12, height = 10, dpi = 300)
message("Saved: output/figures/qq_plots_all_classes.png")

# Example: Single class
plot_qq_by_class(long_data, "Phosphatidylcholines")

# =============================================================================
# 2. HISTOGRAMS WITH DENSITY CURVES
# =============================================================================

#' Plot histograms for metabolites in a class
#'
#' @param data Long format data
#' @param class_name Metabolite class name
plot_histograms <- function(data, class_name) {
  data %>%
    filter(Class == class_name) %>%
    ggplot(aes(x = Value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                   fill = "lightblue", alpha = 0.7, color = "black") +
    geom_density(color = "red", linewidth = 1) +
    facet_wrap(~Metabolite, scales = "free") +
    labs(title = paste("Histograms:", class_name),
         x = "Log2 Concentration",
         y = "Density") +
    theme_minimal()
}

# Example usage
p_hist <- plot_histograms(long_data, "Acylcarnitines")
ggsave("output/figures/histogram_acylcarnitines.png", p_hist, width = 12, height = 8, dpi = 300)

# =============================================================================
# 3. BOXPLOTS BY GROUP AND CLASS
# =============================================================================

# Overall distribution by class
boxplot_by_class <- ggplot(long_data, aes(x = Class, y = log10(Value), fill = Class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of Metabolite Concentrations by Class",
       x = "Metabolite Class",
       y = expression(Log[10]~Concentration~(µM))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill = "none")

ggsave("output/figures/boxplot_by_class.png", boxplot_by_class, 
       width = 10, height = 6, dpi = 300, bg = "white")
message("Saved: output/figures/boxplot_by_class.png")

# Filter to classes with ≥8 metabolites
filtered_data <- long_data %>%
  group_by(Class) %>%
  filter(n_distinct(Metabolite) >= 8) %>%
  ungroup()

boxplot_filtered <- ggplot(filtered_data, aes(x = Class, y = log10(Value), fill = Class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  theme_minimal() +
  labs(title = "Metabolite Classes with ≥8 Members",
       x = "Metabolite Class",
       y = expression(Log[10]~Concentration~(µM))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill = "none")

ggsave("output/figures/boxplot_by_class_filtered.png", boxplot_filtered, 
       width = 10, height = 6, dpi = 300, bg = "white")

# Boxplot by group
boxplot_by_group <- ggplot(long_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  facet_wrap(~Class, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of Metabolite Classes by Group",
       x = "Group",
       y = "Log2 Concentration (µM)") +
  scale_fill_manual(values = group_colors)

ggsave("output/figures/boxplot_by_group.png", boxplot_by_group, 
       width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: output/figures/boxplot_by_group.png")

# =============================================================================
# 4. SPECIFIC METABOLITE PLOTTING
# =============================================================================

#' Create boxplot for a specific metabolite
#'
#' @param data Long format data
#' @param class_name Metabolite class
#' @param metabolite_name Specific metabolite name
plot_specific_metabolite <- function(data, class_name, metabolite_name) {
  subset_data <- data %>% 
    filter(Class == class_name, Metabolite == metabolite_name)
  
  p <- ggplot(subset_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Distribution of", metabolite_name, "by Group"),
         x = NULL,
         y = "Metabolite Concentration (log2 µM)") +
    scale_fill_manual(values = group_colors) +
    theme(
      panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}

# Example: C3 from Acylcarnitines
p_c3 <- plot_specific_metabolite(long_data, "Acylcarnitines", "C3")
ggsave("output/figures/boxplot_C3.png", p_c3, width = 8, height = 6, dpi = 600, bg = "white")

# =============================================================================
# 5. NORMALITY TESTING
# =============================================================================

message("\n=== Running normality tests ===")

# Shapiro-Wilk test for all metabolites
shapiro_results <- sapply(metabolite_names, function(met) {
  values <- rawdata_clean[[met]]
  if (length(values) >= 3 && length(values) <= 5000) {
    shapiro.test(values)$p.value
  } else {
    NA
  }
})

# Kolmogorov-Smirnov test
ks_results <- sapply(metabolite_names, function(met) {
  values <- rawdata_clean[[met]]
  if (length(values) >= 3) {
    ks.test(values, "pnorm", mean(values), sd(values))$p.value
  } else {
    NA
  }
})

# Create summary dataframe
normality_summary <- data.frame(
  Metabolite = metabolite_names,
  Shapiro_p = shapiro_results,
  KS_p = ks_results,
  Shapiro_normal = shapiro_results > 0.05,
  KS_normal = ks_results > 0.05
) %>%
  left_join(metabolite_info, by = "Metabolite")

# Overall summary
cat("\nShapiro-Wilk Test Summary:\n")
cat("  Normally distributed (p > 0.05):", sum(shapiro_results > 0.05, na.rm = TRUE), "\n")
cat("  Non-normal (p ≤ 0.05):", sum(shapiro_results <= 0.05, na.rm = TRUE), "\n")
cat("  Proportion normal:", round(mean(shapiro_results > 0.05, na.rm = TRUE), 3), "\n\n")

cat("Kolmogorov-Smirnov Test Summary:\n")
cat("  Normally distributed (p > 0.05):", sum(ks_results > 0.05, na.rm = TRUE), "\n")
cat("  Non-normal (p ≤ 0.05):", sum(ks_results <= 0.05, na.rm = TRUE), "\n")
cat("  Proportion normal:", round(mean(ks_results > 0.05, na.rm = TRUE), 3), "\n")

# Save results
write.csv(normality_summary, "output/tables/normality_test_results.csv", row.names = FALSE)
message("\nSaved: output/tables/normality_test_results.csv")

# Summary by class
class_normality <- normality_summary %>%
  group_by(Class) %>%
  summarise(
    Total_Metabolites = n(),
    SW_Normal_Count = sum(Shapiro_normal, na.rm = TRUE),
    KS_Normal_Count = sum(KS_normal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Metabolites))

write.csv(class_normality, "output/tables/normality_by_class.csv", row.names = FALSE)
message("Saved: output/tables/normality_by_class.csv")

# =============================================================================
# 6. DENSITY PLOTS BY GROUP
# =============================================================================

#' Generate density plots for metabolites in a class
#'
#' @param data Long format data
#' @param class_name Metabolite class name
plot_density_by_class <- function(data, class_name) {
  data %>%
    filter(Class == class_name) %>%
    ggplot(aes(x = Value, fill = Group, color = Group)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = paste("Density Plots:", class_name),
         x = "Metabolite Level",
         y = "Density") +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors)
}

# Example
p_density <- plot_density_by_class(long_data, "Acylcarnitines")
ggsave("output/figures/density_acylcarnitines.png", p_density, 
       width = 12, height = 8, dpi = 300)

message("\n=== Inspection complete ===")
message("All figures saved to output/figures/")
message("All tables saved to output/tables/")