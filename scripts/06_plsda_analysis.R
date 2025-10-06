# =============================================================================
# 06_plsda_analysis.R
# Partial Least Squares Discriminant Analysis for classification
# =============================================================================

library(tidyverse)
library(mixOmics)
library(caret)
library(flextable)
library(officer)

# Load the data
source("scripts/01_data_loading.R")
data <- load_metabolomics_data()

rawdata_clean <- data$rawdata
metabolite_names <- data$metabolite_names

# Load plot settings
load("output/data/plot_settings.RData")

set.seed(123)

# =============================================================================
# PREPARE DATA FOR PLS-DA
# =============================================================================

message("\n=== Preparing data for PLS-DA ===")

# Extract metabolite matrix
X <- as.matrix(rawdata_clean[, metabolite_names])

# Convert to numeric
X <- apply(X, 2, as.numeric)

# Remove columns with near-zero variance - more robust handling
nzv <- caret::nearZeroVar(X, saveMetrics = FALSE)
if (length(nzv) > 0) {
  message("  - Removing ", length(nzv), " near-zero variance metabolites")
  X <- X[, -nzv, drop = FALSE]
  message("  - Kept ", ncol(X), " metabolites")
} else {
  message("  - No near-zero variance metabolites found")
  message("  - Kept all ", ncol(X), " metabolites")
}

# Impute any remaining NAs with median
na_count <- sum(is.na(X))
if (na_count > 0) {
  message("  - Imputing ", na_count, " missing values")
  for (j in 1:ncol(X)) {
    if (anyNA(X[, j])) {
      X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
    }
  }
}

# Autoscale (mean-center and unit-variance)
X <- scale(X, center = TRUE, scale = TRUE)

# Group vector
y <- factor(rawdata_clean$Group, levels = c("Control", "IHD", "T2D"))

message("  - Final data dimensions: ", nrow(X), " samples × ", ncol(X), " metabolites")
message("  - Groups: Control=", sum(y=="Control"), ", IHD=", sum(y=="IHD"), ", T2D=", sum(y=="T2D"))

# =============================================================================
# PERFORM PLS-DA
# =============================================================================

message("\n=== Running PLS-DA ===")

# Perform PLS-DA with 2 components
ncomp <- 2
plsda_fit <- plsda(X, y, ncomp = ncomp)

# Extract scores and variance explained
scores <- plsda_fit$variates$X
prop_var <- plsda_fit$prop_expl_var$X

message("  - Component 1 explains ", round(prop_var[1] * 100, 1), "% variance")
message("  - Component 2 explains ", round(prop_var[2] * 100, 1), "% variance")

# =============================================================================
# MODEL EVALUATION
# =============================================================================

message("\n=== Evaluating PLS-DA model ===")

# Test different distance metrics
distance_methods <- c("max.dist", "centroids.dist", "mahalanobis.dist")
pred_out <- predict(plsda_fit, newdata = X, ncomp = ncomp)

# Store results for each distance method
summary_list <- list()
predictions_list <- list()  # Store predictions for plotting

for (dist_method in distance_methods) {
  # Extract predictions - THIS IS THE KEY PART
  pred_matrix <- pred_out$class[[dist_method]]
  
  # For component 2
  if (is.matrix(pred_matrix)) {
    pred_cls <- pred_matrix[, ncomp]
  } else {
    pred_cls <- pred_matrix
  }
  pred_cls <- factor(pred_cls, levels = levels(y))
  
  # Store for later plotting
  predictions_list[[dist_method]] <- pred_cls
  
  # Debug: print first 10 predictions to see if they differ
  cat("\n", dist_method, " - First 10 predictions:\n")
  print(pred_cls[1:10])
  
  # Confusion matrix
  cm <- confusionMatrix(pred_cls, y)
  
  # Calculate per-class metrics
  classes <- levels(y)
  per_class_df <- data.frame()
  
  for (cl in classes) {
    TP <- sum(pred_cls == cl & y == cl)
    FN <- sum(pred_cls != cl & y == cl)
    FP <- sum(pred_cls == cl & y != cl)
    TN <- sum(pred_cls != cl & y != cl)
    
    Sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    Specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    Balanced_Accuracy <- mean(c(Sensitivity, Specificity), na.rm = TRUE)
    
    per_class_df <- rbind(
      per_class_df,
      data.frame(
        Distance_Method = dist_method,
        Class = cl,
        Sensitivity = round(Sensitivity, 3),
        Specificity = round(Specificity, 3),
        Balanced_Accuracy = round(Balanced_Accuracy, 3),
        Accuracy = NA_real_,
        Kappa = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Overall metrics row
  overall_row <- data.frame(
    Distance_Method = dist_method,
    Class = "Overall",
    Sensitivity = NA_real_,
    Specificity = NA_real_,
    Balanced_Accuracy = NA_real_,
    Accuracy = round(as.numeric(cm$overall["Accuracy"]), 3),
    Kappa = round(as.numeric(cm$overall["Kappa"]), 3),
    stringsAsFactors = FALSE
  )
  
  summary_list[[dist_method]] <- rbind(per_class_df, overall_row)
  
  cat("\n", dist_method, " Accuracy:", round(cm$overall["Accuracy"], 3), "\n")
}

# Combine all results
summary_table <- bind_rows(summary_list)

# =============================================================================
# CREATE PLS-DA SCORE PLOTS
# =============================================================================

message("\n=== Creating PLS-DA plots ===")

# Create output directory for PLS-DA plots
if (!dir.exists("output/figures/plsda")) {
  dir.create("output/figures/plsda", recursive = TRUE)
}

for (dist_method in distance_methods) {
  # Use stored predictions
  pred_cls <- predictions_list[[dist_method]]
  
  # Create plot dataframe
  plot_df <- data.frame(
    Comp1 = scores[, 1],
    Comp2 = scores[, 2],
    TrueGroup = y,
    PredGroup = pred_cls
  )
  
  # Calculate convex hulls
  hulls <- plot_df %>%
    group_by(TrueGroup) %>%
    slice(chull(Comp1, Comp2)) %>%
    ungroup()
  
  # Create plot
  p <- ggplot(plot_df, aes(x = Comp1, y = Comp2)) +
    geom_polygon(
      data = hulls,
      aes(fill = TrueGroup, group = TrueGroup),
      alpha = 0.2, color = NA
    ) +
    geom_point(aes(color = TrueGroup, shape = PredGroup), size = 3, stroke = 1) +
    scale_color_manual(values = group_colors, name = "True Group") +
    scale_fill_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes, name = "Predicted") +
    labs(
      title = paste("PLS-DA:", dist_method),
      x = paste0("Component 1 (", round(prop_var[1] * 100, 1), "%)"),
      y = paste0("Component 2 (", round(prop_var[2] * 100, 1), "%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    ) +
    guides(fill = "none")
  
  # Save plot
  filename <- paste0("output/figures/plsda/plsda_", gsub("\\.", "_", dist_method), ".png")
  ggsave(filename, p, width = 8, height = 6, dpi = 600, bg = "white")
  message("Saved: ", filename)
}
# =============================================================================
# EXPORT RESULTS
# =============================================================================

message("\n=== Exporting results ===")

# Save summary table as CSV
write.csv(summary_table, "output/tables/plsda_performance_metrics.csv", row.names = FALSE)
message("Saved: output/tables/plsda_performance_metrics.csv")

# Create formatted flextable
ft <- flextable(summary_table) %>%
  set_header_labels(
    Distance_Method = "Distance Method",
    Class = "Class",
    Sensitivity = "Sensitivity",
    Specificity = "Specificity",
    Balanced_Accuracy = "Balanced Accuracy",
    Accuracy = "Overall Accuracy",
    Kappa = "Kappa"
  ) %>%
  colformat_num(j = c("Sensitivity", "Specificity", "Balanced_Accuracy", 
                      "Accuracy", "Kappa"), digits = 3) %>%
  theme_booktabs() %>%
  autofit() %>%
  set_caption("PLS-DA Classification Performance by Distance Method")

# Export to Word
doc <- read_docx() %>%
  body_add_par("PLS-DA Classification Results", style = "heading 1") %>%
  body_add_par("Supervised classification using Partial Least Squares Discriminant Analysis with different distance metrics.", 
               style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "output/tables/plsda_results.docx")
message("Saved: output/tables/plsda_results.docx")

# =============================================================================
# VARIABLE IMPORTANCE
# =============================================================================

message("\n=== Extracting variable importance ===")

# Get VIP scores (Variable Importance in Projection)
vip_scores <- vip(plsda_fit)
vip_df <- data.frame(
  Metabolite = rownames(vip_scores),
  VIP_Comp1 = vip_scores[, 1],
  VIP_Comp2 = vip_scores[, 2]
) %>%
  arrange(desc(VIP_Comp1))

# Save top 20 most important metabolites
write.csv(head(vip_df, 20), "output/tables/plsda_top_metabolites.csv", row.names = FALSE)
message("Saved: output/tables/plsda_top_metabolites.csv")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== PLS-DA analysis complete ===")
cat("\nModel Summary:\n")
cat("  - Number of components:", ncomp, "\n")
cat("  - Variance explained: Comp1 =", round(prop_var[1] * 100, 1), 
    "%, Comp2 =", round(prop_var[2] * 100, 1), "%\n")
cat("  - Best performing distance method:\n")

best_method <- summary_table %>%
  filter(Class == "Overall") %>%
  arrange(desc(Accuracy)) %>%
  slice(1)

cat("    ", best_method$Distance_Method, 
    " (Accuracy =", best_method$Accuracy, 
    ", Kappa =", best_method$Kappa, ")\n")