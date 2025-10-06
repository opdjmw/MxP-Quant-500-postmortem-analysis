# =============================================================================
# 00_setup.R
# Package installation and project setup for cardiac metabolomics analysis
# =============================================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Define all required packages
required_packages <- c(
  # Data manipulation
  "tidyverse", "dplyr", "tidyr", "purrr",
  
  # Data import
  "readxl",
  
  # Visualization
  "ggplot2", "ggrepel", "gplots", "ggfortify", "pheatmap", 
  "ggpubr", "corrplot", "patchwork", "viridis",
  
  # Statistical analysis
  "FactoMineR", "factoextra", "gtools", "FSA", "psych", 
  "performance", "lme4", "rcompanion", "vegan", "caret",
  
  # Report generation
  "flextable", "officer", "knitr"
)

# Bioconductor packages
bioc_packages <- c("missMDA", "mixOmics", "ropls")

# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# Function to install Bioconductor packages
install_bioc_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing Bioconductor package:", pkg, "..."))
      BiocManager::install(pkg, update = FALSE)
    }
  }
}

# Install all packages
message("Installing CRAN packages...")
install_and_load(required_packages)

message("\nInstalling Bioconductor packages...")
install_bioc_packages(bioc_packages)

# Create output directory structure
output_dirs <- c(
  "output/figures",
  "output/tables",
  "output/data"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(paste("Created directory:", dir))
  }
}

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
  message("Created data/ directory - please place your Excel files here")
}

# Set plotting defaults
theme_set(theme_minimal())

# Color palettes for consistent visualization
group_colors <- c(
  "IHD" = "#82B7FF",
  "T2D" = "#FF9289", 
  "Control" = "#00D65C"
)

group_shapes <- c(
  "IHD" = 17,
  "T2D" = 16,
  "Control" = 15
)

# Save these for use in other scripts
save(group_colors, group_shapes, file = "output/data/plot_settings.RData")

message("\n=== Setup Complete ===")
message("Package installation finished")
message("Output directories created")
message("Ready to run analysis scripts")
message("\nNext step: Run 01_data_loading.R")