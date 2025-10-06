# Cardiac Metabolomics Analysis Pipeline

Analysis of targeted metabolomics data from cardiac tissue samples using the Biocrates MxP Quant 500 kit.

## Project Overview

This project analyzes metabolomic profiles from three groups:
- **IHD**: Ischemic Heart Disease
- **T2D**: Type 2 Diabetes (formerly DM2)
- **Control**: Control samples

## Data Source

This pipeline analyzes data from the **Biocrates MxP Quant 500 kit**, a targeted metabolomics platform. Raw data is processed through the **MetIDQ™ (WebIDQ) software** provided by Biocrates.

### WebIDQ Platform Features

The WebIDQ platform allows users to:
- Apply quality control filters (e.g., 75% sample presence threshold)
- Select normalization methods (Median QC, internal standards, etc.)
- Choose concentration units (µM, log2-transformed, etc.)
- Export data with or without quality flags (< LOD, < LLOQ, > ULOQ)
- Generate multiple dataset versions with different filtering criteria

**For this analysis**, two exports were created from WebIDQ:

1. **Dataset1_µM.xlsx** (Analysis-ready)
   - Filtered: 75% sample presence threshold applied
   - Transformed: log2(µM) scale
   - Cleaned: LOD/ULOQ values pre-imputed
   - Export setting: "Concentration [log2(µM)] | Target Normalization [Median QC Level 2]"
   - Result: 463 metabolites retained

2. **Dataset2_µM.xlsx** (Quality reference)
   - Unfiltered: All detected metabolites included
   - Original scale: µM (not log-transformed)
   - Uncleaned: Contains quality flag text ("< LOD", "< LLOQ", "> ULOQ")
   - Export setting: "Concentration [µM] | Target Normalization [Median QC Level 2]"
   - Result: 628 total metabolites

### Creating Your Own Datasets

To replicate this analysis with your own data:

1. **Process samples** using Biocrates MxP Quant 500 kit
2. **Upload to WebIDQ** and complete quality control
3. **Export Dataset1** (for analysis):
   - Apply 75% presence filter
   - Select log2 transformation
   - Enable imputation for < LOD and > ULOQ values
   - Choose your normalization method
4. **Export Dataset2** (optional, for quality stats):
   - Remove all filters
   - Keep original µM scale
   - Include quality flag annotations

**Note**: The specific filtering and export options available in WebIDQ may vary depending on your platform version and license. Consult the [Biocrates documentation](https://biocrates.com) or your WebIDQ manual for details on available export configurations.
## Repository Structure

```
├── README.md
├── scripts/
│   ├── 00_setup.R              # Package installation and setup
│   ├── 01_data_loading.R       # Data loading and preprocessing
│   ├── 02_inspection.R         # Data inspection and QC
│   ├── 03_statistical_tests.R  # Kruskal-Wallis, PERMANOVA
│   ├── 04_volcano_plots.R      # Volcano plot generation
│   ├── 05_pca_analysis.R       # Principal Component Analysis
│   ├── 06_plsda_analysis.R     # PLS-DA classification
│   └── 07_supplementary.R      # Supplementary tables
├── data/
│   ├── Dataset1_µM.xlsx        # Cleaned metabolomics data (463 metabolites)
│   └── Dataset2_µM.xlsx        # Full dataset (628 metabolites)
└── output/
    ├── figures/
    └── tables/
```

## Getting Started

### Prerequisites

R version ≥ 4.0.0 recommended

### Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/cardiac-metabolomics.git
cd cardiac-metabolomics
```

2. Install required packages:
```r
source("scripts/00_setup.R")
```

3. Set up your data directory structure:
```r
# Create output directories
dir.create("output/figures", recursive = TRUE)
dir.create("output/tables", recursive = TRUE)
```

## Analysis Workflow

Run scripts in numerical order:

### 1. Data Loading (`01_data_loading.R`)
- Loads raw metabolomics data
- Extracts metabolite names and class labels
- Extracts LOD (Limit of Detection) values
- Imputes missing values using LOD/2 method
- Creates long-format dataset for analysis

### 2. Data Inspection (`02_inspection.R`)
- Q-Q plots for normality assessment
- Histograms and density plots
- Boxplots by group and metabolite class
- Shapiro-Wilk and Kolmogorov-Smirnov tests

### 3. Statistical Analysis (`03_statistical_tests.R`)
- Kruskal-Wallis tests for group differences
- Dunn's post-hoc tests
- Effect size calculations (epsilon squared)
- PERMANOVA for multivariate analysis

### 4. Volcano Plots (`04_volcano_plots.R`)
- Pairwise comparisons (IHD vs Control, T2D vs Control, T2D vs IHD)
- Fold-change and p-value visualization
- Identification of significantly altered metabolites

### 5. PCA (`05_pca_analysis.R`)
- Principal Component Analysis
- Score plots with convex hulls
- Variance explained tables

### 6. PLS-DA (`06_plsda_analysis.R`)
- Supervised classification
- Model validation with multiple distance metrics
- Performance metrics (sensitivity, specificity, balanced accuracy)

### 7. Supplementary Materials (`07_supplementary.R`)
- Summary tables for publication
- Metabolite concentration ranges
- LOD/LLOQ/ULOQ statistics

## Key Features

### Data Processing
- **Missing value imputation**: LOD/2 method
- **Data transformation**: Log2 transformation for normalization
- **Quality control**: Retention threshold ≥75% sample presence

### Statistical Methods
- Non-parametric tests (Kruskal-Wallis)
- Multiple testing correction (FDR)
- Effect size estimation
- Multivariate methods (PCA, PLS-DA, PERMANOVA)

### Visualization
- Publication-ready plots (600 DPI)
- Color-blind friendly palettes
- Consistent styling across figures

## Output Files

### Figures
- `boxplot_significant.png`: Boxplots of significant metabolites
- `volcano_facet.png`: Three-panel volcano plots
- `PCA_PC1_PC3.png`: PCA score plot
- `PLSDA_*.png`: PLS-DA classification plots
- `heatmap_PMI.png`: PMI effect heatmap
- `heatmap_CauseOfDeath.png`: Cause of death heatmap

### Tables (Word format)
- `kruskal_results.docx`: Statistical test results
- `Dunn_Test_Results.docx`: Post-hoc comparisons
- `Supplementary_Volcano_Table.docx`: Volcano plot data
- `Supplementary_Metabolite_Table.docx`: Concentration ranges
- `PCA_SupplementaryTables.docx`: PCA statistics
- `PLSDA_DistanceMethods_Summary.docx`: PLS-DA performance

## Important Notes

### File Paths
All file paths are relative to the project root. Modify the `data_dir` variable in `01_data_loading.R` if your data is stored elsewhere.

### Data Requirements
- **Dataset1_µM.xlsx**: 463 retained metabolites (≥75% presence)
- **Dataset2_µM.xlsx**: Full 628 metabolites (pre-filtering)
- First 30 rows contain metadata (class labels, LOD values)
- Sample data starts at row 31

### Group Naming
The code recodes "DM2" to "T2D" for consistency. Original data may use either nomenclature.

## Citation

If you use this code, please cite contact oscar.puntervold@sund.ku.dk for citation details as the original paper related to this work has not been submitted.



## Contact

Oscar Emil Puntervold
oscar.puntervold@sund.ku.dk  
djyoloprolo@gmail,com

## Troubleshooting

### Common Issues

**"File not found" errors**: 
- Ensure data files are in the `data/` directory
- Check that file names match exactly (case-sensitive)

**Package installation fails**:
- Update R to latest version
- Install Bioconductor packages separately if needed

**Memory issues with large datasets**:
- Close other R sessions
- Increase memory limit: `memory.limit(size = 16000)` (Windows)

### Getting Help

- Open an issue on GitHub
- Check script comments for specific function documentation
- Review the original paper for methodological details

## Session Info

To ensure reproducibility, record your session info after running analyses:
```r
sessionInfo()
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
```

## Acknowledgments

Analysis pipeline developed for cardiac metabolomics study using Biocrates MxP Quant 500 kit.


Last updated: October 2025
