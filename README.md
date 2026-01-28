# BgNorm

<!-- badges: start -->
[![R-CMD-check](https://github.com/Malvikakh/BgNorm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Malvikakh/BgNorm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

**BgNorm** is a Gaussian mixture model (GMM)-based method for background correction and normalization of spatial imaging data, particularly multiplexed spatial proteomics images. The package models intensity distributions using a three-component mixture model and performs marker- and pixel- or cell-level signal adjustment to remove technical background while preserving biological variation.  

BgNorm provides tools for processing multi-channel imaging intensities and generating background-adjusted measurements suitable for downstream spatial and single-cell analyses.

### Mathematical Background

BgNorm assumes three types of pixels:
- **Background** ($X_1 = U_1$): Only background signals
- **Non-specific binding** ($X_2 = U_1 + U_2$): Background + non-specific binding
- **Signal** ($X_3 = U_1 + U_2 + U_3$): Background + non-specific + true biological signal

The observed intensities are modeled as a 3-component Gaussian mixture:

$$f(X) = \sum_{j=1}^{3} \pi_j \mathcal{N}(X; \mu_j, \sigma_j^2)$$

For signal pixels (component 3), the background-corrected intensity is calculated using conditional expectation:

$$E(U_3 | X, C=3) = (\mu_3 - \mu_2) + w \cdot (X - \mu_3)$$

where $w = 1 - \rho$ is a variance weight derived from the correlation structure, and $\rho = \text{sign} \cdot \frac{\min(\sigma_3^2, \sigma_2^2)}{\sigma_3^2}$.

The final adjusted intensity accounts for uncertainty in component membership by weighting with posterior probabilities:

$$E(U_3 | X) = E(U_3 | X, C=3) \cdot P(C=3 | X)$$

This statistically-principled approach effectively removes background and non-specific binding while preserving true biological signals.

---

## Key Features

- **Multi-channel Image Processing**: Pre-processing of multi-channel TIFF images via Bio-Formats  
- **GMM-based Background Correction**: Three-component Gaussian mixture models to identify background, non-specific signal, and biological signal  
- **Pixel- and Cell-level Adjustment**: Signal correction performed at pixel resolution with optional aggregation to cells. For Direct cell level adjustment, we recommend using 2 GMM components.   
- **Cross-sample Normalization**: Optional quantile normalization for cross-sample comparability  
- **Statistical Diagnostics**: Effect size calcualtion for checking stain quality  (Cohen’s d)  
- **Quality Control Visualization**: Heatmaps and summary plots to assess normalization performance  
- **HPC Compatibility**: Support for task array workflows for large-scale datasets  

---

## Installation

Install the development version from GitHub:

```r
remotes::install_github("Malvikakh/BgNormR")
```

### System Requirements

- R >= 4.5.0
- Java >= 1.8 (required for RBioFormats)

Before using the package, increase Java heap size:

```r
options(java.parameters = "-Xmx200g")
```

## Quick Start

### Processing Pipeline

#### GMM Normalization Pipeline:

1. **Read & Filter**: Load intensity data and apply median filtering
2. **Log Transform**: Apply log2(x / cofactor + 1) transformation
3. **Fit GMM**: Fit 3-component Gaussian mixture model to identify:
   - Component 1: Background pixels
   - Component 2: Non-specific binding pixels
   - Component 3: Signal pixels
4. **Classify**: Assign each pixel to a mixture component
5. **Adjust**: Apply variance-weighted deconvolution to remove background from singal
6. **Normalize**: Apply quantile normalization for cross-sample comparability
7. **Aggregate**: Compute per-cell median intensities


### High-Level Workflow (Recommended)

The easiest way to run a complete GMM normalization workflow:

```r
library(BgNorm)

# Set Java heap size
options(java.parameters = "-Xmx200g")

# Run complete workflow with one function
results <- run_gmm_workflow(
  image_file = "sample.tif",
  channel_names = c("CD3", "CD8", "CD20"),
  cell_mask = "segmentation.rds",
  output_dir = "results",
  create_qc = TRUE,
  save_results = TRUE
)

# View results
print(results)
summary(results)

# Access QC plots
print(results$qc_plots$cohen_d)
```

### GMM-based Normalization (Step-by-Step)

Process spatial imaging data with advanced Gaussian Mixture Model normalization:

```r
library(BgNorm)

# Set Java heap size
options(java.parameters = "-Xmx200g")

# For this step, the channel data and segmentation mask should be pre-loaded or available as rds objects.

# Process a single channel with GMM normalization
result <- processChannelGMM(
  file = "path/to/channel_data.rds",
  channel_name = "CD45",
  cell_mask = "path/to/segmentation_mask.rds",
  cofactor = 150,
  bit_depth = 16,
  G = 3  # Number of mixture components
)

# Access results
pixel_data <- result$pixel_data  # Per-pixel intensities and classifications
cell_data <- result$cell_data    # Per-cell median intensities
gmm_summary <- result$gmm_summary  # GMM fit parameters
```


