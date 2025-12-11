# spatialAtomizeR <img src="man/figures/logo.png" align="right" height="120" alt="" />

Bayesian Spatial Regression with Misaligned Data

## Overview

`spatialAtomizeR` implements Bayesian atom-based regression methods (ABRM) for assessing associations between spatially-misaligned variables, i.e., variables measured over two distinct and non-nested sets of spatial areas. The ABRM approach does not require any a priori re-alignment of the variables. This package uses Nimble under the hood for flexible and efficient Bayesian implementation. The package handles situations where:

- Outcome data and some covariates are measured on one spatial scale (called the "Y-grid"), while the remaining covariates are measured on a different spatial scale (called the "X-grid")
- The areas comprising the two spatial scales are misaligned, i.e., have mismatched boundaries, and neither scale is fully nested within the other
- Variables follow different distributions (normal, Poisson, binomial)

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bellayqian/spatialAtomizeR")
```

## Quick Start

### Basic Workflow

```r
library(spatialAtomizeR)
library(nimble)  # Required for ABRM models

# 1. Simulate misaligned spatial data with full parameter specification
sim_data <- simulate_misaligned_data(
  seed = 42,
  dist_covariates_x = c('normal', 'poisson', 'binomial'),
  dist_covariates_y = c('normal', 'poisson', 'binomial'),
  dist_y = 'poisson',
  x_intercepts = c(4, -1, -1),      # Intercepts for X covariates
  y_intercepts = c(4, -1, -1),      # Intercepts for Y covariates
  x_correlation = 0.5,              # Spatial correlation for X
  y_correlation = 0.5,              # Spatial correlation for Y
  beta0_y = -1,                     # Outcome intercept
  beta_x = c(-0.03, 0.1, -0.2),    # Coefficients for X covariates
  beta_y = c(0.03, -0.1, 0.2)      # Coefficients for Y covariates
)

# 2. Get NIMBLE model code
model_code <- get_abrm_model()

# 3. Run ABRM analysis
results <- run_abrm(
  gridx = sim_data$gridx,
  gridy = sim_data$gridy,
  atoms = sim_data$atoms,
  model_code = model_code,
  true_params = sim_data$true_params, # optional vector of true outcome model coefficient parameters
  norm_idx_x = 1,   # Index of normal-distributed X covariate
  pois_idx_x = 2,   # Index of Poisson-distributed X covariate
  binom_idx_x = 3,  # Index of binomial-distributed X covariate
  norm_idx_y = 1,   # Index of normal-distributed Y covariate
  pois_idx_y = 2,   # Index of Poisson-distributed Y covariate
  binom_idx_y = 3,  # Index of binomial-distributed Y covariate
  dist_y = 2,       # Outcome distribution: 1=normal, 2=poisson, 3=binomial
  niter = 50000,    # MCMC iterations
  nburnin = 30000,  # Burn-in iterations
  nchains = 2       # Number of chains
)

# 4. View results
print(results$parameter_estimates)
```

## Main Features

### Data Simulation
- Create two spatial grids ("X-grid" and "Y-grid") with non-nested spatial misalignment
- Generate synthetic spatially correlated variables with customizable distributions over each spatial grid
- Specify true parameter values for validation

### Model Fitting
- Atom-based Bayesian regression with NIMBLE
- Support for mixed-type variables (normal, Poisson, binomial)
- Multivariate CAR models to allow for information-sharing over space and across variables
- Automatic convergence diagnostics

### Method Comparison
- Compare ABRM with dasymetric mapping
- Calculate bias, RMSE, and coverage rates
- Generate comparison plots

### Sensitivity Analysis
- Test across different correlation structures
- Multiple simulations per setting
- Automated result summarization

### S3 Object System

All main functions return S3 objects with dedicated print, summary, and plot methods:
```r
# Create simulated data
sim_data <- simulate_misaligned_data(...)
class(sim_data)  # "misaligned_data"

# View results with clean formatting
print(sim_data)   # Clean overview
summary(sim_data) # Detailed information

# Run ABRM analysis
results <- run_abrm(...)
class(results)    # "abrm"

print(results)    # Shows parameter count, bias, coverage
summary(results)  # Shows full parameter table
plot(results)     # Shows MCMC diagnostic plots

# Compare methods
comparison <- run_both_methods(...)
class(comparison) # "abrm_comparison"

print(comparison)   # Shows method comparison summary
summary(comparison) # Shows detailed metrics by method
```

## S3 Methods Examples

The package provides intuitive S3 methods for all major output types:
```r
# Simulated data
sim_data <- simulate_misaligned_data(seed = 123, ...)
print(sim_data)
# Output:
# Simulated Misaligned Spatial Data
# ==================================
# Y-grid cells: 25
# X-grid cells: 100
# Atoms: 200
# ...

# ABRM results
results <- run_abrm(...)
print(results)
# Output:
# ABRM Model Results
# ==================
# Number of parameters estimated: 6
# Mean absolute bias: 0.0234
# Coverage rate: 95.00%
# Use summary() for detailed parameter estimates

summary(results)  # Shows full parameter table
```

## Key Functions

| Function | Description |
|----------|-------------|
| `simulate_misaligned_data()` | Generate simulated spatial data with full parameter control |
| `get_abrm_model()` | Get NIMBLE model specification |
| `run_abrm()` | Run ABRM analysis (wrapper function) |
| `run_nimble_model()` | Run NIMBLE MCMC with diagnostics |
| `run_both_methods()` | Compare ABRM and dasymetric mapping |
| `run_sensitivity_analysis()` | Conduct sensitivity analysis |
| `prepare_spatial_bookkeeping()` | Prepare spatial indices |
| `prepare_adjacency_matrices()` | Create spatial adjacency structures |
| `prepare_nimble_inputs()` | Prepare NIMBLE model inputs |

## Data Simulation Parameters

The `simulate_misaligned_data()` function accepts the following parameters:

**Reproducibility Parameters:**
- `seed`: Random seed for reproducibility

**Covariate Distributions:**
- `dist_covariates_x`: Vector of distribution types for X-grid covariates (e.g., `c('normal', 'poisson', 'binomial')`)
- `dist_covariates_y`: Vector of distribution types for Y-grid covariates
- `dist_y`: Distribution type for outcome variable (`'normal'`, `'poisson'`, or `'binomial'`)

**Data Generation Parameters:**
- `x_intercepts`: Intercepts for X-grid covariates (length must match `dist_covariates_x`)
- `y_intercepts`: Intercepts for Y-grid covariates (length must match `dist_covariates_y`)
- `beta0_y`: Intercept for the outcome model
- `beta_x`: True coefficients for X-grid covariates in outcome model
- `beta_y`: True coefficients for Y-grid covariates in outcome model

**Between-Variable Correlation:**
- `x_correlation`: Correlation between X-grid covariates (0 to 1)
- `y_correlation`: Correlation between Y-grid covariates (0 to 1)

## Distribution Type Indices

When running ABRM models, you need to specify which covariates follow which distributions:

- `norm_idx_x`, `norm_idx_y`: Indices of normally-distributed covariates
- `pois_idx_x`, `pois_idx_y`: Indices of Poisson-distributed covariates
- `binom_idx_x`, `binom_idx_y`: Indices of binomially-distributed covariates
- `dist_y`: Outcome distribution type (1=normal, 2=poisson, 3=binomial)

**Example:** If `dist_covariates_x = c('normal', 'poisson', 'binomial')`, then:
- `norm_idx_x = 1` (first covariate)
- `pois_idx_x = 2` (second covariate)
- `binom_idx_x = 3` (third covariate)

## Example: Comprehensive Sensitivity Analysis

```r
library(spatialAtomizeR)
library(nimble)

# Define base parameters
base_params <- list(
  dist_covariates_x = c('normal','poisson','binomial'),
  dist_covariates_y = c('normal','poisson','binomial'),
  dist_y = 'poisson',
  x_intercepts = c(4, -1, -1),
  y_intercepts = c(4, -1, -1),
  beta0_y = -1,
  beta_x = c(-0.03, 0.1, -0.2),
  beta_y = c(0.03, -0.1, 0.2)
)

# Get model code
model_code <- get_abrm_model()

# Run sensitivity analysis across correlation structures
sensitivity_results <- run_sensitivity_analysis(
  correlation_grid = c(0.2, 0.6),
  n_sims_per_setting = 3,
  base_params = base_params,
  model_code = model_code,
  base_seed = 123
)

# View summary by correlation
print(sensitivity_results$summary_by_correlation)

# Access detailed results
write.csv(
  sensitivity_results$combined_results,
  "sensitivity_analysis_full_results.csv"
)
```

## Requirements

- R >= 4.0.0
- **nimble** for MCMC sampling (must be loaded)
- Spatial packages: **sp**, **sf**, **spdep**, **raster**
- **BiasedUrn** for multivariate hypergeometric sampling
- MASS for multivariate normal generation
- dplyr, tidyr for data manipulation

## Funding and Project Information

This work was funded by the Robert Wood Johnson Foundation, Grant 81746. Project details are provided below.

**Project Title:** Aligning spatially misaligned data for health equity analysis, action, and accountability

**Principal Investigators:** Dr. Nancy Krieger (PI) and Dr. Rachel Nethery (co-PI)

**Start Date:** July 2024

**Project Team and Collaborators:**
- Yunzhe Qian (Bella), MS (Research Assistant, Dept of Biostatistics, HSPH)
- Rachel Nethery, PhD (Associate Professor, Dept of Biostatistics, HSPH)
- Nancy Krieger, PhD (Professor, Department of Social and Behavioral Sciences (SBS), HSPH)
- Nykesha Johnson, MPH (Statistical Data Analyst/Data Manager, SBS, HSPH)

## Citation

If you use this package, please cite:

Qian Y, Nethery R, Krieger N, Johnson N (2025). spatialAtomizeR: Spatial Analysis with Misaligned Data Using Atom-Based Regression Models. R package version 0.2.4, https://github.com/bellayqian/spatialAtomizeR.

## About

This work is an extension of:

Nethery, R. C., Testa, C., Tabb, L. P., Hanage, W. P., Chen, J. T., & Krieger, N. (2023). Addressing spatial misalignment in population health research: a case study of US congressional district political metrics and county health data. MedRxiv.

**Spatial misalignment**—which occurs when data on multiple variables are collected using mismatched geographic boundary definitions—is a longstanding challenge in public health research. For instance, congressional districts can cut across multiple counties, and environmental hazard zones may cross census tract boundaries, in both cases creating intersecting areas that complicate efforts to study the relationships between health outcomes and their social, political, and environmental determinants.

**Atom-based regression models (ABRM)** offer a promising alternative by using atoms—the intersecting areas of all relevant units—as the fundamental units of analysis. By preserving the original spatial resolution of the data, ABRM account for uncertainty in statistical relationships while offering a robust method for handling misaligned data.

## Getting Help

- **Report bugs**: https://github.com/bellayqian/spatialAtomizeR/issues
- **Documentation**: `?spatialAtomizeR`
- **Vignette**: `vignette("getting-started", package = "spatialAtomizeR")`
- **Function help**: `?simulate_misaligned_data`, `?run_abrm`, etc.

## License

MIT License
