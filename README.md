# spatialAtomizeR

Bayesian Spatial Regression with Misaligned Data

## Overview

`spatialAtomizeR` implements atom-based Bayesian regression methods (ABRM) for spatial data with misaligned grids. The package handles situations where:

- Outcome data and covariates are available on different spatial grids
- The spatial grids are non-nested (misaligned)
- Variables follow different distributions (normal, Poisson, binomial)

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bellayqian/spatialAtomizeR")
```

## Quick Start

```r
library(spatialAtomizeR)

# 1. Simulate misaligned spatial data
sim_data <- simulate_misaligned_data(
  seed = 42,
  res1 = c(5, 5),  # Y grid resolution
  res2 = c(10, 10), # X grid resolution
  dist_covariates_x = c('normal', 'poisson', 'binomial'),
  dist_covariates_y = c('normal', 'poisson', 'binomial'),
  dist_y = 'poisson',
  x_correlation = 0.5,
  y_correlation = 0.5,
  beta_x = c(-0.03, 0.1, -0.2),
  beta_y = c(0.03, -0.1, 0.2)
)

# 2. Get NIMBLE model code
model_code <- get_abrm_model()

# 3. Run ABRM
results <- run_abrm(
  sim_data = sim_data,
  model_code = model_code,
  norm_idx_x = 1,
  pois_idx_x = 2,
  binom_idx_x = 3,
  norm_idx_y = 1,
  pois_idx_y = 2,
  binom_idx_y = 3,
  dist_y = 2,  # 1=normal, 2=poisson, 3=binomial
  niter = 50000,
  nburnin = 30000,
  nchains = 2
)

# 4. View results
print(results$parameter_estimates)
```

## Main Features

### Data Simulation
- Generate spatially correlated variables with customizable distributions
- Create non-nested misaligned spatial grids
- Specify true parameter values for validation

### Model Fitting
- Atom-based Bayesian regression with NIMBLE
- Support for mixed-type variables (normal, Poisson, binomial)
- Multivariate CAR models for spatial correlation
- Automatic convergence diagnostics

### Method Comparison
- Compare ABRM with dasymetric mapping
- Calculate bias, RMSE, and coverage rates
- Generate comparison plots

### Sensitivity Analysis
- Test across different correlation structures
- Multiple simulations per setting
- Automated result summarization

## Key Functions

| Function | Description |
|----------|-------------|
| `simulate_misaligned_data()` | Generate simulated spatial data |
| `get_abrm_model()` | Get NIMBLE model specification |
| `run_abrm()` | Run ABRM analysis |
| `run_both_methods()` | Compare ABRM and dasymetric mapping |
| `run_sensitivity_analysis()` | Conduct sensitivity analysis |
| `prepare_spatial_bookkeeping()` | Prepare spatial indices |
| `prepare_nimble_inputs()` | Prepare NIMBLE model inputs |

## Example Workflow

```r
# Sensitivity analysis across correlation structures
sensitivity_results <- run_sensitivity_analysis(
  correlation_grid = c(0.2, 0.6),
  n_sims_per_setting = 3,
  model_code = get_abrm_model(),
  base_seed = 123
)

# View summary
print(sensitivity_results$summary_by_correlation)
```

## Requirements

- R >= 4.0.0
- NIMBLE for MCMC sampling
- Spatial packages: sp, sf, spdep
- BiasedUrn for multivariate hypergeometric sampling

## Funding and Project Information

This work was funded by the Robert Wood Johnson Foundation, Grant 81746. Project details are provided below.

**Project Title:** Aligning spatially misaligned data for health equity analysis, action, and accountability

**Principal Investigators:** Dr. Nancy Krieger (PI) and Dr. Rachel Nethery (co-PI)

**Start Date:** July 2024

**Project Team and Collaborators:**
- Yunzhe Qian (Bella), MS (Research Assistant, Dept of Biostatistics, HSPH)
- Rachel Nethery, PhD (Assistant Professor, Dept of Biostatistics, HSPH)
- Nancy Krieger, PhD (Professor, Department of Social and Behavioral Sciences (SBS), HSPH)
- Nykesha Johnson, MPH (Statistical Data Analyst/Data Manager, SBS, HSPH)

## Citation

If you use this package, please cite:

Qian, Y., & Nethery, R. (2025). spatialAtomizeR: Atom-Based Regression Models for Misaligned Spatial Data. R package version 0.1.0.

## About

This work is an extension of:

Nethery, R. C., Testa, C., Tabb, L. P., Hanage, W. P., Chen, J. T., & Krieger, N. (2023). Addressing spatial misalignment in population health research: a case study of US congressional district political metrics and county health data. MedRxiv.

Spatial misalignment—which occurs when data on multiple variables are collected using mismatched geographic boundary definitions—is a longstanding challenge in public health research. For instance, congressional districts can cut across multiple counties, and environmental hazard zones may cross census tract boundaries, in both cases creating intersecting areas that complicate efforts to study the relationships between health outcomes and their social, political, and environmental determinants.

Atom-based regression models (ABRM) offer a promising alternative by using atoms, the intersecting areas of all relevant units, as the fundamental units of analysis. By preserving the original spatial resolution of the data, ABRM account for uncertainty in statistical relationships while offering a robust method for handling misaligned data.

## Getting Help

- Report bugs at: https://github.com/bellayqian/spatialAtomizeR/issues
- Check documentation: `?spatialAtomizeR`
- See vignette: `vignette("getting-started", package = "spatialAtomizeR")`

## License

MIT License