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
  true_params = sim_data$true_params, # optional: true coefficients (see ?abrm_params)
  x_types = c("norm", "pois", "binom"),  # distribution of each X covariate
  y_types = c("norm", "pois", "binom"),  # distribution of each Y covariate
  dist_y = "pois",  # outcome distribution: "norm", "pois" or "binom"
  niter = 50000,    # MCMC iterations
  nburnin = 30000,  # Burn-in iterations
  nchains = 2       # Number of chains
)

# 4. View results
print(results$parameter_estimates)
plot_trace(results)         # trace plots (a ggplot object)
plot_density(results)       # posterior density plots (a ggplot object)
mcmc_diagnostics(results)   # R-hat and effective sample sizes
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

## Covariate Distributions and Column Names

When running ABRM models, say which distribution each covariate follows:

- `x_types`: character vector with one element per X-grid covariate
- `y_types`: character vector with one element per Y-grid covariate
- `dist_y`: outcome distribution

Each takes the values `"norm"`, `"pois"` or `"binom"`.

**Example:** If `dist_covariates_x = c("norm", "pois", "binom")`, then
`x_types = c("norm", "pois", "binom")`.

By default `run_abrm()` looks for columns named `covariate_x_1`, `covariate_x_2`,
..., `covariate_y_1`, ... and `y` (outcome). To use your own column names, pass
them with `x_vars`, `y_vars` and `outcome`; `id_x`, `id_y`, `atom_id_x`,
`atom_id_y` and `population` name the ID and population columns. Grid IDs can be
integer, numeric or character.

The `norm_idx_*`/`pois_idx_*`/`binom_idx_*` arguments and numeric `dist_y`
(1, 2, 3) from earlier versions still work but are deprecated.

## Creating Atoms

Atoms are the pieces of the intersection of the two grids; `create_atoms(gridy, gridx)`
computes them with the **sf** package alone.

## True Parameter Values

To compare estimates against known coefficients, build `true_params` with
`abrm_params(beta0_y = , beta_x = , beta_y = )`. Vectors can be unnamed (matched by
position) or named (matched by covariate name).

## Requirements

- R >= 4.1.0
- **nimble** for MCMC sampling (must be loaded)
- Spatial packages: **sp**, **sf**, **spdep**
- **BiasedUrn** for multivariate hypergeometric sampling
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

Qian Y, Nethery R, Krieger N, Johnson N (2026). spatialAtomizeR: Spatial Analysis with Misaligned Data Using Atom-Based Regression Models. R package version 0.2.8, https://github.com/bellayqian/spatialAtomizeR.

## About

Atom-based regression models (ABRM) were originally proposed by:

Mugglin, A. S., Carlin, B. P., & Gelfand, A. E. (2000). Fully model-based approaches for spatially misaligned data. Journal of the American Statistical Association, 95(451), 877-887.

The methods were further investigated in:

Trevisani, M., & Gelfand, A. E. (2013). Spatial misalignment models for small area estimation: A simulation study. In Advances in Theoretical and Applied Statistics (pp. 269-279). Springer.

This package builds upon the applied work in:

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
