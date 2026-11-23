# Run this locally, NOT inside the vignette
# e.g., in data-raw/vignette_cache.R

library(spatialAtomizeR)
set.seed(500)

sim_data <- simulate_misaligned_data(
  seed = 42,
  
  # Distribution specifications for covariates
  dist_covariates_x = c('norm', 'pois', 'binom'),
  dist_covariates_y = c('norm', 'pois', 'binom'),
  dist_y = 'pois',  # Outcome distribution
  
  # Intercepts for data generation (REQUIRED)
  x_intercepts = c(4, -1, -1),      # One per X covariate
  y_intercepts = c(4, -1, -1),      # One per Y covariate
  beta0_y = -1,                     # Outcome model intercept
  
  # Spatial correlation parameters
  x_correlation = 0.5,  # Correlation between X covariates
  y_correlation = 0.5,  # Correlation between Y covariates
  
  # True effect sizes for outcome model
  beta_x = c(-0.03, 0.1, -0.2),    # Effects of X covariates
  beta_y = c(0.03, -0.1, 0.2)      # Effects of Y covariates
)

model_code <- get_abrm_model()

abrm_results <- run_abrm(
  gridx = sim_data$gridx,
  gridy = sim_data$gridy,
  atoms = sim_data$atoms,
  model_code = model_code,
  true_params = sim_data$true_params,  # Optional: for validation
  
  # Distribution of each covariate (same order as dist_covariates_x/y)
  x_types = c("norm", "pois", "binom"),
  y_types = c("norm", "pois", "binom"),
  
  # Outcome distribution: "norm", "pois" or "binom"
  dist_y = "pois",
  
  # MCMC parameters
  niter = 50000,    # Total iterations per chain
  nburnin = 30000,  # Burn-in iterations
  nchains = 2,      # Number of chains
  seed = 123,
  compute_waic = TRUE
)

# Verify fitted() works BEFORE saving
fitted(abrm_results)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(abrm_results, "inst/extdata/abrm_vignette_results.rds")
saveRDS(sim_data,     "inst/extdata/sim_vignette_data.rds")
