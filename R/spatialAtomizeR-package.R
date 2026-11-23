#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Declare global variables used in NIMBLE models to avoid R CMD check NOTEs
utils::globalVariables(c(
  # NIMBLE model parameters and indices
  "D", "J_x", "J_y", "S_x", "S_y", "p_x", "p_y", "p_a",
  "covars_atom_x", "covars_atom_y", "beta_a_x", "beta_a_yx",
  "norm_n_x", "norm_n_y", "pois_n_x", "pois_n_y", "binom_n_x", "binom_n_y",
  "norm_idx_x", "norm_idx_y", "pois_idx_x", "pois_idx_y", "binom_idx_x", "binom_idx_y",
  "vartype_y",
  
  # Model parameters
  "beta_0_x", "beta_0_y", "beta_0_yx", "beta_y",
  "sigma2_x", "sigma2_y", "sigma2_yx",
  "tau_x", "tau_y", "tau_yx",
  "phi_x", "phi_y", "phi_yx",
  "psi_x", "psi_yx",
  "Prec_x", "Prec_yx",
  
  # Spatial structure
  "expand_x", "expand_y", "xlatent_ind", "ylatent_ind",
  "pop_atoms_x", "pop_atoms_y", "pop_atoms_outcome", "y_to_atom",
  
  # Data arrays
  "x", "x_reorder", "x_latent", "yx_obs", "yx_latent", "y_obs",
  
  # Working variables
  "w_x", "w_y", "w_yx",
  
  # NIMBLE functions
  "inverse", "inprod",
  
  # NIMBLE custom distribution functions
  "dmfnchypg", "Rmfnchypg", "rmfnchypg", "returnType"
))

#' spatialAtomizeR: Bayesian Spatial Regression with Misaligned Data
#'
#' Implements atom-based Bayesian regression methods (ABRM) for spatial data 
#' with misaligned grids.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{simulate_misaligned_data}}}{Generate simulated spatial data}
#'   \item{\code{\link{create_atoms}}}{Intersect two grids to create atoms}
#'   \item{\code{\link{get_abrm_model}}}{Get NIMBLE model code for ABRM}
#'   \item{\code{\link{run_abrm}}}{Run atom-based Bayesian regression model}
#'   \item{\code{\link{abrm_params}}}{Specify true coefficient values for \code{run_abrm()}}
#'   \item{\code{\link{plot_trace}}, \code{\link{plot_density}}}{MCMC trace and density plots}
#'   \item{\code{\link{mcmc_diagnostics}}}{R-hat and effective sample sizes}
#' }
#'
#' @name spatialAtomizeR-package
#' @aliases spatialAtomizeR
NULL

# Create package environment
.pkg_env <- new.env(parent = emptyenv())

# Package load hook
.onLoad <- function(libname, pkgname) {
  # Initialize the environment variable
  .pkg_env$distributions_registered <- FALSE
}

# Package attach hook
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("spatialAtomizeR loaded. Use get_abrm_model() to access the model.")
}