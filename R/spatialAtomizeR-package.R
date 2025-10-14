#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' spatialAtomizeR: Bayesian Spatial Regression with Misaligned Data
#'
#' Implements atom-based Bayesian regression methods (ABRM) for spatial data 
#' with misaligned grids.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{simulate_misaligned_data}}}{Generate simulated spatial data}
#'   \item{\code{\link{get_abrm_model}}}{Get NIMBLE model code for ABRM}
#'   \item{\code{\link{run_abrm}}}{Run atom-based Bayesian regression model}
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