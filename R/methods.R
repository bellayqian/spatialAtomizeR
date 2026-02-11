#' Print method for abrm objects
#' @param x An abrm object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the ABRM model results including convergence status, number of parameters estimated, and key fit statistics.
#' @export
print.abrm <- function(x, ...) {
  cat("\nABRM Model Results\n")
  cat("==================\n")
  cat("Number of parameters estimated:", nrow(x$parameter_estimates), "\n")
  if (is.data.frame(x$parameter_estimates) && "true_beta" %in% names(x$parameter_estimates)) {
    cat("Mean absolute bias:", 
        round(mean(abs(x$parameter_estimates$bias), na.rm = TRUE), 4), "\n")
    cat("Coverage rate:", 
        round(mean(x$parameter_estimates$within_ci, na.rm = TRUE) * 100, 2), "%\n")
  }
  cat("\nUse summary() for detailed parameter estimates\n")
  invisible(x)
}

#' Summary method for abrm objects
#' @param object An abrm object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{object}. The function is called for its side effect of printing the ABRM model summary including detailed parameter estimates.
#' @export
summary.abrm <- function(object, ...) {
  cat("ABRM Model Summary\n")
  cat("==================\n\n")
  # Remove redundant row names in parameter estimates
  param_table <- object$parameter_estimates
  rownames(param_table) <- NULL
  
  # Better column names
  if("estimated_beta" %in% names(param_table)) {
    names(param_table)[names(param_table) == "estimated_beta"] <- "Estimate"
  }
  if("std_error" %in% names(param_table)) {
    names(param_table)[names(param_table) == "std_error"] <- "Std.Error"
  }
  if("ci_lower" %in% names(param_table)) {
    names(param_table)[names(param_table) == "ci_lower"] <- "CI.Lower"
  }
  if("ci_upper" %in% names(param_table)) {
    names(param_table)[names(param_table) == "ci_upper"] <- "CI.Upper"
  }
  
  print(param_table, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' Plot method for abrm objects
#' 
#' @param x An object of class "abrm"
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of displaying MCMC diagnostic plots (trace plots and density plots) if they are available in the abrm object.
#' @importFrom graphics par
#' @export
plot.abrm <- function(x, ...) {
  if(!is.null(x$mcmc_results$convergence$plots)) {
    print(x$mcmc_results$convergence$plots$trace)
    print(x$mcmc_results$convergence$plots$density)
  } else {
    message("No diagnostic plots available.\n")
  }
  invisible(x)
}

#' Print method for misaligned_data objects
#' @param x A misaligned_data object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the simulated misaligned spatial data including grid dimensions and number of atoms.
#' @export
print.misaligned_data <- function(x, ...) {
  message("Simulated Misaligned Spatial Data\n")
  message("==================================\n")
  message("Y-grid cells: ", nrow(x$gridy), "\n")
  message("X-grid cells: ", nrow(x$gridx), "\n")
  message("Atoms: ", nrow(x$atoms), "\n")
  message("Number of X covariates: ", length(grep("covariate_x", names(x$gridx))), "\n")
  message("Number of Y covariates: ", length(grep("covariate_y", names(x$gridy))), "\n")
  if (!is.null(x$true_params)) {
    message("True parameters available\n")
  }
  invisible(x)
}

#' Summary method for misaligned_data objects
#' @param object A misaligned_data object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{object}. The function is called for its side effect of printing the misaligned data summary including grid information and true parameter values (beta_x and beta_y).
#' @export
summary.misaligned_data <- function(object, ...) {
  print.misaligned_data(object)
  message("True beta_x: ", paste(object$true_params$beta_x, collapse = ", "), "\n")
  message("True beta_y: ", paste(object$true_params$beta_y, collapse = ", "), "\n")
  invisible(object)
}

#' Variance-covariance method for abrm objects
#' Extracts variance-covariance matrices for regression coefficients from 
#' MCMC posterior samples. Returns separate matrices for X-grid and Y-grid
#' coefficients.
#' @param object An abrm object from run_abrm()
#' @param ... Additional arguments (unused)
#' @return A list with class "vcov.abrm" containing:
#'   \item{vcov_beta_x}{Variance-covariance matrix for X-grid coefficients}
#'   \item{vcov_beta_y}{Variance-covariance matrix for Y-grid coefficients}
#'   \item{vcov_beta_0}{Variance of the intercept (scalar)}
#'   \item{vcov_all}{Combined variance-covariance matrix for all parameters}
#' @details 
#' The variance-covariance matrices are computed from the posterior samples
#' of the MCMC chains. If multiple chains were run, samples are combined
#' across chains before computing covariances.
#' 
#' @examples
#' \dontrun{
#' # Fit model
#' results <- run_abrm(...)
#' 
#' # Get variance-covariance matrices
#' vcov_mats <- vcov(results)
#' 
#' # Access specific matrices
#' vcov_mats$vcov_beta_x  # Covariance for X-grid coefficients
#' vcov_mats$vcov_beta_y  # Covariance for Y-grid coefficients
#' 
#' # Compute standard errors from diagonal
#' sqrt(diag(vcov_mats$vcov_beta_x))
#' 
#' # Compute correlation matrix
#' cov2cor(vcov_mats$vcov_beta_y)
#' }
#' @importFrom stats var cov
#' @export
vcov.abrm <- function(object, ...) {
  
  # Check if MCMC samples are available
  if (is.null(object$mcmc_results) || is.null(object$mcmc_results$samples)) {
    stop("MCMC samples not found in abrm object. Cannot compute variance-covariance matrix.")
  }
  
  # Extract MCMC samples
  samples <- object$mcmc_results$samples
  
  # Convert to matrix - handle different storage formats
  if (inherits(samples, "mcmc.list")) {
    # coda mcmc.list object
    samples_matrix <- do.call(rbind, samples)
  } else if (inherits(samples, "mcmc")) {
    # Single coda mcmc object
    samples_matrix <- as.matrix(samples)
  } else if (is.list(samples) && !is.data.frame(samples)) {
    # List of chains (chain1, chain2, etc.) - YOUR FORMAT
    # Check if list elements are matrices
    if (all(sapply(samples, is.matrix))) {
      samples_matrix <- do.call(rbind, samples)
    } else {
      stop("List elements are not matrices")
    }
  } else if (is.matrix(samples)) {
    # Already a matrix
    samples_matrix <- samples
  } else {
    stop("Unexpected format for MCMC samples: ", class(samples))
  }
  
  # Get parameter names
  param_names <- colnames(samples_matrix)
  
  if (is.null(param_names)) {
    stop("MCMC samples do not have column names. Cannot identify parameters.")
  }
  
  # Identify parameter indices
  beta_0_idx <- grep("^beta_0_y$|^beta0_y$", param_names)
  beta_x_idx <- grep("^beta_x\\[", param_names)
  beta_y_idx <- grep("^beta_y\\[", param_names)
  
  # Create list to store results
  vcov_list <- list()
  
  # Compute variance for intercept (scalar)
  if (length(beta_0_idx) > 0) {
    vcov_list$vcov_beta_0 <- var(samples_matrix[, beta_0_idx, drop = TRUE])
    names(vcov_list$vcov_beta_0) <- param_names[beta_0_idx]
  } else {
    vcov_list$vcov_beta_0 <- NULL
    warning("No intercept parameter (beta_0_y) found in MCMC samples")
  }
  
  # Compute variance-covariance for X-grid coefficients
  if (length(beta_x_idx) > 0) {
    vcov_list$vcov_beta_x <- cov(samples_matrix[, beta_x_idx, drop = FALSE])
    rownames(vcov_list$vcov_beta_x) <- param_names[beta_x_idx]
    colnames(vcov_list$vcov_beta_x) <- param_names[beta_x_idx]
  } else {
    vcov_list$vcov_beta_x <- NULL
    warning("No X-grid coefficients (beta_x[]) found in MCMC samples")
  }
  
  # Compute variance-covariance for Y-grid coefficients
  if (length(beta_y_idx) > 0) {
    vcov_list$vcov_beta_y <- cov(samples_matrix[, beta_y_idx, drop = FALSE])
    rownames(vcov_list$vcov_beta_y) <- param_names[beta_y_idx]
    colnames(vcov_list$vcov_beta_y) <- param_names[beta_y_idx]
  } else {
    vcov_list$vcov_beta_y <- NULL
    warning("No Y-grid coefficients (beta_y[]) found in MCMC samples")
  }
  
  # Compute combined variance-covariance for all beta parameters
  all_beta_idx <- c(beta_0_idx, beta_x_idx, beta_y_idx)
  if (length(all_beta_idx) > 0) {
    vcov_list$vcov_all <- cov(samples_matrix[, all_beta_idx, drop = FALSE])
    rownames(vcov_list$vcov_all) <- param_names[all_beta_idx]
    colnames(vcov_list$vcov_all) <- param_names[all_beta_idx]
  } else {
    vcov_list$vcov_all <- NULL
    warning("No beta parameters found in MCMC samples")
  }
  
  # Add class for potential future print/summary methods
  class(vcov_list) <- c("vcov.abrm", "list")
  
  return(vcov_list)
}

#' Print method for vcov.abrm objects
#' @param x A vcov.abrm object
#' @param ... Additional arguments (currently unused)
#' @return Invisibly returns the input object
#' @export
print.vcov.abrm <- function(x, ...) {
  cat("\nVariance-Covariance Matrices for ABRM Model\n")
  cat("============================================\n\n")
  
  if (!is.null(x$vcov_beta_0)) {
    cat("Intercept Variance:\n")
    cat("  ", names(x$vcov_beta_0), ": ", round(x$vcov_beta_0, 6), "\n\n", sep = "")
  }
  
  if (!is.null(x$vcov_beta_x)) {
    cat("X-Grid Coefficients Variance-Covariance Matrix:\n")
    cat("  Dimensions:", nrow(x$vcov_beta_x), "x", ncol(x$vcov_beta_x), "\n")
    print(round(x$vcov_beta_x, 6))
    cat("\n")
  }
  
  if (!is.null(x$vcov_beta_y)) {
    cat("Y-Grid Coefficients Variance-Covariance Matrix:\n")
    cat("  Dimensions:", nrow(x$vcov_beta_y), "x", ncol(x$vcov_beta_y), "\n")
    print(round(x$vcov_beta_y, 6))
    cat("\n")
  }
  
  cat("Use vcov_object$vcov_beta_x, vcov_object$vcov_beta_y, etc. to access matrices\n")
  
  invisible(x)
}