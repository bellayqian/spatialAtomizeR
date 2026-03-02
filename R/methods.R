#' Print Method for Misaligned_data Objects
#' @param x A misaligned_data object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the simulated misaligned spatial data including grid dimensions and number of atoms.
#' @examples
#' ## Simulate misaligned spatial data
#' sim_data <- simulate_misaligned_data(
#'   seed              = 1,
#'   dist_covariates_x = c("normal", "poisson"),
#'   dist_covariates_y = c("normal", "poisson"),
#'   dist_y            = "poisson",
#'   x_intercepts      = c(0, -1),
#'   y_intercepts      = c(0, -1),
#'   beta0_y           = -1,
#'   beta_x            = c(0.1, -0.05),
#'   beta_y            = c(-0.1, 0.05)
#' )
#' print(sim_data)
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

#' Summary Method for Misaligned Data Objects
#'
#' Prints grid dimensions, atom counts, and the true regression coefficients
#' (\code{beta_x} and \code{beta_y}) used to generate the data.
#'
#' @param object An object of class \code{"misaligned_data"} returned by
#'   \code{\link{simulate_misaligned_data}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{object}. Called for its side effect of
#'   printing the data summary including the true parameter values.
#'
#' @examples
#' ## Simulate misaligned spatial data
#' sim_data <- simulate_misaligned_data(
#'   seed              = 1,
#'   dist_covariates_x = c("normal", "poisson"),
#'   dist_covariates_y = c("normal", "poisson"),
#'   dist_y            = "poisson",
#'   x_intercepts      = c(0, -1),
#'   y_intercepts      = c(0, -1),
#'   beta0_y           = -1,
#'   beta_x            = c(0.1, -0.05),
#'   beta_y            = c(-0.1, 0.05)
#' )
#' summary(sim_data)
#'
#' @export
summary.misaligned_data <- function(object, ...) {
  print.misaligned_data(object)
  message("True beta_x: ", paste(object$true_params$beta_x, collapse = ", "), "\n")
  message("True beta_y: ", paste(object$true_params$beta_y, collapse = ", "), "\n")
  invisible(object)
}

#' Plot Misaligned Data Object
#'
#' Visualizes the spatial layout of a \code{misaligned_data} object. By default,
#' both grids are overlaid to illustrate the misalignment between them.
#'
#' @param x A \code{misaligned_data} object.
#' @param which Character string specifying what to plot. One of \code{"both"}
#'   (default), \code{"gridy"}, \code{"gridx"}, or \code{"atoms"}.
#' @param color_y Color for the Y-grid boundaries. Default \code{"blue"}.
#' @param color_x Color for the X-grid boundaries. Default \code{"orange"}.
#' @param title Optional character string for the plot title. If \code{NULL},
#'   a default title is used.
#' @param ... Additional arguments passed to \code{ggplot2::geom_sf()} when
#'   \code{which} is not \code{"both"}.
#'
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' ## Simulate misaligned spatial data
#' sim_data <- simulate_misaligned_data(
#'   seed              = 1,
#'   dist_covariates_x = c("normal", "poisson"),
#'   dist_covariates_y = c("normal", "poisson"),
#'   dist_y            = "poisson",
#'   x_intercepts      = c(0, -1),
#'   y_intercepts      = c(0, -1),
#'   beta0_y           = -1,
#'   beta_x            = c(0.1, -0.05),
#'   beta_y            = c(-0.1, 0.05)
#' )
#' # Default: overlay both grids to visualise misalignment
#' plot(sim_data)
#'
#' # Plot a single component
#' plot(sim_data, which = "gridy")
#' plot(sim_data, which = "gridx")
#' plot(sim_data, which = "atoms")
#'
#' # Custom colours
#' plot(sim_data, color_y = "steelblue", color_x = "firebrick")
#'
#' # Custom title
#' plot(sim_data, title = "Utah Spatial Misalignment")
#'
#' @importFrom ggplot2 ggplot geom_sf labs theme_void
#' @export
plot.misaligned_data <- function(x, which   = c("both", "gridy", "gridx", "atoms"),
                                 color_y = "blue", color_x = "orange", title   = NULL, ...) {
  
  which <- match.arg(which)
  
  if (which == "both") {
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        data      = x[["gridx"]],
        fill      = NA,
        color     = color_x,
        linewidth = 1.2
      ) +
      ggplot2::geom_sf(
        data      = x[["gridy"]],
        fill      = NA,
        color     = color_y,
        linetype  = "dashed",
        linewidth = 0.5
      ) +
      ggplot2::labs(
        title    = if (!is.null(title)) title else "Spatial Misalignment",
        subtitle = paste0("Orange: gridx  |  Blue (dashed): gridy")
      ) +
      ggplot2::theme_void()
    
  } else {
    
    sf_obj <- x[[which]]
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        data  = sf_obj,
        fill  = NA,
        color = switch(which, gridy = color_y, gridx = color_x, "black"),
        ...
      ) +
      ggplot2::labs(
        title = if (!is.null(title)) title else paste("Spatial layout:", which)
      ) +
      ggplot2::theme_void()
    
  }
  
  print(p)
  invisible(x)
}


#' Summary Method for ABRM Objects
#'
#' Returns a structured summary of a fitted ABRM, including posterior means,
#' standard deviations, and 95% credible intervals for all regression
#' coefficients. When true parameter values are available (i.e., the model
#' was fitted on simulated data via \code{\link{simulate_misaligned_data}}),
#' mean absolute bias and credible-interval coverage are also reported.
#'
#' @param object An object of class \code{"abrm"} returned by
#'   \code{\link{run_abrm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.abrm"} (a list) with components:
#'   \describe{
#'     \item{\code{n_params}}{Number of regression parameters estimated.}
#'     \item{\code{mean_abs_bias}}{Mean absolute bias across parameters, or
#'       \code{NA} if true values were not supplied.}
#'     \item{\code{coverage}}{Percentage of parameters whose true value falls
#'       within its 95\% credible interval, or \code{NA} if unavailable.}
#'     \item{\code{estimates}}{Data frame of posterior summaries.}
#'   }
#'   The object is printed via \code{\link{print.summary.abrm}}.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "normal", dist_covariates_y = "normal",
#'     dist_y = "normal", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), norm_idx_x = 1, norm_idx_y = 1,
#'     dist_y = 1, niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   summary(results)    # full parameter table with bias and coverage
#' }
#'
#' @export
summary.abrm <- function(object, ...) {
  retlist <- list(
    n_params     = nrow(object$parameter_estimates),
    has_true     = "true_beta" %in% names(object$parameter_estimates),
    mean_abs_bias = if ("bias" %in% names(object$parameter_estimates))
      round(mean(abs(object$parameter_estimates$bias), na.rm = TRUE), 4)
    else NA,
    coverage     = if ("within_ci" %in% names(object$parameter_estimates))
      round(mean(object$parameter_estimates$within_ci, na.rm = TRUE) * 100, 2)
    else NA,
    estimates    = object$parameter_estimates
  )
  class(retlist) <- c("summary.abrm", "list")
  retlist
}

#' Print Method for summary.abrm Objects
#'
#' Prints the full parameter table from a \code{\link{summary.abrm}} object,
#' including posterior means, standard deviations, credible intervals, and
#' (when available) bias and coverage.
#'
#' @param x   An object of class \code{"summary.abrm"} returned by
#'   \code{\link{summary.abrm}}.
#' @param digits Integer; number of significant digits. Default \code{4}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "normal", dist_covariates_y = "normal",
#'     dist_y = "normal", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), norm_idx_x = 1, norm_idx_y = 1,
#'     dist_y = 1, niter = 1000, nburnin = 500, nchains = 2, seed = 1,
#'     save_plots = FALSE
#'   )
#'   s <- summary(results)
#'   print(s)              # same as typing summary(results) interactively
#'   print(s, digits = 2)  # fewer decimal places
#' }
#'
#' @export
print.summary.abrm <- function(x, digits = 4, ...) {
  cat("ABRM Model Summary\n")
  cat("==================\n")
  cat("Parameters estimated:", x$n_params, "\n")
  if (!is.na(x$mean_abs_bias)) cat("Mean absolute bias: ", x$mean_abs_bias, "\n")
  if (!is.na(x$coverage))      cat("Coverage rate: ",      x$coverage, "%\n\n")
  print(x$estimates, row.names = FALSE, digits = digits)
  invisible(x)
}

#' Print Method for ABRM Objects
#' 
#' Prints a concise summary of a fitted ABRM, including the number of
#' parameters estimated and, when true parameter values are available, the
#' mean absolute bias and 95\% credible-interval coverage rate.
#' 
#' @param x An abrm object
#' @param ... Additional arguments (unused)
#' 
#' @return Invisibly returns \code{x}. Called for its side effect of printing
#'   a brief summary of the simulated misaligned spatial data.
#'
#' @examples
#' \donttest{
#'   ## Step 1: Simulate misaligned spatial data
#'   sim_data <- simulate_misaligned_data(
#'     seed              = 1,
#'     dist_covariates_x = c("normal", "poisson"),
#'     dist_covariates_y = c("normal", "poisson"),
#'     dist_y            = "poisson",
#'     x_intercepts      = c(0, -1),
#'     y_intercepts      = c(0, -1),
#'     beta0_y           = -1,
#'     beta_x            = c(0.1, -0.05),
#'     beta_y            = c(-0.1, 0.05)
#'   )
#'
#'   ## Step 2: Retrieve the NIMBLE model code
#'   model_code <- get_abrm_model()
#'
#'   ## Step 3: Fit the ABRM
#'   ## (niter and nburnin are intentionally small for illustration only;
#'   ##  use larger values, e.g. niter = 50000, nburnin = 30000, in practice)
#'   results <- run_abrm(
#'     gridx       = sim_data$gridx,
#'     gridy       = sim_data$gridy,
#'     atoms       = sim_data$atoms,
#'     model_code  = model_code,
#'     true_params = sim_data$true_params,
#'     norm_idx_x  = 1,
#'     pois_idx_x  = 2,
#'     norm_idx_y  = 1,
#'     pois_idx_y  = 2,
#'     dist_y      = 2,
#'     niter       = 1000,
#'     nburnin     = 500,
#'     nchains     = 2,
#'     seed        = 1,
#'     save_plots  = FALSE
#'   )
#'
#'   ## Step 4: Print a concise model summary
#'   print(results)
#'   # Reports: number of parameters, mean absolute bias, coverage rate
#' }
#'
#' @export
print.abrm <- function(x, ...) {
  pe <- x$parameter_estimates
  cat("\nABRM Model Results\n==================\n")
  cat("Parameters estimated:", nrow(pe), "\n")
  
  # Show bias and coverage if true values were available
  if ("bias" %in% names(pe))
    cat("Mean absolute bias: ",
        round(mean(abs(pe$bias), na.rm = TRUE), 4), "\n")
  if ("within_ci" %in% names(pe))
    cat("Coverage rate:      ",
        round(mean(pe$within_ci, na.rm = TRUE) * 100, 2), "%\n")
  
  cat("\nUse summary() for the full parameter table.\n")
  invisible(x)
}

#' Plot Method for ABRM Objects
#'
#' Displays MCMC diagnostic plots — trace plots and posterior density plots —
#' for the parameters monitored during model fitting. If no diagnostic plots
#' are stored in the object, a message is issued instead.
#'
#' @param x   An object of class \code{"abrm"} returned by
#'   \code{\link{run_abrm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}. Called for its side effect of rendering
#'   the MCMC diagnostic plots.
#'
#' @examples
#' \donttest{
#'   ## Step 1: Simulate misaligned spatial data
#'   sim_data <- simulate_misaligned_data(
#'     seed              = 1,
#'     dist_covariates_x = c("normal", "poisson"),
#'     dist_covariates_y = c("normal", "poisson"),
#'     dist_y            = "poisson",
#'     x_intercepts      = c(0, -1),
#'     y_intercepts      = c(0, -1),
#'     beta0_y           = -1,
#'     beta_x            = c(0.1, -0.05),
#'     beta_y            = c(-0.1, 0.05)
#'   )
#'
#'   ## Step 2: Retrieve the NIMBLE model code
#'   model_code <- get_abrm_model()
#'
#'   ## Step 3: Fit the ABRM
#'   results <- run_abrm(
#'     gridx       = sim_data$gridx,
#'     gridy       = sim_data$gridy,
#'     atoms       = sim_data$atoms,
#'     model_code  = model_code,
#'     true_params = sim_data$true_params,
#'     norm_idx_x  = 1,
#'     pois_idx_x  = 2,
#'     norm_idx_y  = 1,
#'     pois_idx_y  = 2,
#'     dist_y      = 2,
#'     niter       = 1000,
#'     nburnin     = 500,
#'     nchains     = 2,
#'     seed        = 1,
#'     save_plots  = FALSE
#'   )
#'
#'   ## Step 4: Display MCMC trace and posterior density plots
#'   plot(results)
#' }
#'
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

#' Variance-Covariance Method for ABRM Objects
#' Extracts variance-covariance matrices for regression coefficients from 
#' MCMC posterior samples. Returns separate matrices for X-grid and Y-grid
#' coefficients.
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
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
#' \donttest{
#'   ## Step 1: Simulate misaligned spatial data
#'   sim_data <- simulate_misaligned_data(
#'     seed              = 1,
#'     dist_covariates_x = c("normal", "poisson"),
#'     dist_covariates_y = c("normal", "poisson"),
#'     dist_y            = "poisson",
#'     x_intercepts      = c(0, -1),
#'     y_intercepts      = c(0, -1),
#'     beta0_y           = -1,
#'     beta_x            = c(0.1, -0.05),
#'     beta_y            = c(-0.1, 0.05)
#'   )
#'
#'   ## Step 2: Retrieve the NIMBLE model code
#'   model_code <- get_abrm_model()
#'
#'   ## Step 3: Fit the ABRM
#'   results <- run_abrm(
#'     gridx       = sim_data$gridx,
#'     gridy       = sim_data$gridy,
#'     atoms       = sim_data$atoms,
#'     model_code  = model_code,
#'     true_params = sim_data$true_params,
#'     norm_idx_x  = 1,
#'     pois_idx_x  = 2,
#'     norm_idx_y  = 1,
#'     pois_idx_y  = 2,
#'     dist_y      = 2,
#'     niter       = 1000,
#'     nburnin     = 500,
#'     nchains     = 2,
#'     seed        = 1,
#'     save_plots  = FALSE
#'   )
#'
#'   ## Step 4: Extract posterior variance-covariance matrices
#'   vcov_mats <- vcov(results)
#'
#'   # Posterior standard errors for X-grid and Y-grid coefficients
#'   sqrt(diag(vcov_mats$vcov_beta_x))
#'   sqrt(diag(vcov_mats$vcov_beta_y))
#'
#'   # Posterior correlation matrix across all beta parameters
#'   cov2cor(vcov_mats$vcov_all)
#' }
#'
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

#' Print Method for vcov.abrm Objects
#'
#' Prints all variance-covariance matrices stored in a \code{"vcov.abrm"}
#' object returned by \code{\link{vcov.abrm}}.
#'
#' @param x   An object of class \code{"vcov.abrm"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}. Called for its side effect of printing
#'   the intercept variance and the X-grid and Y-grid covariance matrices.
#'
#' @examples
#' \donttest{
#'   ## Step 1: Simulate misaligned spatial data
#'   sim_data <- simulate_misaligned_data(
#'     seed              = 1,
#'     dist_covariates_x = c("normal", "poisson"),
#'     dist_covariates_y = c("normal", "poisson"),
#'     dist_y            = "poisson",
#'     x_intercepts      = c(0, -1),
#'     y_intercepts      = c(0, -1),
#'     beta0_y           = -1,
#'     beta_x            = c(0.1, -0.05),
#'     beta_y            = c(-0.1, 0.05)
#'   )
#'
#'   ## Step 2: Retrieve the NIMBLE model code
#'   model_code <- get_abrm_model()
#'
#'   ## Step 3: Fit the ABRM
#'   results <- run_abrm(
#'     gridx       = sim_data$gridx,
#'     gridy       = sim_data$gridy,
#'     atoms       = sim_data$atoms,
#'     model_code  = model_code,
#'     true_params = sim_data$true_params,
#'     norm_idx_x  = 1,
#'     pois_idx_x  = 2,
#'     norm_idx_y  = 1,
#'     pois_idx_y  = 2,
#'     dist_y      = 2,
#'     niter       = 1000,
#'     nburnin     = 500,
#'     nchains     = 2,
#'     seed        = 1,
#'     save_plots  = FALSE
#'   )
#'
#'   ## Step 4: Extract and print the variance-covariance matrices
#'   vcov_mats <- vcov(results)
#'   print(vcov_mats)
#'   # Prints: intercept variance, X-grid and Y-grid covariance matrices
#' }
#'
#' @export
print.vcov.abrm <- function(x, ...) {
  cat("\nVariance-Covariance Matrices for ABRM Model\n")
  cat("============================================\n\n")
  
  if (!is.null(x$vcov_beta_0)) {
    cat("Intercept Variance:\n")
    cat("  ", names(x$vcov_beta_0), ": ", round(x$vcov_beta_0, 6), "\n", sep = "")
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

# ───────────────── New S3 methods for abrm ────────────────────────────────

#' Extract Coefficients from an ABRM Object
#'
#' Returns the posterior mean estimates for all regression coefficients as a
#' named numeric vector.
#'
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector of posterior mean estimates. Names correspond
#'   to parameter labels (e.g., \code{"beta_x[1]"}, \code{"beta_y[1]"},
#'   \code{"beta_0_y"}).
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "normal", dist_covariates_y = "normal",
#'     dist_y = "normal", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), norm_idx_x = 1, norm_idx_y = 1,
#'     dist_y = 1, niter = 1000, nburnin = 500, nchains = 2, seed = 1,
#'     save_plots = FALSE
#'   )
#'   coef(results)
#' }
#'
#' @export
coef.abrm <- function(object, ...) {
  pe   <- object$all_parameters
  keep <- grep("^beta_0_y$|^beta_x\\[|^beta_y\\[", pe$variable)
  pe   <- pe[keep, , drop = FALSE]
  structure(pe$estimated_beta, names = pe$variable)
}

#' Credible Intervals for an ABRM Object
#'
#' Returns the 95\% posterior credible intervals for all estimated parameters.
#' Note: these are Bayesian credible intervals from MCMC quantiles, not
#' frequentist confidence intervals.
#'
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param parm  Optional character vector of parameter names to subset. If
#'   \code{NULL} (default), all parameters are returned.
#' @param level Ignored. Intervals are always the 95\% credible intervals stored
#'   during model fitting. Included only for S3 compatibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix with two columns, \code{CI.Lower} and \code{CI.Upper}, and
#'   one row per parameter.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "normal", dist_covariates_y = "normal",
#'     dist_y = "normal", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), norm_idx_x = 1, norm_idx_y = 1,
#'     dist_y = 1, niter = 1000, nburnin = 500, nchains = 2, seed = 1,
#'     save_plots = FALSE
#'   )
#'   confint(results)
#'   confint(results, parm = "beta_x[1]")
#' }
#'
#' @importFrom stats confint quantile
#' @export
confint.abrm <- function(object, parm = NULL, level = 0.95, ...) {
  # Fast path: if level == 0.95, use the stored CI columns directly
  if (isTRUE(all.equal(level, 0.95))) {
    pe  <- object$parameter_estimates
    mat <- matrix(c(pe$ci_lower, pe$ci_upper), ncol = 2,
                  dimnames = list(pe$variable, c("CI.Lower", "CI.Upper")))
    if (!is.null(parm)) mat <- mat[rownames(mat) %in% parm, , drop = FALSE]
    return(mat)
  }
  
  # General path: recompute from raw MCMC samples (mirrors vcov.abrm logic)
  samples <- object$mcmc_results$samples
  if (inherits(samples, "mcmc.list")) {
    smat <- do.call(rbind, samples)
  } else if (is.list(samples) && all(sapply(samples, is.matrix))) {
    smat <- do.call(rbind, samples)
  } else {
    smat <- as.matrix(samples)
  }
  
  if (!is.null(parm)) smat <- smat[, colnames(smat) %in% parm, drop = FALSE]
  
  alpha <- 1 - level
  lo <- apply(smat, 2, quantile, probs = alpha / 2)
  hi <- apply(smat, 2, quantile, probs = 1 - alpha / 2)
  mat <- cbind(CI.Lower = lo, CI.Upper = hi)
  mat
}

#' Fitted Values for an ABRM Object
#'
#' Computes posterior-mean fitted values at the Y-grid level by applying the
#' estimated regression coefficients to the supplied covariate data. Because
#' the \code{abrm} object does not store the original grids, \code{gridy} and
#' \code{gridx} must be supplied by the user.
#'
#' @param object  An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param gridy   The Y-grid \code{sf} data frame used to fit the model.
#'   Must contain columns named \code{covariate_y_1}, \code{covariate_y_2}, etc.
#' @param gridx   The X-grid \code{sf} data frame used to fit the model.
#'   Must contain columns named \code{covariate_x_1}, \code{covariate_x_2}, etc.
#' @param atoms   The atoms \code{sf} data frame containing \code{ID_y} and
#'   \code{ID_x} linking atoms to their parent Y- and X-grid cells, plus a
#'   \code{population} column for offset weighting.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector of posterior-mean fitted values (on the linear predictor scale), one per Y-grid cell.
#' 
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "normal", dist_covariates_y = "normal",
#'     dist_y = "normal", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), norm_idx_x = 1, norm_idx_y = 1,
#'     dist_y = 1, niter = 1000, nburnin = 500, nchains = 2, seed = 1,
#'     save_plots = FALSE
#'   )
#'   fitted(results,
#'          gridy = sim_data$gridy,
#'          gridx = sim_data$gridx,
#'          atoms = sim_data$atoms)
#' }
#'
#' @importFrom stats coef fitted
#' @export
fitted.abrm <- function(object, gridy = NULL, gridx = NULL, atoms = NULL, ...) {
  
  # Guard: all three grids must be supplied
  if (is.null(gridy) || is.null(gridx) || is.null(atoms)) {
    stop(
      "fitted.abrm requires the original grid data.\n",
      "Please supply: fitted(object, gridy = gridy, gridx = gridx, atoms = atoms)\n",
      "Hint: if you used sim_data, pass sim_data$gridy, sim_data$gridx, sim_data$atoms."
    )
  }
  
  # Extract posterior mean coefficients
  betas     <- coef(object)                               # named vector
  beta_0    <- betas[grep("^beta_0_y$|^beta0_y$", names(betas))]
  beta_x    <- betas[grep("^beta_x\\[",            names(betas))]
  beta_y    <- betas[grep("^beta_y\\[",            names(betas))]
  
  # Drop geometry for arithmetic
  gridy_df  <- sf::st_drop_geometry(gridy)
  gridx_df  <- sf::st_drop_geometry(gridx)
  atoms_df  <- sf::st_drop_geometry(atoms)
  
  # Identify covariate columns
  x_cov_cols <- sort(grep("^covariate_x_", names(gridx_df), value = TRUE))
  y_cov_cols <- sort(grep("^covariate_y_", names(gridy_df), value = TRUE))
  
  if (length(beta_x) != length(x_cov_cols))
    warning("Number of beta_x coefficients (", length(beta_x),
            ") does not match X-grid covariate columns (", length(x_cov_cols), ").")
  if (length(beta_y) != length(y_cov_cols))
    warning("Number of beta_y coefficients (", length(beta_y),
            ") does not match Y-grid covariate columns (", length(y_cov_cols), ").")
  
  # Aggregate X-grid covariates to Y-grid via atoms
  # For each Y-grid cell, compute population-weighted mean of X covariates
  # across all atoms belonging to that cell.
  y_ids <- gridy_df$ID
  
  x_contrib <- matrix(0, nrow = length(y_ids), ncol = length(beta_x))
  
  for (i in seq_along(y_ids)) {
    atom_idx  <- which(atoms_df$ID_y == y_ids[i])
    if (length(atom_idx) == 0) next
    
    atom_x_ids  <- atoms_df$ID_x[atom_idx]
    atom_pops   <- atoms_df$population[atom_idx]
    total_pop   <- sum(atom_pops)
    if (total_pop == 0) next
    
    for (k in seq_along(x_cov_cols)) {
      # Get the X-grid value for each atom's parent X-cell, weighted by population
      x_vals <- gridx_df[[x_cov_cols[k]]][match(atom_x_ids, gridx_df$ID)]
      x_contrib[i, k] <- sum(x_vals * atom_pops, na.rm = TRUE) / total_pop
    }
  }
  
  # Y-grid covariate matrix
  y_contrib <- as.matrix(gridy_df[, y_cov_cols, drop = FALSE])
  
  # Compute linear predictor: beta_0 + X*beta_x + Z*beta_y
  linear_pred <- as.numeric(beta_0) +
    as.numeric(x_contrib %*% beta_x) +
    as.numeric(y_contrib %*% beta_y)
  
  # Return named vector
  cell_names  <- if (!is.null(rownames(gridy_df))) rownames(gridy_df) else paste0("y_cell_", y_ids)
  structure(linear_pred, names = cell_names)
}

#' Generic Function for WAIC
#'
#' Extracts the Widely Applicable Information Criterion (WAIC) from a
#' fitted model object.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A \code{"waic_abrm"} object for \code{abrm} models. See
#'   \code{\link{waic.abrm}} for details.
#'
#' @export
waic <- function(object, ...) UseMethod("waic")

#' WAIC for an ABRM Object
#'
#' Returns the Widely Applicable Information Criterion (WAIC) for a fitted
#' ABRM. WAIC is the recommended information criterion for Bayesian models and
#' is computed by NIMBLE when \code{compute_waic = TRUE} is passed to
#' \code{\link{run_abrm}}.
#'
#' @param object An object of class \code{"abrm"} returned by
#'   \code{\link{run_abrm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"waic_abrm"}, which is a named list with
#'   components:
#'   \describe{
#'     \item{\code{waic}}{The scalar WAIC value. Lower values indicate better
#'       predictive accuracy penalised for model complexity.}
#'     \item{\code{lppd}}{Log pointwise predictive density — the
#'       goodness-of-fit component of WAIC.}
#'     \item{\code{penalty}}{The effective number of parameters
#'       (\eqn{p_{\mathrm{WAIC}}}) — the complexity penalty component.}
#'     \item{\code{n_params}}{The number of regression parameters estimated.}
#'   }
#'   Returns an error if the model was not fitted with
#'   \code{compute_waic = TRUE}.
#'
#' @details
#' WAIC (Watanabe 2010) is defined as:
#' \deqn{\mathrm{WAIC} = -2 \cdot \mathrm{lppd} + 2 \cdot p_{\mathrm{WAIC}}}
#' where \eqn{\mathrm{lppd}} is the log pointwise predictive density and
#' \eqn{p_{\mathrm{WAIC}}} is the effective number of parameters. Unlike AIC,
#' WAIC is fully Bayesian and averages over the posterior distribution rather
#' than evaluating at a point estimate.
#'
#' WAIC values are only meaningful when compared between models fitted on
#' \emph{identical} data with the same outcome distribution. A true marginal
#' log-likelihood is not available for ABRM, so \code{\link[stats]{AIC}} and
#' \code{\link[stats]{logLik}} are not supported; use \code{waic()} instead.
#'
#' @references
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' \emph{Journal of Machine Learning Research}, 11, 3571--3594.
#'
#' @seealso \code{\link{waic}}, \code{\link{print.waic_abrm}},
#'   \code{\link{run_abrm}}
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed              = 1,
#'     dist_covariates_x = "normal",
#'     dist_covariates_y = "normal",
#'     dist_y            = "normal",
#'     x_intercepts      = 0,
#'     y_intercepts      = 0,
#'     beta0_y           = 0,
#'     beta_x            = 0.1,
#'     beta_y            = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx        = sim_data$gridx,
#'     gridy        = sim_data$gridy,
#'     atoms        = sim_data$atoms,
#'     model_code   = get_abrm_model(),
#'     norm_idx_x   = 1,
#'     norm_idx_y   = 1,
#'     dist_y       = 1,
#'     niter        = 1000,
#'     nburnin      = 500,
#'     nchains      = 2,
#'     seed         = 1,
#'     compute_waic = TRUE    # required; re-fit if this was FALSE
#'   )
#'
#'   w <- waic(results)
#'   print(w)          # formatted table: WAIC, lppd, pWAIC, n_params
#'   w$waic            # raw scalar for comparing models
#'   w$penalty         # effective number of parameters (pWAIC)
#' }
#'
#' @export
waic.abrm <- function(object, ...) {
  
  raw <- object$mcmc_results$WAIC
  
  if (is.null(raw)) {
    stop(
      "WAIC is not available: the model was fitted without compute_waic = TRUE.\n",
      "Re-fit using: run_abrm(..., compute_waic = TRUE)"
    )
  }
  
  # NIMBLE returns a list; field names vary slightly by version.
  # Defensively extract with fallbacks.
  waic_val  <- raw$WAIC
  lppd_val  <- if (!is.null(raw$lppd)) raw$lppd else NA_real_
  pen_val   <- if (!is.null(raw$pWAIC)) raw$pWAIC else
    if (!is.null(raw$penalty)) raw$penalty else NA_real_
  
  result <- list(
    waic     = waic_val,
    lppd     = lppd_val,
    penalty  = pen_val,
    n_params = nrow(object$parameter_estimates)
  )
  class(result) <- "waic_abrm"
  result
}

#' Print Method for waic_abrm Objects
#'
#' Displays a formatted WAIC summary for a fitted ABRM.
#'
#' @param x   An object of class \code{"waic_abrm"} returned by
#'   \code{\link{waic.abrm}}.
#' @param digits Integer; number of decimal places. Default \code{3}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.waic_abrm <- function(x, digits = 3, ...) {
  cat("WAIC for ABRM\n")
  cat("=============\n")
  cat("WAIC   :", round(x$waic, digits), "\n")
  if (!is.na(x$lppd))
    cat("lppd   :", round(x$lppd, digits),
        " (log pointwise predictive density)\n")
  if (!is.na(x$penalty))
    cat("pWAIC  :", round(x$penalty, digits),
        " (effective number of parameters)\n")
  cat("n_params:", x$n_params,
      " (regression coefficients estimated)\n")
  cat("\nNote: lower WAIC indicates better predictive fit.\n",
      "      Compare only models fitted on identical data.\n", sep = "")
  invisible(x)
}