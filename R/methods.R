#' Print Method for Misaligned Data Objects
#'
#' @param x A misaligned_data object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the simulated misaligned spatial data including grid dimensions and number of atoms.
#' @examples
#' ## Simulate misaligned spatial data
#' sim_data <- simulate_misaligned_data(
#'   seed              = 1,
#'   dist_covariates_x = c("norm", "pois"),
#'   dist_covariates_y = c("norm", "pois"),
#'   dist_y            = "pois",
#'   x_intercepts      = c(0, -1),
#'   y_intercepts      = c(0, -1),
#'   beta0_y           = -1,
#'   beta_x            = c(0.1, -0.05),
#'   beta_y            = c(-0.1, 0.05)
#' )
#' print(sim_data)
#' @export
print.misaligned_data <- function(x, ...) {
  # Count only the covariate columns themselves, not the per-covariate
  # population columns (e.g. population_covariate_x_1) that can sit beside them.
  n_x <- length(grep("^covariate_x_[0-9]+$", names(x$gridx)))
  n_y <- length(grep("^covariate_y_[0-9]+$", names(x$gridy)))
  message("Simulated Misaligned Spatial Data")
  message("==================================")
  message("Y-grid cells: ", nrow(x$gridy))
  message("X-grid cells: ", nrow(x$gridx))
  message("Atoms: ", if (is.null(x$atoms)) "not computed" else nrow(x$atoms))
  message("Number of X covariates: ", n_x)
  message("Number of Y covariates: ", n_y)
  if (!is.null(x$true_params)) {
    message("True parameters available")
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
#' sim_data <- simulate_misaligned_data(
#'   seed = 1, dist_covariates_x = c("norm", "pois"),
#'   dist_covariates_y = c("norm", "pois"), dist_y = "pois",
#'   x_intercepts = c(0, -1), y_intercepts = c(0, -1),
#'   beta0_y = -1, beta_x = c(0.1, -0.05), beta_y = c(-0.1, 0.05)
#' )
#' summary(sim_data)
#'
#' @export
summary.misaligned_data <- function(object, ...) {
  print.misaligned_data(object)
  message("True beta_x: ", paste(object$true_params$beta_x, collapse = ", "))
  message("True beta_y: ", paste(object$true_params$beta_y, collapse = ", "))
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
#'   \code{which = "atoms"} is an error if \code{x$atoms} has not been
#'   computed (see \code{\link{create_atoms}}).
#' @param color_y Color for the Y-grid boundaries. Default \code{"blue"}.
#' @param color_x Color for the X-grid boundaries. Default \code{"orange"}.
#' @param title Optional character string for the plot title. If \code{NULL},
#'   a default title is used.
#' @param ... Additional arguments passed to \code{ggplot2::geom_sf()} when
#'   \code{which} is not \code{"both"}.
#'
#' @return A \code{ggplot} object. It is drawn when printed (for example at
#'   the console), and can be modified with further \code{ggplot2} layers or
#'   saved with \code{ggplot2::ggsave()}.
#'
#' @examples
#' sim_data <- simulate_misaligned_data(
#'   seed = 1, dist_covariates_x = c("norm", "pois"),
#'   dist_covariates_y = c("norm", "pois"), dist_y = "pois",
#'   x_intercepts = c(0, -1), y_intercepts = c(0, -1),
#'   beta0_y = -1, beta_x = c(0.1, -0.05), beta_y = c(-0.1, 0.05)
#' )
#' # Default: overlay both grids to visualise misalignment
#' plot(sim_data)
#'
#' # Plot a single component
#' plot(sim_data, which = "gridy")
#' plot(sim_data, which = "atoms")
#'
#' # The result is a ggplot object that can be extended
#' p <- plot(sim_data, color_y = "steelblue", color_x = "firebrick")
#' p + ggplot2::labs(caption = "Simulated data")
#'
#' @importFrom ggplot2 ggplot geom_sf labs theme_void
#' @export
plot.misaligned_data <- function(x, which = c("both", "gridy", "gridx", "atoms"),
                                 color_y = "blue", color_x = "orange", title = NULL, ...) {
  which <- match.arg(which)

  layer_ok <- function(obj) inherits(obj, "sf") && nrow(obj) > 0L
  if (which == "atoms" && !layer_ok(x$atoms)) {
    stop("The atoms have not been computed for this object. Create them with ",
         "x$atoms <- create_atoms(x$gridy, x$gridx) before plotting which = \"atoms\".",
         call. = FALSE)
  }
  needed <- if (which == "both") c("gridx", "gridy") else which
  for (comp in needed) {
    if (!layer_ok(x[[comp]])) {
      stop(sprintf("x$%s is missing or is not a non-empty sf data frame.", comp),
           call. = FALSE)
    }
  }

  if (which == "both") {
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = x[["gridx"]], fill = NA, color = color_x, linewidth = 1.2) +
      ggplot2::geom_sf(data = x[["gridy"]], fill = NA, color = color_y,
                       linetype = "dashed", linewidth = 0.5) +
      ggplot2::labs(
        title    = if (!is.null(title)) title else "Spatial Misalignment",
        subtitle = sprintf("%s: gridx  |  %s (dashed): gridy", color_x, color_y)
      ) +
      ggplot2::theme_void()
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_sf(
        data  = x[[which]],
        fill  = NA,
        color = switch(which, gridy = color_y, gridx = color_x, "black"),
        ...
      ) +
      ggplot2::labs(title = if (!is.null(title)) title else paste("Spatial layout:", which)) +
      ggplot2::theme_void()
  }
}


#' Summary Method for ABRM Objects
#'
#' Returns a structured summary of a fitted ABRM, including posterior means,
#' standard deviations, and 95% credible intervals for all regression
#' coefficients. When true parameter values are available (i.e., the model
#' was fitted with \code{true_params}), mean absolute bias and
#' credible-interval coverage are also reported.
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
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", true_params = sim_data$true_params,
#'     niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   summary(results)    # full parameter table with bias and coverage
#' }
#'
#' @export
summary.abrm <- function(object, ...) {
  pe <- object$parameter_estimates
  retlist <- list(
    n_params      = nrow(pe),
    has_true      = "true_beta" %in% names(pe),
    mean_abs_bias = if ("bias" %in% names(pe))
      round(mean(abs(pe$bias), na.rm = TRUE), 4) else NA,
    coverage      = if ("within_ci" %in% names(pe))
      round(mean(pe$within_ci, na.rm = TRUE) * 100, 2) else NA,
    estimates     = pe
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
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   s <- summary(results)
#'   print(s, digits = 2)
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
#' Prints a concise summary of a fitted ABRM: the number of parameters
#' estimated, a one-line convergence check and, when true parameter values
#' are available, the mean absolute bias and 95\% credible-interval coverage.
#'
#' @param x An abrm object
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = c("norm", "pois"),
#'     dist_covariates_y = c("norm", "pois"), dist_y = "pois",
#'     x_intercepts = c(0, -1), y_intercepts = c(0, -1),
#'     beta0_y = -1, beta_x = c(0.1, -0.05), beta_y = c(-0.1, 0.05)
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(),
#'     x_types = c("norm", "pois"), y_types = c("norm", "pois"), dist_y = "pois",
#'     true_params = sim_data$true_params,
#'     niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   print(results)
#' }
#'
#' @export
print.abrm <- function(x, ...) {
  pe <- x$parameter_estimates
  cat("\nABRM Model Results\n==================\n")
  cat("Parameters estimated:", nrow(pe), "\n")

  d <- x$mcmc_results$convergence$diagnostics
  if (!is.null(d) && any(!is.na(d$converged))) {
    cat(sprintf("Converged (R-hat < 1.1, ESS > 400): %d of %d parameters\n",
                sum(d$converged, na.rm = TRUE), nrow(d)))
  }
  if ("bias" %in% names(pe))
    cat("Mean absolute bias: ", round(mean(abs(pe$bias), na.rm = TRUE), 4), "\n")
  if ("within_ci" %in% names(pe))
    cat("Coverage rate:      ", round(mean(pe$within_ci, na.rm = TRUE) * 100, 2), "%\n")

  cat("\nUse summary() for the full parameter table and mcmc_diagnostics() for\n",
      "convergence details.\n", sep = "")
  invisible(x)
}

# ----------------------------------------------------------------------------
# MCMC plots
# ----------------------------------------------------------------------------

#' MCMC Trace and Density Plots for ABRM Objects
#'
#' \code{plot_trace()} draws MCMC trace plots and \code{plot_density()} draws
#' posterior density plots of the regression coefficients of a fitted ABRM,
#' one panel per parameter and one colour per chain. Panels are labelled with
#' the covariate names. \code{plot()} on an \code{abrm} object calls one of
#' these, chosen by \code{type}.
#'
#' @param x An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param parm Optional character vector of parameters to plot, given either
#'   as covariate names (e.g. \code{"intercept"}, \code{"covariate_x_1"}) or
#'   as internal names (e.g. \code{"beta_x[1]"}). Default: all coefficients.
#' @param title Plot title.
#' @param type For \code{plot()}: \code{"trace"} (default) or \code{"density"}.
#' @param ... For \code{plot()}: further arguments passed to
#'   \code{plot_trace()} or \code{plot_density()}.
#'
#' @return A \code{ggplot} object, drawn when printed. It can be modified
#'   with additional \code{ggplot2} layers or saved with
#'   \code{ggplot2::ggsave()}.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   plot_trace(results)
#'   plot_density(results, parm = c("covariate_x_1", "covariate_y_1"))
#'   plot(results, type = "density")   # same as plot_density(results)
#' }
#'
#' @name abrm_plots
NULL

#' @rdname abrm_plots
#' @export
plot_trace <- function(x, parm = NULL, title = "MCMC trace plots") {
  args <- .abrm_plot_args(x, parm)
  create_trace_plot(args$chains, params = args$raw, labels = args$labels, title = title)
}

#' @rdname abrm_plots
#' @export
plot_density <- function(x, parm = NULL, title = "Posterior density plots") {
  args <- .abrm_plot_args(x, parm)
  create_density_plot(args$chains, params = args$raw, labels = args$labels, title = title)
}

#' @rdname abrm_plots
#' @export
plot.abrm <- function(x, type = c("trace", "density"), ...) {
  type <- match.arg(type)
  switch(type,
         trace   = plot_trace(x, ...),
         density = plot_density(x, ...))
}

#' Collect chains, parameter names and labels for the plotting functions
#' @keywords internal
#' @noRd
.abrm_plot_args <- function(x, parm) {
  if (!inherits(x, "abrm")) stop("'x' must be an abrm object from run_abrm().", call. = FALSE)
  chains <- x$mcmc_results$samples
  if (is.null(chains)) stop("No MCMC samples are stored in this object.", call. = FALSE)
  if (is.matrix(chains)) chains <- list(chain1 = chains)
  raw <- .abrm_beta_names(x, colnames(chains[[1]]), parm)
  if (length(raw) == 0L) stop("No parameters to plot.", call. = FALSE)
  list(chains = chains, raw = raw, labels = .abrm_label(x, raw))
}

# ----------------------------------------------------------------------------
# Convergence diagnostics accessor
# ----------------------------------------------------------------------------

#' MCMC Convergence Diagnostics for ABRM Objects
#'
#' Returns the convergence diagnostics computed after fitting: posterior
#' mean and standard deviation, Gelman-Rubin potential scale reduction factor
#' (R-hat) and effective sample size (ESS) for each regression coefficient.
#'
#' @param object An object of class \code{"abrm"} returned by
#'   \code{\link{run_abrm}}.
#' @param parm Optional character vector of parameters, given as covariate
#'   names or internal names (as for \code{\link{plot_trace}}). Default: all.
#' @param rhat_threshold,ess_threshold A parameter is flagged as converged
#'   when \code{rhat < rhat_threshold} and \code{ess > ess_threshold}.
#'   Defaults 1.1 and 400.
#'
#' @return A data frame with one row per parameter and columns
#'   \code{variable} (covariate name or \code{"intercept"}), \code{parameter}
#'   (internal name), \code{mean}, \code{sd}, \code{rhat}, \code{ess} and
#'   \code{converged}. \code{rhat} (and therefore \code{converged}) is
#'   \code{NA} when only one chain was run. Use
#'   \code{all(mcmc_diagnostics(fit)$converged)} for an overall check.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   mcmc_diagnostics(results)
#'   mcmc_diagnostics(results, ess_threshold = 100)
#'   all(mcmc_diagnostics(results)$converged)
#' }
#'
#' @export
mcmc_diagnostics <- function(object, parm = NULL, rhat_threshold = 1.1,
                             ess_threshold = 400) {
  if (!inherits(object, "abrm")) {
    stop("'object' must be an abrm object from run_abrm().", call. = FALSE)
  }
  d <- object$mcmc_results$convergence$diagnostics
  if (is.null(d)) {
    stop("No convergence diagnostics are stored in this object.", call. = FALSE)
  }
  betas <- .abrm_beta_names(object, d$parameter, parm)
  rows <- if (is.null(parm)) {
    c(match(betas, d$parameter), setdiff(seq_len(nrow(d)), match(betas, d$parameter)))
  } else {
    match(betas, d$parameter)
  }
  d <- d[rows, , drop = FALSE]
  out <- data.frame(
    variable  = .abrm_label(object, d$parameter),
    parameter = d$parameter,
    mean      = d$mean,
    sd        = d$sd,
    rhat      = d$rhat,
    ess       = d$ess,
    converged = ifelse(is.na(d$rhat), NA, d$rhat < rhat_threshold & d$ess > ess_threshold),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out
}

# ----------------------------------------------------------------------------
# vcov / coef / confint
# ----------------------------------------------------------------------------

#' Variance-Covariance Method for ABRM Objects
#'
#' Extracts variance-covariance matrices for regression coefficients from
#' MCMC posterior samples. Returns separate matrices for X-grid, Y-grid and
#' (if present) atom-level coefficients, plus the combined matrix. Rows and
#' columns are labelled with covariate names.
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param ... Additional arguments (unused)
#' @return A list with class "vcov.abrm" containing:
#'   \item{vcov_beta_0}{Posterior variance of the intercept (named scalar)}
#'   \item{vcov_beta_x}{Variance-covariance matrix for X-grid coefficients}
#'   \item{vcov_beta_y}{Variance-covariance matrix for Y-grid coefficients}
#'   \item{vcov_beta_a}{Variance-covariance matrix for atom-level
#'     coefficients (only when atom-level covariates were used)}
#'   \item{vcov_all}{Combined variance-covariance matrix for all coefficients}
#' @details
#' The matrices are computed from the posterior samples pooled across chains.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = c("norm", "pois"),
#'     dist_covariates_y = c("norm", "pois"), dist_y = "pois",
#'     x_intercepts = c(0, -1), y_intercepts = c(0, -1),
#'     beta0_y = -1, beta_x = c(0.1, -0.05), beta_y = c(-0.1, 0.05)
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(),
#'     x_types = c("norm", "pois"), y_types = c("norm", "pois"), dist_y = "pois",
#'     niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   vcov_mats <- vcov(results)
#'   sqrt(diag(vcov_mats$vcov_beta_x))
#'   cov2cor(vcov_mats$vcov_all)
#' }
#'
#' @importFrom stats var cov
#' @export
vcov.abrm <- function(object, ...) {
  smat <- .abrm_sample_matrix(object)
  if (is.null(colnames(smat))) {
    stop("MCMC samples do not have column names. Cannot identify parameters.", call. = FALSE)
  }
  raw <- .abrm_beta_names(object, colnames(smat))
  if (length(raw) == 0L) stop("No regression coefficients found in the MCMC samples.",
                              call. = FALSE)

  block <- function(pattern) {
    cols <- raw[grepl(pattern, raw)]
    if (length(cols) == 0L) return(NULL)
    m <- stats::cov(smat[, cols, drop = FALSE])
    dimnames(m) <- list(.abrm_label(object, cols), .abrm_label(object, cols))
    m
  }

  vcov_list <- list()
  b0 <- raw[raw == "beta_0_y"]
  if (length(b0)) {
    vcov_list$vcov_beta_0 <- stats::setNames(stats::var(smat[, b0]), .abrm_label(object, b0))
  }
  vcov_list$vcov_beta_x <- block("^beta_x\\[")
  vcov_list$vcov_beta_y <- block("^beta_y\\[")
  vcov_list$vcov_beta_a <- block("^beta_a\\[")
  vcov_list$vcov_all    <- block(".")

  class(vcov_list) <- c("vcov.abrm", "list")
  vcov_list
}

#' Print Method for vcov.abrm Objects
#'
#' Prints all variance-covariance matrices stored in a \code{"vcov.abrm"}
#' object returned by \code{\link{vcov.abrm}}.
#'
#' @param x   An object of class \code{"vcov.abrm"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.vcov.abrm <- function(x, ...) {
  cat("\nVariance-Covariance Matrices for ABRM Model\n")
  cat("============================================\n\n")

  if (!is.null(x$vcov_beta_0)) {
    cat("Intercept Variance:\n")
    cat("  ", names(x$vcov_beta_0), ": ", round(x$vcov_beta_0, 6), "\n\n", sep = "")
  }
  show_block <- function(m, heading) {
    if (is.null(m)) return(invisible())
    cat(heading, "\n", sep = "")
    cat("  Dimensions:", nrow(m), "x", ncol(m), "\n")
    print(round(m, 6))
    cat("\n")
  }
  show_block(x$vcov_beta_x, "X-Grid Coefficients Variance-Covariance Matrix:")
  show_block(x$vcov_beta_y, "Y-Grid Coefficients Variance-Covariance Matrix:")
  show_block(x$vcov_beta_a, "Atom-Level Coefficients Variance-Covariance Matrix:")

  cat("Use vcov_object$vcov_beta_x, vcov_object$vcov_beta_y, etc. to access matrices\n")
  invisible(x)
}

#' Extract Coefficients from ABRM Objects
#'
#' Returns the posterior mean estimates for all regression coefficients as a
#' named numeric vector.
#'
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector of posterior mean estimates, named
#'   \code{"intercept"} followed by the covariate names.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   coef(results)
#' }
#'
#' @export
coef.abrm <- function(object, ...) {
  pe <- object$parameter_estimates
  stats::setNames(pe$estimated_beta, pe$variable)
}

#' Credible Intervals for ABRM Objects
#'
#' Returns equal-tailed posterior credible intervals for the regression
#' coefficients, computed from the MCMC samples pooled across chains.
#' Note: these are Bayesian credible intervals from MCMC quantiles, not
#' frequentist confidence intervals.
#'
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param parm Optional character vector of parameters, given as covariate
#'   names (e.g. \code{"covariate_x_1"}, \code{"intercept"}) or internal
#'   names (e.g. \code{"beta_x[1]"}). If \code{NULL} (default), all
#'   coefficients are returned.
#' @param level Credible level (default \code{0.95}).
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix with two columns, \code{CI.Lower} and \code{CI.Upper}, and
#'   one row per parameter, labelled with covariate names.
#'
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   confint(results)
#'   confint(results, parm = "covariate_x_1", level = 0.9)
#' }
#'
#' @importFrom stats confint quantile
#' @export
confint.abrm <- function(object, parm = NULL, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)
  }
  smat <- .abrm_sample_matrix(object)
  raw <- .abrm_beta_names(object, colnames(smat), parm)
  alpha <- 1 - level
  lo <- apply(smat[, raw, drop = FALSE], 2, stats::quantile, probs = alpha / 2)
  hi <- apply(smat[, raw, drop = FALSE], 2, stats::quantile, probs = 1 - alpha / 2)
  matrix(c(lo, hi), ncol = 2,
         dimnames = list(.abrm_label(object, raw), c("CI.Lower", "CI.Upper")))
}

#' Fitted Values for ABRM Objects
#'
#' Computes posterior-mean fitted values at the Y-grid level by applying the
#' estimated regression coefficients to the supplied covariate data. 
#' @param object An object of class \code{"abrm"} returned by \code{\link{run_abrm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per Y-grid cell and four columns:
#'   \describe{
#'     \item{\code{y_grid_id}}{The Y-grid cell ID, as supplied in \code{gridy}.}
#'     \item{\code{observed}}{The observed outcome value \code{y} for each
#'       Y-grid cell. Note: the outcome is observed at the Y-grid level for
#'       all cells. This is distinct from covariates, which may be
#'       \emph{latent} (unobserved at the Y-grid level) for cells that span
#'       multiple atoms.}
#'     \item{\code{fitted}}{The model-predicted outcome on the original
#'       outcome scale (counts for Poisson/binomial; continuous for normal).}
#'     \item{\code{residual}}{Observed minus fitted.}
#'   }
#' 
#' @examples
#' \donttest{
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", x_intercepts = 0, y_intercepts = 0,
#'     beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   results <- run_abrm(
#'     gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
#'     model_code = get_abrm_model(), x_types = "norm", y_types = "norm",
#'     dist_y = "norm", niter = 1000, nburnin = 500, nchains = 2, seed = 1
#'   )
#'   fitted(results)
#' }
#'
#' @importFrom stats coef fitted
#' @export
fitted.abrm <- function(object, ...) {
  if (is.null(object$fitted_values)) {
    stop("No fitted values found. Re-fit the model with the current version of run_abrm().")
  }
  
  observed <- object$y_observed
  fitted_vals <- object$fitted_values
  
  data.frame(
    y_grid_id = object$y_grid_ids,
    observed  = observed,
    fitted    = round(fitted_vals, 4),
    residual  = round(observed - fitted_vals, 4),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
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

#' WAIC for ABRM Objects
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
#'     dist_covariates_x = "norm",
#'     dist_covariates_y = "norm",
#'     dist_y            = "norm",
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
#'     x_types      = "norm",
#'     y_types      = "norm",
#'     dist_y       = "norm",
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
  cat("\nNote: lower WAIC indicates better predictive fit.\n",
      "      Compare only models fitted on identical data.\n", sep = "")
  invisible(x)
}