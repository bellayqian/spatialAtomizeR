#' Print method for abrm objects
#' @param x An abrm object
#' @param ... Additional arguments (unused)
#' @export
print.abrm <- function(x, ...) {
  cat("ABRM Model Results\n")
  cat("==================\n\n")
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
#' @export
summary.abrm <- function(object, ...) {
  cat("ABRM Model Summary\n")
  cat("==================\n\n")
  print(object$parameter_estimates)
  invisible(object)
}

#' Plot method for abrm objects
#' 
#' @param x An object of class "abrm"
#' @param ... Additional arguments (ignored)
#' @export
#' @importFrom graphics par
plot.abrm <- function(x, ...) {
  if(!is.null(x$mcmc_results$convergence$plots)) {
    print(x$mcmc_results$convergence$plots$trace)
    print(x$mcmc_results$convergence$plots$density)
  } else {
    cat("No diagnostic plots available.\n")
  }
  invisible(x)
}

#' Print method for abrm_comparison objects
#' @param x A abrm_comparison object
#' @param ... Additional arguments (unused)
#' @export
print.abrm_comparison <- function(x, ...) {
  cat("Method Comparison Results\n")
  cat("=========================\n\n")
  methods <- unique(x$combined_comparison$method)
  for (m in methods) {
    cat(m, "method:\n")
    subset_data <- x$combined_comparison[x$combined_comparison$method == m, ]
    cat("  Mean absolute bias:", round(mean(abs(subset_data$bias)), 4), "\n")
    cat("  RMSE:", round(sqrt(mean(subset_data$bias^2)), 4), "\n")
  }
  cat("\nUse summary() for detailed comparison\n")
  invisible(x)
}

#' Summary method for abrm_comparison objects
#' 
#' @param object An object of class "abrm_comparison"
#' @param ... Additional arguments (ignored)
#' @export
summary.abrm_comparison <- function(object, ...) {
  cat("Method Comparison Summary\n")
  cat("=========================\n\n")
  
  # Calculate summary by method
  for(method in unique(object$combined_comparison$method)) {
    method_data <- object$combined_comparison[object$combined_comparison$method == method, ]
    cat(sprintf("\n%s Method:\n", method))
    cat(sprintf("  Mean Absolute Bias: %.4f\n", mean(abs(method_data$bias))))
    cat(sprintf("  RMSE: %.4f\n", sqrt(mean(method_data$bias^2))))
    cat(sprintf("  Coverage Rate: %.1f%%\n", mean(method_data$within_ci) * 100))
  }
  
  invisible(object)
}

#' Plot method for abrm_comparison objects
#' 
#' @param x An object of class "abrm_comparison"
#' @param ... Additional arguments (ignored)
#' @export
plot.abrm_comparison <- function(x, ...) {
  if (!is.null(x$abrm_results$parameter_estimates) && 
      "true_beta" %in% names(x$abrm_results$parameter_estimates)) {
    true_params <- list(
      beta_x = x$abrm_results$parameter_estimates$true_beta[
        grep("covariate_x", x$abrm_results$parameter_estimates$variable)
      ],
      beta_y = x$abrm_results$parameter_estimates$true_beta[
        grep("covariate_y", x$abrm_results$parameter_estimates$variable)
      ]
    )
    create_comparison_plots(x$combined_comparison, tempdir(), true_params)
  } else {
    cat("Cannot create plots without true parameters.\n")
  }
  invisible(x)
}

#' Print method for sensitivity_analysis objects
#' @param x A sensitivity_analysis object
#' @param ... Additional arguments (unused)
#' @export
print.sensitivity_analysis <- function(x, ...) {
  cat("Sensitivity Analysis Results\n")
  cat("============================\n\n")
  cat("Number of simulations:", nrow(x$combined_results) / 
        length(unique(x$combined_results$variable)), "\n")
  cat("Correlation values tested:", 
      unique(x$combined_results$x_correlation), "\n")
  cat("Output directory:", x$output_dir, "\n\n")
  cat("Use summary() for detailed statistics\n")
  invisible(x)
}

#' Summary method for sensitivity_analysis objects
#' @param object A sensitivity_analysis object
#' @param ... Additional arguments (unused)
#' @export
summary.sensitivity_analysis <- function(object, ...) {
  cat("Sensitivity Analysis Summary\n")
  cat("============================\n\n")
  print(object$summary_by_correlation)
  invisible(object)
}

#' Print method for misaligned_data objects
#' @param x A misaligned_data object
#' @param ... Additional arguments (unused)
#' @export
print.misaligned_data <- function(x, ...) {
  cat("Simulated Misaligned Spatial Data\n")
  cat("==================================\n\n")
  cat("Y-grid cells:", nrow(x$gridy), "\n")
  cat("X-grid cells:", nrow(x$gridx), "\n")
  cat("Atoms:", nrow(x$atoms), "\n")
  cat("Number of X covariates:", length(grep("covariate_x", names(x$gridx))), "\n")
  cat("Number of Y covariates:", length(grep("covariate_y", names(x$gridy))), "\n")
  if (!is.null(x$true_params)) {
    cat("\nTrue parameters available\n")
  }
  invisible(x)
}

#' Summary method for misaligned_data objects
#' @param object A misaligned_data object
#' @param ... Additional arguments (unused)
#' @export
summary.misaligned_data <- function(object, ...) {
  print.misaligned_data(object)
  cat("\nTrue beta_x:", object$true_params$beta_x, "\n")
  cat("True beta_y:", object$true_params$beta_y, "\n")
  invisible(object)
}