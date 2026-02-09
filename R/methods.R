#' Print method for abrm objects
#' @param x An abrm object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the ABRM model results including convergence status, number of parameters estimated, and key fit statistics.
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
#' @importFrom graphics par
plot.abrm <- function(x, ...) {
  if(!is.null(x$mcmc_results$convergence$plots)) {
    print(x$mcmc_results$convergence$plots$trace)
    print(x$mcmc_results$convergence$plots$density)
  } else {
    message("No diagnostic plots available.\n")
  }
  invisible(x)
}

#' Print method for abrm_comparison objects
#' @param x A abrm_comparison object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of method comparison results including the number of simulations and methods compared.
#' @keywords internal
#' @noRd
print.abrm_comparison <- function(x, ...) {
  message("Method Comparison Results\n")
  message("=========================\n\n")
  methods <- unique(x$combined_comparison$method)
  for (m in methods) {
    message(m, "method:\n")
    subset_data <- x$combined_comparison[x$combined_comparison$method == m, ]
    message("  Mean absolute bias:", round(mean(abs(subset_data$bias)), 4), "\n")
    message("  RMSE:", round(sqrt(mean(subset_data$bias^2)), 4), "\n")
  }
  message("\nUse summary() for detailed comparison\n")
  invisible(x)
}

#' Summary method for abrm_comparison objects
#' 
#' @param object An object of class "abrm_comparison"
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object \code{object}. The function is called for its side effect of printing the method comparison summary including combined comparison results across all simulations.
#' @keywords internal
#' @noRd
summary.abrm_comparison <- function(object, ...) {
  message("Method Comparison Summary\n")
  message("=========================\n\n")
  
  # Calculate summary by method
  for(method in unique(object$combined_comparison$method)) {
    method_data <- object$combined_comparison[object$combined_comparison$method == method, ]
    message(sprintf("\n%s Method:\n", method))
    message(sprintf("  Mean Absolute Bias: %.4f\n", mean(abs(method_data$bias))))
    message(sprintf("  RMSE: %.4f\n", sqrt(mean(method_data$bias^2))))
    message(sprintf("  Coverage Rate: %.1f%%\n", mean(method_data$within_ci) * 100))
  }
  
  invisible(object)
}

#' Plot method for abrm_comparison objects
#' 
#' @param x An object of class "abrm_comparison"
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of displaying comparison plots between methods if true parameters are available in the comparison object.
#' @keywords internal
#' @noRd
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
    message("Cannot create plots without true parameters.\n")
  }
  invisible(x)
}

#' Print method for sensitivity_analysis objects
#' @param x A sensitivity_analysis object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the sensitivity analysis results including the number of simulations and scenarios tested.
#' @keywords internal
#' @noRd
print.sensitivity_analysis <- function(x, ...) {
  message("Sensitivity Analysis Results\n")
  message("============================\n\n")
  message("Number of simulations:", nrow(x$combined_results) / 
        length(unique(x$combined_results$variable)), "\n")
  message("Correlation values tested:", 
      unique(x$combined_results$x_correlation), "\n")
  message("Output directory:", x$output_dir, "\n\n")
  message("Use summary() for detailed statistics\n")
  invisible(x)
}

#' Summary method for sensitivity_analysis objects
#' @param object A sensitivity_analysis object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{object}. The function is called for its side effect of printing the sensitivity analysis summary including comprehensive results across all simulation scenarios.
#' @keywords internal
#' @noRd
summary.sensitivity_analysis <- function(object, ...) {
  message("Sensitivity Analysis Summary\n")
  message("============================\n\n")
  message(object$summary_by_correlation)
  invisible(object)
}

#' Print method for misaligned_data objects
#' @param x A misaligned_data object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object \code{x}. The function is called for its side effect of printing a summary of the simulated misaligned spatial data including grid dimensions and number of atoms.
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
summary.misaligned_data <- function(object, ...) {
  print.misaligned_data(object)
  message("True beta_x: ", paste(object$true_params$beta_x, collapse = ", "), "\n")
  message("True beta_y: ", paste(object$true_params$beta_y, collapse = ", "), "\n")
  invisible(object)
}