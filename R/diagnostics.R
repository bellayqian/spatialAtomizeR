# Declare global variables to avoid R CMD check notes
utils::globalVariables(c(
  "iteration", "value", "chain", "variable", "method", "estimated_beta",
  "ci_lower", "ci_upper", "true_beta", "sim_number", "relative_bias",
  "within_ci", "covered", "bias", "x_correlation", "y_correlation",
  "mean_estimate", "mean_lower", "mean_upper", "true_value", "mean_rel_bias",
  "coverage_rate", "converged", "mean_estimated_beta", "mean_bias", "mean_relative_bias"
))

#' Check MCMC Diagnostics
#'
#' Calculates convergence diagnostics including Gelman-Rubin statistics
#' and effective sample sizes for MCMC output
#'
#' @param mcmc_output Output from NIMBLE MCMC containing samples
#' @param sim_metadata Optional list with simulation metadata
#' @return List with diagnostic results including convergence statistics and plots
#' @importFrom coda mcmc mcmc.list gelman.diag effectiveSize
#' @importFrom stats cor sd median
#' @importFrom utils head
#' @keywords internal
#' @noRd
check_mcmc_diagnostics <- function(mcmc_output, sim_metadata = NULL) {
  # Validate chain structure
  chains_list <- mcmc_output$samples
  if (!is.list(chains_list) || length(chains_list) < 2) {
    warning("Insufficient number of chains for convergence diagnostics")
    return(list(converged = FALSE, error = "insufficient_chains"))
  }
  
  # Convert to mcmc objects
  chains_mcmc <- lapply(chains_list, function(x) {
    if (is.null(colnames(x))) colnames(x) <- paste0("param", 1:ncol(x))
    mcmc(x)
  })
  
  mcmc_chains <- mcmc.list(chains_mcmc)
  
  # Calculate Gelman-Rubin statistics
  gelman_diag <- try(gelman.diag(mcmc_chains, autoburnin = FALSE, multivariate = FALSE)$psrf, silent = TRUE)
  
  # Calculate effective sample sizes
  eff_size <- try(effectiveSize(mcmc_chains), silent = TRUE)
  
  # Parameter means and SDs
  param_names <- colnames(chains_list[[1]])
  param_means <- colMeans(do.call(rbind, chains_list))
  param_sds <- apply(do.call(rbind, chains_list), 2, sd)
  
  # Create diagnostics data frame
  diagnostics <- data.frame(
    parameter = param_names,
    mean = param_means,
    sd = param_sds,
    rhat = if(!inherits(gelman_diag, "try-error")) gelman_diag[,1] else rep(NA, length(param_names)),
    ess = if(!inherits(eff_size, "try-error")) eff_size else rep(NA, length(param_names)),
    stringsAsFactors = FALSE
  )
  
  # Add convergence assessment
  diagnostics$converged <- with(diagnostics, rhat < 1.1 & ess > 400)
  
  # Group parameters
  param_groups <- list(
    beta = grep("^beta", param_names, value = TRUE),
    tau = grep("^tau", param_names, value = TRUE),
    prec = grep("^Prec", param_names, value = TRUE)
  )
  
  # Check group convergence
  group_convergence <- sapply(param_groups, function(group_params) {
    group_diag <- diagnostics[diagnostics$parameter %in% group_params,]
    all(group_diag$converged)
  })
  
  # Create diagnostic plots
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if(is.null(sim_metadata)) {
      # Create default metadata if none provided
      sim_metadata <- list(
        sim_number = 1,
        x_correlation = NA,
        y_correlation = NA
      )
    }
    
    # Generate plots using updated function
    plots <- create_diagnostic_plots(chains_list, sim_metadata)
  } else {
    plots <- NULL
  }
  
  return(list(
    diagnostics = diagnostics,
    group_convergence = group_convergence,
    overall_converged = all(diagnostics$converged),
    plots = plots,
    chain_correlation = cor(do.call(rbind, chains_list))
  ))
}

#' Create Diagnostic Plots
#'
#' @param chains_list List of MCMC chains
#' @param sim_metadata Simulation metadata
#'
#' @return List with trace and density plots
#' @importFrom ggplot2 ggplot aes geom_line geom_density facet_wrap theme_minimal labs theme element_text
#' @importFrom reshape2 melt
#' @keywords internal
#' @noRd
create_diagnostic_plots <- function(chains_list, sim_metadata) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting but is not installed.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required for plotting but is not installed.")
  }
  
  # Get parameters from chain data
  params <- colnames(chains_list[[1]])
  
  # Create title with metadata
  plot_title <- sprintf("MCMC Diagnostics (x_cor=%.1f, y_cor=%.1f, sim=%d)",
                        sim_metadata$x_correlation, 
                        sim_metadata$y_correlation, 
                        sim_metadata$sim_number)
  
  # Prepare data for plotting
  plot_data <- lapply(seq_along(chains_list), function(chain) {
    data <- as.data.frame(chains_list[[chain]][, params, drop = FALSE])
    data$iteration <- 1:nrow(data)
    data$chain <- as.factor(chain)
    reshape2::melt(data, id.vars = c("iteration", "chain"))
  })
  
  # Combine and ensure it's a proper data frame
  plot_data <- as.data.frame(do.call(rbind, plot_data))
  
  # Ensure variable is a factor
  plot_data$variable <- as.factor(plot_data$variable)
  
  # Create trace plots - use .data pronoun for NSE variables
  trace_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = value, color = chain)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::facet_wrap(~variable, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = plot_title,
                  subtitle = "Trace Plots for All Parameters",
                  x = "Iteration", y = "Parameter Value") +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 11))
  
  # Create density plots
  density_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = chain)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = plot_title,
                  subtitle = "Density Plots for All Parameters",
                  x = "Parameter Value", y = "Density") +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 11))
  
  return(list(trace = trace_plot, density = density_plot))
}

#' Print Convergence Summary
#'
#' @param convergence_results Results from check_mcmc_diagnostics
#' @return No return value, called for side effects. Prints convergence diagnostics including overall convergence status, minimum effective sample size, median ESS, maximum Rhat statistic, and identifies parameters with high relative variance.
#' @keywords internal
#' @noRd
print_convergence_summary <- function(convergence_results) {
  message("Convergence Summary:\n")
  message("===================\n")
  
  # Overall convergence
  message(sprintf("Overall convergence: %s\n", 
              ifelse(convergence_results$overall_converged, "Achieved", "Not achieved")))
  
  # Group convergence
  message("Parameter Group Convergence:\n")
  for(group in names(convergence_results$group_convergence)) {
    message(sprintf("%s parameters: %s\n", 
                group,
                ifelse(convergence_results$group_convergence[[group]], 
                       "Converged", "Not converged")))
  }
  
  # Print detailed summary for problematic parameters
  problem_params <- subset(convergence_results$diagnostics, 
                           !converged)
  
  if(nrow(problem_params) > 0) {
    message("Parameters requiring attention:\n")
    print(problem_params[, c("parameter", "rhat", "ess")])
  }
  
  # Print general diagnostics summary
  message("Diagnostic Summary Statistics:\n")
  diagnostics <- convergence_results$diagnostics
  message(sprintf("Median ESS: %.1f\n", median(diagnostics$ess, na.rm = TRUE)))
  message(sprintf("Max Rhat: %.3f\n", max(diagnostics$rhat, na.rm = TRUE)))
  
  # Check for any unusually high variances
  high_var_params <- subset(diagnostics, sd/abs(mean) > 0.5)
  if(nrow(high_var_params) > 0) {
    message("Parameters with high relative variance:\n")
    print(high_var_params[, c("parameter", "mean", "sd")])
  }
}

#' Create Comparison Plots
#'
#' @param comparison_data Data frame with comparison results
#' @param output_dir Output directory for plots
#' @param true_params True parameter values (optional)
#' @return No return value, called for side effects. Creates and displays/saves comparison plots of coefficient estimates, bias, and coverage rates across different methods. If output_dir is provided, plots are saved to a PDF file.
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_bar geom_hline coord_flip theme_minimal labs position_dodge facet_wrap .data
#' @importFrom dplyr select mutate
#' @importFrom grDevices pdf dev.off
#' @keywords internal
#' @noRd
create_comparison_plots <- function(comparison_data, output_dir, true_params = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting but is not installed.")
  }
  
  # Ensure we have the method column
  if(!"method" %in% names(comparison_data)) {
    warning("No 'method' column found in comparison_data")
    comparison_data$method <- "ABRM"
  }
  
  # Pivot the data for side by side comparison
  plot_data <- comparison_data %>%
    dplyr::select(.data$sim_number, .data$method, .data$variable, .data$estimated_beta, 
                  .data$ci_lower, .data$ci_upper, .data$true_beta) %>%
    dplyr::mutate(variable = factor(.data$variable, levels = unique(.data$variable)))
  
  # Create coefficient plot
  coef_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable, y = .data$estimated_beta, color = .data$method)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper), 
                           position = ggplot2::position_dodge(width = 0.5), width = 0.2) +
    ggplot2::geom_point(ggplot2::aes(y = .data$true_beta), color = "black", shape = 4, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Comparison of Coefficient Estimates by Method",
                  subtitle = "Points = Estimated values with 95% CI\nX = True values",
                  x = "Variable",
                  y = "Coefficient Value")
  
  # Create bias plot
  bias_data <- comparison_data %>%
    dplyr::select(.data$method, .data$variable, .data$relative_bias) %>%
    dplyr::mutate(variable = factor(.data$variable, levels = unique(.data$variable)))
  
  bias_plot <- ggplot2::ggplot(bias_data, ggplot2::aes(x = .data$variable, y = .data$relative_bias, fill = .data$method)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Relative Bias by Method (%)",
                  x = "Variable",
                  y = "Relative Bias (%)")
  
  # Create coverage plot (within_ci)
  coverage_data <- comparison_data %>%
    dplyr::select(.data$method, .data$variable, .data$within_ci) %>%
    dplyr::mutate(variable = factor(.data$variable, levels = unique(.data$variable)),
                  covered = ifelse(.data$within_ci, "Yes", "No"))
  
  coverage_plot <- ggplot2::ggplot(coverage_data, ggplot2::aes(x = .data$variable, fill = .data$covered)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::facet_wrap(~method) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "95% CI Coverage by Method",
                  subtitle = "Whether true parameter value falls within estimated 95% CI",
                  x = "Variable",
                  y = "Count")
  
  # Save plots
  if(!is.null(output_dir)) {
    grDevices::pdf(file.path(output_dir, "coefficient_comparison.pdf"), width = 10, height = 8)
    print(coef_plot)
    print(bias_plot)
    print(coverage_plot)
    grDevices::dev.off()
    message("Comparison plots saved to ", file.path(output_dir, "coefficient_comparison.pdf"))
  } else {
    # Just display plots without saving
    print(coef_plot)
    print(bias_plot)
    print(coverage_plot)
  }
}

#' Create Summary Statistics
#'
#' @param all_results Combined results data frame
#' @param output_dir Output directory
#' @return A list with two components: \code{method_summary} (data frame containing overall summary statistics for each method, including mean bias, relative bias, and coverage rates) and \code{param_summary} (data frame containing parameter-specific comparisons across methods). The function also saves these summaries to CSV files in the output directory.
#' @importFrom dplyr group_by summarize select
#' @importFrom tidyr pivot_wider
#' @importFrom utils write.csv
#' @keywords internal
#' @noRd
create_summary_statistics <- function(all_results, output_dir) {
  # Calculate summary statistics by method
  summary_stats <- all_results %>%
    group_by(method) %>%
    summarize(
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(relative_bias)),
      rmse = sqrt(mean(bias^2)),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Save summary
  write.csv(summary_stats, file = file.path(output_dir, "method_summary_stats.csv"), 
            row.names = FALSE)
  
  # Create a summary table for each parameter
  # FIXED: Aggregate first, then pivot
  param_summary <- all_results %>%
    group_by(method, variable) %>%
    summarize(
      mean_estimated_beta = mean(estimated_beta),
      mean_true_beta = mean(true_beta),
      mean_bias = mean(bias),
      mean_relative_bias = mean(relative_bias),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    ) %>%
    pivot_wider(
      id_cols = variable,
      names_from = method,
      values_from = c(mean_estimated_beta, mean_bias, mean_relative_bias, coverage_rate),
      names_sep = "_"
    )
  
  # Save parameter summary
  write.csv(param_summary, file = file.path(output_dir, "parameter_summary.csv"), 
            row.names = FALSE)
  
  # Print summary to console
  message("Method Comparison Summary:\n")
  print(summary_stats)
  
  message("Parameter-Specific Comparison:\n")
  print(param_summary)
  
  return(list(method_summary = summary_stats, param_summary = param_summary))
}

#' Create Sensitivity Summary Plots
#'
#' @param combined_results Combined results data frame
#' @param output_dir Output directory
#' @return No return value, called for side effects. Creates and displays/saves summary plots for sensitivity analysis including correlation effects, bias patterns, and coverage rates. If output_dir is provided, plots are saved to a PDF file.
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_bar geom_hline coord_flip facet_grid theme_minimal labs position_dodge
#' @importFrom grDevices pdf dev.off
#' @importFrom dplyr group_by summarize
#' @keywords internal
#' @noRd
create_sensitivity_summary_plots <- function(combined_results, output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting but is not installed.")
  }
  
  # Add correlation values to the combined results
  plot_data <- combined_results %>%
    dplyr::group_by(.data$method, .data$variable, .data$x_correlation, .data$y_correlation) %>%
    dplyr::summarize(
      mean_estimate = mean(.data$estimated_beta),
      mean_lower = mean(.data$ci_lower),
      mean_upper = mean(.data$ci_upper),
      true_value = mean(.data$true_beta),
      mean_bias = mean(.data$bias),
      mean_rel_bias = mean(.data$relative_bias),
      coverage_rate = mean(.data$within_ci) * 100,
      .groups = 'drop'
    )
  
  # Create faceted plot by correlation values
  corr_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable, y = .data$mean_estimate, color = .data$method)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_lower, ymax = .data$mean_upper), 
                           position = ggplot2::position_dodge(width = 0.5), width = 0.2) +
    ggplot2::geom_point(ggplot2::aes(y = .data$true_value), color = "black", shape = 4, size = 3) +
    ggplot2::facet_grid(x_correlation ~ y_correlation, labeller = ggplot2::label_both) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Parameter Estimates by Correlation Structure",
                  subtitle = "Averaged across all simulations",
                  x = "Variable",
                  y = "Coefficient Value")
  
  # Create bias plot by correlation
  bias_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable, y = .data$mean_rel_bias, fill = .data$method)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(x_correlation ~ y_correlation, labeller = ggplot2::label_both) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Relative Bias by Correlation Structure (%)",
                  subtitle = "Averaged across all simulations",
                  x = "Variable",
                  y = "Relative Bias (%)")
  
  # Create coverage rate plot
  coverage_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable, y = .data$coverage_rate, fill = .data$method)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    ggplot2::facet_grid(x_correlation ~ y_correlation, labeller = ggplot2::label_both) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "95% CI Coverage Rate by Correlation Structure",
                  subtitle = "Percentage of simulations where true value falls within CI",
                  x = "Variable", 
                  y = "Coverage Rate (%)")
  
  # Save plots
  if(!is.null(output_dir)) {
    grDevices::pdf(file.path(output_dir, "sensitivity_summary_plots.pdf"), width = 12, height = 10)
    print(corr_plot)
    print(bias_plot)
    print(coverage_plot)
    grDevices::dev.off()
    message("Sensitivity summary plots saved to ", file.path(output_dir, "sensitivity_summary_plots.pdf"))
  } else {
    # Just display plots without saving
    print(corr_plot)
    print(bias_plot)
    print(coverage_plot)
  }
}

