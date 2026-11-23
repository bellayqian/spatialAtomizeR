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
#' Calculates convergence diagnostics (Gelman-Rubin R-hat and effective
#' sample size) for MCMC output. Parameter names are used as they appear in
#' the samples (\code{run_nimble_model()} has already renamed them).
#'
#' @param mcmc_output Output from NIMBLE MCMC containing \code{samples}, a
#'   list of per-chain matrices.
#' @param rhat_threshold,ess_threshold Thresholds used for the
#'   \code{converged} flag.
#' @return List with \code{diagnostics} (data frame: parameter, mean, sd,
#'   rhat, ess, converged), \code{group_convergence}, \code{overall_converged},
#'   \code{n_chains} and \code{chain_correlation}.
#' @importFrom coda mcmc mcmc.list gelman.diag effectiveSize
#' @importFrom stats cor sd median
#' @keywords internal
#' @noRd
check_mcmc_diagnostics <- function(mcmc_output, rhat_threshold = 1.1,
                                   ess_threshold = 400) {
  chains_list <- mcmc_output$samples
  if (is.matrix(chains_list)) chains_list <- list(chain1 = chains_list)
  n_chains <- length(chains_list)

  mcmc_chains <- coda::mcmc.list(lapply(chains_list, coda::mcmc))
  pooled <- do.call(rbind, chains_list)
  param_names <- colnames(pooled)

  if (n_chains >= 2L) {
    gelman <- tryCatch(
      coda::gelman.diag(mcmc_chains, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1],
      error = function(e) rep(NA_real_, length(param_names)))
  } else {
    warning("R-hat requires at least 2 chains; only effective sample sizes are reported.",
            call. = FALSE)
    gelman <- rep(NA_real_, length(param_names))
  }
  eff_size <- tryCatch(coda::effectiveSize(mcmc_chains),
                       error = function(e) rep(NA_real_, length(param_names)))

  diagnostics <- data.frame(
    parameter = param_names,
    mean = colMeans(pooled),
    sd = apply(pooled, 2, stats::sd),
    rhat = unname(gelman),
    ess = unname(eff_size),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  diagnostics$converged <- ifelse(is.na(diagnostics$rhat), NA,
                                  diagnostics$rhat < rhat_threshold &
                                    diagnostics$ess > ess_threshold)

  param_groups <- list(
    beta = grep("^beta", param_names, value = TRUE),
    tau = grep("^tau", param_names, value = TRUE),
    prec = grep("^Prec", param_names, value = TRUE)
  )
  group_convergence <- vapply(param_groups, function(group_params) {
    all(diagnostics$converged[diagnostics$parameter %in% group_params], na.rm = FALSE)
  }, logical(1))

  list(
    diagnostics = diagnostics,
    group_convergence = group_convergence,
    overall_converged = all(diagnostics$converged),
    n_chains = n_chains,
    chain_correlation = stats::cor(pooled)
  )
}

#' Reshape a list of MCMC chains to long format for ggplot2
#' @param chains_list List of per-chain matrices (iterations x parameters).
#' @param params Columns to keep (default: all).
#' @param labels Optional facet labels, same length as \code{params}.
#' @keywords internal
#' @noRd
.chains_long <- function(chains_list, params = NULL, labels = NULL) {
  if (is.null(params)) params <- colnames(chains_list[[1]])
  if (is.null(labels)) labels <- params
  do.call(rbind, lapply(seq_along(chains_list), function(ch) {
    m <- chains_list[[ch]][, params, drop = FALSE]
    data.frame(
      iteration = rep(seq_len(nrow(m)), times = ncol(m)),
      chain = factor(ch, levels = seq_along(chains_list)),
      variable = factor(rep(labels, each = nrow(m)), levels = unique(labels)),
      value = as.vector(m)
    )
  }))
}

#' Create a Trace Plot from MCMC Chains
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap theme_minimal labs theme element_text
#' @keywords internal
#' @noRd
create_trace_plot <- function(chains_list, params = NULL, labels = NULL,
                              title = "MCMC trace plots") {
  plot_data <- .chains_long(chains_list, params, labels)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = value, color = chain)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::facet_wrap(~variable, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Iteration (after burn-in and thinning)",
                  y = "Parameter value", color = "Chain") +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 11))
}

#' Create a Posterior Density Plot from MCMC Chains
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap theme_minimal labs theme element_text
#' @keywords internal
#' @noRd
create_density_plot <- function(chains_list, params = NULL, labels = NULL,
                                title = "Posterior density plots") {
  plot_data <- .chains_long(chains_list, params, labels)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = value, fill = chain)) +
    ggplot2::geom_density(alpha = 0.3) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Parameter value", y = "Density", fill = "Chain") +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 11))
}

#' Print Convergence Summary
#'
#' @param convergence_results Results from check_mcmc_diagnostics
#' @return No return value, called for side effects (messages with median ESS,
#'   maximum R-hat, and parameters with high relative variance).
#' @keywords internal
#' @noRd
print_convergence_summary <- function(convergence_results) {
  diagnostics <- convergence_results$diagnostics
  if (is.null(diagnostics) || nrow(diagnostics) == 0L) return(invisible(NULL))
  msg <- sprintf("Diagnostic summary: median ESS %.1f", median(diagnostics$ess, na.rm = TRUE))
  if (any(!is.na(diagnostics$rhat))) {
    msg <- sprintf("%s, max R-hat %.3f", msg, max(diagnostics$rhat, na.rm = TRUE))
  }
  message(msg)

  high_var <- diagnostics[diagnostics$sd / abs(diagnostics$mean) > 0.5, , drop = FALSE]
  if (nrow(high_var) > 0L) {
    message("Parameters with high relative variance (sd/|mean| > 0.5):\n",
            paste(utils::capture.output(
              print(high_var[, c("parameter", "mean", "sd")], row.names = FALSE)),
              collapse = "\n"))
  }
  invisible(NULL)
}

#' Create Summary Statistics
#'
#' @param all_results Combined results data frame
#' @param output_dir Output directory
#' @return A list with two components: \code{method_summary} (data frame containing overall summary statistics for each method, including mean bias, relative bias, and coverage rates) and \code{param_summary} (data frame containing parameter-specific comparisons across methods). The function also saves these summaries to CSV files in the output directory.
#' @importFrom dplyr group_by summarize
#' @importFrom tidyr pivot_wider
#' @importFrom utils write.csv
#' @keywords internal
#' @noRd
create_summary_statistics <- function(all_results, output_dir) {
  summary_stats <- all_results |>
    dplyr::group_by(method) |>
    dplyr::summarize(
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(relative_bias)),
      rmse = sqrt(mean(bias^2)),
      coverage_rate = mean(within_ci) * 100,
      .groups = "drop"
    )

  utils::write.csv(summary_stats, file = file.path(output_dir, "method_summary_stats.csv"),
                   row.names = FALSE)

  param_summary <- all_results |>
    dplyr::group_by(method, variable) |>
    dplyr::summarize(
      mean_estimated_beta = mean(estimated_beta),
      mean_true_beta = mean(true_beta),
      mean_bias = mean(bias),
      mean_relative_bias = mean(relative_bias),
      coverage_rate = mean(within_ci) * 100,
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = variable,
      names_from = method,
      values_from = c(mean_estimated_beta, mean_bias, mean_relative_bias, coverage_rate),
      names_sep = "_"
    )

  utils::write.csv(param_summary, file = file.path(output_dir, "parameter_summary.csv"),
                   row.names = FALSE)

  message("Method Comparison Summary:")
  print(summary_stats)
  message("Parameter-Specific Comparison:")
  print(param_summary)

  list(method_summary = summary_stats, param_summary = param_summary)
}
