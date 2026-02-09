#' Run ABRM Analysis
#'
#' Runs the Atom-Based Regression Model on simulated data
#'
#' @param gridx The X-grid sf dataframe, containing a numeric area ID variable named 'ID' and covariates named 'covariate_x_1','covariate_x_2',... 
#' @param gridy The Y-grid sf dataframe, containing a numeric area ID variable named 'ID', covariates named 'covariate_y_1','covariate_y_2',...., and an outcome named 'y'. 
#' @param atoms The atom sf dataframe, which should contain numeric variables named 'ID_x' and 'ID_y' holding the X-grid and Y-grid cell IDs for each atom, as well as an atom-level population count named 'population'.
#' @param model_code NIMBLE model code from get_abrm_model()
#' @param true_params The true outcome model regression coefficient parameters, if known (e.g., from simulate_misaligned_data())
#' @param norm_idx_x Vector of numeric indices of X-grid covariates (ordered as 'covariate_x_1','covariate_x_2',...) that should be treated as normally-distributed
#' @param pois_idx_x Vector of numeric indices of X-grid covariates (ordered as 'covariate_x_1','covariate_x_2',...) that should be treated as Poisson-distributed
#' @param binom_idx_x Vector of numeric indices of X-grid covariates (ordered as 'covariate_x_1','covariate_x_2',...) that should be treated as binomial-distributed
#' @param norm_idx_y Vector of numeric indices of Y-grid covariates (ordered as 'covariate_y_1','covariate_y_2',...) that should be treated as normally-distributed
#' @param pois_idx_y Vector of numeric indices of Y-grid covariates (ordered as 'covariate_y_1','covariate_y_2',...) that should be treated as Poisson-distributed
#' @param binom_idx_y Vector of numeric indices of Y-grid covariates (ordered as 'covariate_y_1','covariate_y_2',...) that should be treated as binomial-distributed
#' @param dist_y Distribution type for outcome (1=normal, 2=poisson, 3=binomial)
#' @param niter Number of MCMC iterations (default: 50000)
#' @param nburnin Number of burn-in iterations (default: 30000)
#' @param nchains Number of MCMC chains (default: 2)
#' @param thin Thinning interval (default: 10)
#' @param sim_metadata Optional simulation metadata list
#' @param save_plots Logical, whether to save diagnostic plots (default: TRUE)
#' @param output_dir Directory for saving outputs (default: NULL)
#'
#' @return List containing MCMC results and parameter estimates
#' @export
#' @importFrom dplyr %>% group_by summarize select mutate
#' @importFrom tidyr pivot_wider
run_abrm <- function(gridx,
                     gridy,
                     atoms,
                     model_code,
                     true_params=NULL,
                     norm_idx_x = NULL, 
                     pois_idx_x = NULL, 
                     binom_idx_x = NULL,
                     norm_idx_y = NULL, 
                     pois_idx_y = NULL, 
                     binom_idx_y = NULL,
                     dist_y = 2,
                     niter = 50000,
                     nburnin = 30000,
                     nchains = 2,
                     thin = 10,
                     sim_metadata = NULL,
                     save_plots = TRUE,
                     output_dir = NULL) {
  
  if (!('ID' %in% names(gridx))) stop("gridx must contain an area ID variable named 'ID'")
  
  if (!('ID' %in% names(gridy))) stop("gridy must contain an area ID variable named 'ID'")
  
  if (!('ID_y' %in% names(atoms) & 'ID_x' %in% names(atoms))) stop("atoms must contain an X-grid ID variable named 'ID_x' and a Y-grid ID variable named 'ID_y'")
  
  if (!('covariate_x_1' %in% names(gridx))) stop("gridx must contain at least 1 covariate named 'covariate_x_1'")
  
  if (!('covariate_y_1' %in% names(gridy))) stop("gridy must contain at least 1 covariate named 'covariate_y_1'") 
  
  if (!('y' %in% names(gridy))) stop("gridy must contain an outcome named 'y'")
  
  if (!('population' %in% names(atoms))) stop("atoms must contain a population count named 'population'")  
  
  sim_data<-list('gridx'=gridx,'gridy'=gridy,'atoms'=atoms,'true_params'=true_params)
  
  message("Preparing ABRM inputs...\n")
  
  # Prepare spatial bookkeeping
  bookkeeping <- prepare_spatial_bookkeeping(sim_data)
  
  # Prepare adjacency matrices
  adjacency <- prepare_adjacency_matrices(bookkeeping$gridy_yorder, bookkeeping$gridx_xorder)
  
  # Prepare NIMBLE inputs
  nimble_inputs <- prepare_nimble_inputs(
    bookkeeping, adjacency, sim_data,
    norm_idx_x = norm_idx_x,
    pois_idx_x = pois_idx_x,
    binom_idx_x = binom_idx_x,
    norm_idx_y = norm_idx_y,
    pois_idx_y = pois_idx_y,
    binom_idx_y = binom_idx_y,
    dist_y = dist_y
  )
  
  # Run NIMBLE model
  abrm_results <- run_nimble_model(
    constants = nimble_inputs$constants,
    data = nimble_inputs$data,
    inits = nimble_inputs$inits,
    sim_metadata = sim_metadata,
    model_code = model_code,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    thin = thin,
    save_plots = save_plots,
    output_dir = output_dir
  )
  
  # Extract parameter estimates
  abrm_parameters <- data.frame(
    variable = rownames(abrm_results$summary$all.chains),
    estimated_beta = abrm_results$summary$all.chains[, "Mean"],
    std_error = abrm_results$summary$all.chains[, "St.Dev."],
    ci_lower = abrm_results$summary$all.chains[, "95%CI_low"],
    ci_upper = abrm_results$summary$all.chains[, "95%CI_upp"],
    stringsAsFactors = FALSE
  )
  
  # If true parameters are available, calculate bias metrics
  if(!is.null(sim_data$true_params)) {
    true_beta_x <- sim_data$true_params$beta_x
    true_beta_y <- sim_data$true_params$beta_y
    true_values <- c(true_beta_x, true_beta_y)
    
    # Filter for beta_y parameters
    beta_params <- grep("^beta_y\\[", abrm_parameters$variable)
    abrm_betas <- abrm_parameters[beta_params, ]
    
    abrm_betas$true_beta <- true_values
    abrm_betas$variable <- c(
      paste0("covariate_x_", 1:length(true_beta_x)),
      paste0("covariate_y_", 1:length(true_beta_y))
    )
    
    # Calculate bias metrics
    abrm_betas$bias <- abrm_betas$estimated_beta - abrm_betas$true_beta
    abrm_betas$relative_bias <- ((abrm_betas$estimated_beta - abrm_betas$true_beta) / 
                                   abrm_betas$true_beta) * 100
    abrm_betas$within_ci <- (abrm_betas$true_beta >= abrm_betas$ci_lower & 
                               abrm_betas$true_beta <= abrm_betas$ci_upper)
    
    result <- list(
      mcmc_results = abrm_results,
      parameter_estimates = abrm_betas,
      all_parameters = abrm_parameters
    )
    class(result) <- "abrm"
    return(result)
  } else {
    result <- list(
      mcmc_results = abrm_results,
      parameter_estimates = abrm_parameters,
      all_parameters = abrm_parameters
    )
    class(result) <- "abrm"
    return(result)
  }
}

#' Run Both Methods and Compare
#'
#' Runs both ABRM and dasymetric mapping methods and compares results
#'
#' @param sim_data List of data elements to be used in the ABRM, structured like the output from the simulate_misaligned_data() function. The first element of this list is the Y-grid sf dataframe (named 'gridy'), containing a numeric area ID variable named 'ID_y', covariates named 'covariate_y_1','covariate_y_2',...., and an outcome named 'y'. The second element of this list is the X-grid sf dataframe (named 'gridx'), containing a numeric area ID variable named 'ID_x' and covariates named 'covariate_x_1','covariate_x_2',... The third element of the list is the atom sf dataframe (named 'atoms'), which should contain variables named 'ID_x' and 'ID_y' holding the X-grid and Y-grid cell IDs for each atom, as well as an atom-level population count named 'population'.
#' @param sim_metadata Simulation metadata
#' @param model_code NIMBLE model code
#' @param nimble_params List of NIMBLE parameters (niter, nburnin, thin, nchains)
#' @param output_dir Output directory
#' @param norm_idx_x Indices of normal X covariates
#' @param pois_idx_x Indices of Poisson X covariates
#' @param binom_idx_x Indices of binomial X covariates
#' @param norm_idx_y Indices of normal Y covariates
#' @param pois_idx_y Indices of Poisson Y covariates
#' @param binom_idx_y Indices of binomial Y covariates
#' @param dist_y Distribution type for outcome (1=normal, 2=poisson, 3=binomial)
#' @param outcome_type Outcome distribution name
#'
#' @return List with combined comparison, ABRM results, and dasymetric results
#' @export
#' @importFrom dplyr %>%
run_both_methods <- function(sim_data, sim_metadata, model_code, 
                             nimble_params, output_dir,
                             norm_idx_x, pois_idx_x, binom_idx_x,
                             norm_idx_y, pois_idx_y, binom_idx_y,
                             dist_y, outcome_type) {
  
  # Get true parameters
  true_beta_x <- sim_data$true_params$beta_x
  true_beta_y <- sim_data$true_params$beta_y
  
  # Run ABRM method
  message("  Running ABRM method...\n")
  abrm_full_results <- run_abrm(
    gridx = sim_data$gridx,
    gridy = sim_data$gridy,
    atoms = sim_data$atoms,
    model_code = model_code,
    true_params = sim_data$true_params,
    norm_idx_x = norm_idx_x,
    pois_idx_x = pois_idx_x,
    binom_idx_x = binom_idx_x,
    norm_idx_y = norm_idx_y,
    pois_idx_y = pois_idx_y,
    binom_idx_y = binom_idx_y,
    dist_y = dist_y,
    niter = nimble_params$niter,
    nburnin = nimble_params$nburnin,
    nchains = nimble_params$nchains,
    thin = nimble_params$thin,
    sim_metadata = sim_metadata,
    save_plots = TRUE,
    output_dir = output_dir
  )
  
  abrm_betas <- abrm_full_results$parameter_estimates
  abrm_betas$method <- "ABRM"
  
  # Run Dasymetric method
  message("  Running Dasymetric mapping method...\n")
  mapped_data <- dasymetric_mapping(sim_data)
  dasymetric_results <- fit_dasymetric_model(mapped_data, outcome_type)
  
  # Extract dasymetric results (exclude intercept)
  dasy_betas <- dasymetric_results[-1, ]  # Remove intercept row
  
  # Add metadata
  true_values <- c(true_beta_x, true_beta_y)
  dasy_betas$variable <- c(
    paste0("covariate_x_", 1:length(true_beta_x)),
    paste0("covariate_y_", 1:length(true_beta_y))
  )
  names(dasy_betas)[names(dasy_betas) == "beta_hat"] <- "estimated_beta"
  dasy_betas$std_error <- (dasy_betas$ci_upper - dasy_betas$ci_lower) / (2 * 1.96)
  dasy_betas$true_beta <- true_values
  dasy_betas$bias <- dasy_betas$estimated_beta - dasy_betas$true_beta
  dasy_betas$relative_bias <- ((dasy_betas$estimated_beta - dasy_betas$true_beta) / 
                                 dasy_betas$true_beta) * 100
  dasy_betas$within_ci <- (dasy_betas$true_beta >= dasy_betas$ci_lower & 
                             dasy_betas$true_beta <= dasy_betas$ci_upper)
  dasy_betas$method <- "Dasymetric"
  
  # Combine results
  combined_comparison <- rbind(abrm_betas, dasy_betas)
  
  # Create result with S3 class
  result <- list(
    combined_comparison = combined_comparison,
    abrm_results = abrm_full_results,
    dasymetric_results = dasymetric_results
  )
  class(result) <- "abrm_comparison"
  return(result)
}

#' Run Sensitivity Analysis
#'
#' Performs sensitivity analysis across different correlation structures
#'
#' @param correlation_grid Vector of correlation values to test
#' @param n_sims_per_setting Number of simulations per correlation setting
#' @param base_params List of base simulation parameters
#' @param mcmc_params List of MCMC parameters
#' @param model_code NIMBLE model code
#' @param base_seed Base random seed
#' @param output_dir Output directory for results (default: NULL, uses tempdir())
#'
#' @return List with combined results, summary statistics, and output directory
#' @export
#' @importFrom dplyr %>% group_by summarize
run_sensitivity_analysis <- function(
    correlation_grid = c(0.2, 0.6),
    n_sims_per_setting = 3,
    base_params = list(
      dist_covariates_x = c('normal','poisson','binomial'),
      dist_covariates_y = c('normal','poisson','binomial'),
      dist_y = 'poisson',
      x_intercepts = c(4, -1, -1),
      y_intercepts = c(4, -1, -1),
      beta0_y = -1,
      beta_x = c(-0.03, 0.1, -0.2),
      beta_y = c(0.03, -0.1, 0.2)
    ),
    mcmc_params = list(
      niter = 50000,
      nburnin = 30000,
      thin = 10,
      nchains = 2
    ),
    model_code,
    base_seed = 123,
    output_dir = NULL
) {
  
  # Create output directory
  if(is.null(output_dir)) {
    output_dir <- tempdir()
    message("Using temporary directory: ", output_dir)
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize results storage
  all_results <- list()
  counter <- 1
  
  # Determine outcome type from dist_y
  outcome_type <- base_params$dist_y
  
  # Run through correlation grid
  for(x_cor in correlation_grid) {
    for(y_cor in correlation_grid) {
      message(sprintf("\n=== Running analysis for x_cor = %.2f, y_cor = %.2f ===\n", x_cor, y_cor))
      
      for(sim in 1:n_sims_per_setting) {
        sim_seed <- base_seed + ((sim-1) * 100) + counter
        
        message(sprintf("\nSimulation %d of %d (seed: %d)\n", sim, n_sims_per_setting, sim_seed))
        
        # Create sub-directory for this simulation
        sim_dir <- file.path(output_dir, sprintf("xcor%.1f_ycor%.1f_sim%d", x_cor, y_cor, sim))
        dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
        
        # Generate simulated data
        message("  Generating data...\n")
        sim_data <- simulate_misaligned_data(
          seed = sim_seed,
          dist_covariates_x = base_params$dist_covariates_x,
          dist_covariates_y = base_params$dist_covariates_y,
          dist_y = base_params$dist_y,
          x_intercepts = base_params$x_intercepts,
          y_intercepts = base_params$y_intercepts,
          x_correlation = x_cor,
          y_correlation = y_cor,
          beta0_y = base_params$beta0_y,
          beta_x = base_params$beta_x,
          beta_y = base_params$beta_y
        )
        
        # Save the simulated dataset
        saveRDS(sim_data, file = file.path(sim_dir, "simulated_data.rds"))
        
        # Create metadata
        sim_metadata <- list(
          sim_number = sim,
          x_correlation = x_cor,
          y_correlation = y_cor,
          output_dir = sim_dir
        )
        
        # Run both methods
        results <- run_both_methods(
          sim_data = sim_data,
          sim_metadata = sim_metadata,
          model_code = model_code,
          nimble_params = mcmc_params,
          output_dir = sim_dir,
          norm_idx_x = which(base_params$dist_covariates_x == 'normal'),
          pois_idx_x = which(base_params$dist_covariates_x == 'poisson'),
          binom_idx_x = which(base_params$dist_covariates_x == 'binomial'),
          norm_idx_y = which(base_params$dist_covariates_y == 'normal'),
          pois_idx_y = which(base_params$dist_covariates_y == 'poisson'),
          binom_idx_y = which(base_params$dist_covariates_y == 'binomial'),
          dist_y = which(c('normal', 'poisson', 'binomial') == base_params$dist_y),
          outcome_type = outcome_type
        )
        
        # Add simulation metadata
        combined_comparison <- results$combined_comparison
        combined_comparison$sim_number <- sim
        combined_comparison$x_correlation <- x_cor
        combined_comparison$y_correlation <- y_cor
        
        # Save results for this simulation
        write.csv(combined_comparison, 
                  file = file.path(sim_dir, "method_comparison.csv"), 
                  row.names = FALSE)
        
        # Create comparison plots for this simulation
        create_comparison_plots(combined_comparison, sim_dir, sim_data$true_params)
        
        # Add to overall results
        all_results[[counter]] <- combined_comparison
        counter <- counter + 1
        
        message(sprintf("  Completed simulation %d\n", sim))
      }
    }
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_results)
  
  # Save combined results
  write.csv(combined_results, 
            file = file.path(output_dir, "sensitivity_analysis_results.csv"), 
            row.names = FALSE)
  
  # Create summary statistics
  message("\nCreating summary statistics...\n")
  summary_stats <- create_summary_statistics(combined_results, output_dir)
  
  # Create sensitivity summary plots
  message("\nCreating sensitivity summary plots...\n")
  create_sensitivity_summary_plots(combined_results, output_dir)
  
  # Aggregate results by correlation and method
  summary_by_correlation <- combined_results %>%
    group_by(method, x_correlation, y_correlation) %>%
    summarize(
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(relative_bias)),
      rmse = sqrt(mean(bias^2)),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  write.csv(summary_by_correlation, 
            file = file.path(output_dir, "sensitivity_summary_by_correlation.csv"), 
            row.names = FALSE)
  
  message("\n=== SENSITIVITY ANALYSIS COMPLETE ===\n")
  message("Results saved to:", output_dir, "\n")
  
  result <- list(
    combined_results = combined_results,
    summary_stats = summary_stats,
    summary_by_correlation = summary_by_correlation,
    output_dir = output_dir
  )
  class(result) <- "sensitivity_analysis"
  return(result)
}