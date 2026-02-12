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
#' @param seed Integer seed for reproducibility. Each chain uses seed+(chain_number-1) (default: NULL)
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
                     seed = NULL,
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
    seed = seed,
    save_plots = save_plots,
    output_dir = output_dir
  )
  
  # Rename MCMC sample columns to distinguish beta_x and beta_y
  p_x <- nimble_inputs$constants$p_x
  p_y <- nimble_inputs$constants$p_y
  
  for (chain_name in names(abrm_results$samples)) {
    old_names <- colnames(abrm_results$samples[[chain_name]])
    new_names <- old_names
    
    # Rename beta_y[1] through beta_y[p_x] to beta_x[1] through beta_x[p_x]
    for (i in 1:p_x) {
      old_pattern <- paste0("beta_y\\[", i, "\\]")
      new_pattern <- paste0("beta_x[", i, "]")
      new_names <- gsub(old_pattern, new_pattern, new_names, fixed = FALSE)
    }
    
    # Rename beta_y[p_x+1] through beta_y[p_x+p_y] to beta_y[1] through beta_y[p_y]
    for (i in 1:p_y) {
      old_idx <- p_x + i
      old_pattern <- paste0("beta_y\\[", old_idx, "\\]")
      new_pattern <- paste0("beta_y[", i, "]")
      new_names <- gsub(old_pattern, new_pattern, new_names, fixed = FALSE)
    }
    
    colnames(abrm_results$samples[[chain_name]]) <- new_names
  }
  
  # Also update summary table row names if it exists
  if (!is.null(abrm_results$summary$all.chains)) {
    old_rownames <- rownames(abrm_results$summary$all.chains)
    new_rownames <- old_rownames
    
    for (i in 1:p_x) {
      old_pattern <- paste0("beta_y\\[", i, "\\]")
      new_pattern <- paste0("beta_x[", i, "]")
      new_rownames <- gsub(old_pattern, new_pattern, new_rownames, fixed = FALSE)
    }
    
    for (i in 1:p_y) {
      old_idx <- p_x + i
      old_pattern <- paste0("beta_y\\[", old_idx, "\\]")
      new_pattern <- paste0("beta_y[", i, "]")
      new_rownames <- gsub(old_pattern, new_pattern, new_rownames, fixed = FALSE)
    }
    
    rownames(abrm_results$summary$all.chains) <- new_rownames
  }
  
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
    
    # Filter for beta parameters (both X and Y)
    beta_x_params <- grep("^beta_x\\[", abrm_parameters$variable)
    beta_y_params <- grep("^beta_y\\[", abrm_parameters$variable)
    beta_params <- c(beta_x_params, beta_y_params)
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
