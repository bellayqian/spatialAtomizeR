#' R Wrapper Function for BiasedUrn Sampling
#'
#' Wraps the BiasedUrn::rMFNCHypergeo function for use in NIMBLE models
#'
#' @param total Integer, total number of items to sample
#' @param odds Numeric vector of odds for each category
#' @param ni Integer vector of population sizes
#' @return Numeric vector of sampled counts
#' @export
#' @importFrom BiasedUrn rMFNCHypergeo
biasedUrn_rmfnc <- function(total, odds, ni) {
  total <- as.integer(round(total))
  ni <- as.integer(round(ni))
  
  odds[odds <= 0] <- 1e-10
  
  sum_ni <- sum(ni)
  if(total > sum_ni) {
    warning(sprintf("Total (%d) exceeds sum of population sizes (%d). Adjusting total.",
                    total, sum_ni))
    total <- sum_ni
  }
  
  if(length(ni) == 1) {
    return(total)
  }
  
  if(any(is.na(ni)) || any(is.na(odds)) || is.na(total)) {
    stop("NA values in inputs")
  }
  if(any(ni < 0) || any(odds < 0) || total < 0) {
    stop("Negative values in inputs")
  }
  
  result <- BiasedUrn::rMFNCHypergeo(1, ni, total, odds)
  return(as.numeric(result))
}

#' Register Custom NIMBLE Distributions
#'
#' Registers the custom multivariate non-central hypergeometric distribution
#' for use in NIMBLE models. This function is called automatically when needed.
#' 
#' @note The <<- operator is used intentionally to create package-level 
#' nimbleFunctions accessible across the package environment.
#'
#' @return Invisible TRUE if successful
#' @export
#' @importFrom nimble nimbleFunction nimbleRcall registerDistributions
register_nimble_distributions <- function() {
  if (.pkg_env$distributions_registered) {
    return(invisible(TRUE))
  }
  
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("nimble package is required but not installed")
  }
  
  # Define the dmfnchypg density function
  .pkg_env$dmfnchypg <- nimble::nimbleFunction(
    run = function(x = double(1), total = double(0), odds = double(1), ni = double(1),
                   log = integer(0)) {
      returnType(double(0))
      
      x_round <- round(x)
      total_round <- round(total)
      
      sumX <- 0
      for(i in 1:length(x)) {
        sumX <- sumX + x[i]
      }
      
      if(abs(sumX - total) > 0.1) {
        if(log == 1) return(-Inf) else return(0)
      }
      
      for(i in 1:length(x)) {
        if(x[i] < 0) {
          if(log == 1) return(-Inf) else return(0)
        }
        if(x_round[i] > round(ni[i])) {
          if(log == 1) return(-Inf) else return(0)
        }
      }
      
      logProb <- 0
      
      for(i in 1:length(ni)) {
        logProb <- logProb + lfactorial(round(ni[i]))
      }
      
      for(i in 1:length(x_round)) {
        logProb <- logProb - lfactorial(x_round[i])
      }
      
      for(i in 1:length(x_round)) {
        logProb <- logProb - lfactorial(round(ni[i]) - x_round[i])
      }
      
      for(i in 1:length(x_round)) {
        safe_odds <- max(odds[i], 1e-10)
        logProb <- logProb + x_round[i] * log(safe_odds)
      }
      
      if(log == 1) return(logProb) else return(exp(logProb))
    }
  )
  
  # Create nimbleRcall wrapper
  .pkg_env$Rmfnchypg <- nimble::nimbleRcall(
    prototype = function(total = double(0), odds = double(1), ni = double(1)) {},
    returnType = double(1),
    Rfun = "biasedUrn_rmfnc"
  )
  
  # Create the nimbleFunction wrapper for random generation
  .pkg_env$rmfnchypg <- nimble::nimbleFunction(
    run = function(n = integer(0), total = double(0), odds = double(1), ni = double(1)) {
      returnType(double(1))
      return(.pkg_env$Rmfnchypg(total, odds, ni))
    }
  )
  
  # Register the custom distribution
  nimble::registerDistributions(list(
    dmfnchypg = list(
      BUGSdist = "dmfnchypg(total, odds, ni)",
      discrete = TRUE,
      types = c('value = double(1)', 'total = double(0)', 'odds = double(1)', 'ni = double(1)'),
      mixedSizes = TRUE,
      pqAvail = FALSE,
      range = c(0, Inf)
    )
  ))
  
  .pkg_env$distributions_registered <- TRUE
  message("NIMBLE custom distributions registered successfully")
  invisible(TRUE)
}

#' Get ABRM Model Code for NIMBLE
#'
#' Returns the NIMBLE code for the Atom-Based Regression Model with mixed-type variables.
#' Automatically registers custom distributions if not already registered.
#'
#' @return A nimbleCode object containing the model specification
#' @export
get_abrm_model <- function() {
  
  # Ensure distributions are registered
  register_nimble_distributions()
  
  mixed_abrm <- nimble::nimbleCode({
    
    ## PRIORS
    for(j in 1:p_x) {
      beta_0_x[j] ~ dnorm(0, sd = 2)
    }
    
    for(j in 1:p_y) {
      beta_0_yx[j] ~ dnorm(0, sd = 2)
    }
    
    beta_0_y ~ dnorm(0, sd = 2)
    
    for (j in 1:(p_x+p_y)){
      beta_y[j] ~ dnorm(0, sd = 1)
    }
    
    if (norm_n_x>0){
      for (j in 1:norm_n_x){
        sigma2_x[j] ~ dinvgamma(2, 1)
      }
    }
    
    if (norm_n_y>0){
      for(j in 1:norm_n_y) {
        sigma2_yx[j] ~ dinvgamma(2, 1)
      }
    }
    
    if(vartype_y == 1) {
      sigma2_y ~ dinvgamma(2, 1)
    }
    
    # Spatial random effects for outcome Y
    tau_y ~ dgamma(2, 2)
    phi_y[1:S_y] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_y)
    
    # Spatial random effects for X-grid covariates
    if (p_x>1){
      for (j in 1:p_x){
        tau_x[j] ~ dgamma(2, 2)
        psi_x[1:S_x,j] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x[j])
      }
      
      # Correlation in X spatial effects
      Prec_x[1:p_x,1:p_x] ~ dwish(R_x[1:p_x,1:p_x], df_x)
      Achol_x[1:p_x,1:p_x] <- chol(inverse(Prec_x[1:p_x,1:p_x]))
      
      for(i in 1:S_x){
        for(j in 1:p_x) {
          phi_x[i,j] <- inprod(Achol_x[j,1:p_x], psi_x[i,1:p_x])
        }
      }
    } else{
      for(j in 1:p_x) {
        tau_x[j] ~ dgamma(2, 2)
        phi_x[1:S_x,j] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x[j])
      }
    }
    
    
    if (p_y>1){
      # Spatial random effects for Y-grid covariates
      for (j in 1:p_y){
        tau_yx[j] ~ dgamma(2, 2)
        psi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j])
      }
      
      # Correlation in Y spatial effects
      Prec_yx[1:p_y,1:p_y] ~ dwish(R_yx[1:p_y,1:p_y], df_yx)
      Achol_yx[1:p_y,1:p_y] <- chol(inverse(Prec_yx[1:p_y,1:p_y]))
      
      for(i in 1:S_y){
        for(j in 1:p_y) {
          phi_yx[i,j] <- inprod(Achol_yx[j,1:p_y], psi_yx[i,1:p_y])
        }
      }
    } else{
      for(j in 1:p_y) {
        tau_yx[j] ~ dgamma(2, 2)
        phi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j])
      }
    }
    
    
    ################################
    ## MODELING X-GRID COVARIATES ##
    ################################
    
    # Calculate type-specific parameters for X variables
    for(j in 1:p_x) {
      for(d in 1:D) {
        # Linear predictor (common to all types)
        linear_pred_x[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j]
      }
    }
    
    # Observed X values (atom-equivalent observations)
    if (norm_n_x>0){
      for (j in 1:norm_n_x){
        for(i in 1:J_x) {
          # Normal
          x[i,norm_idx_x[j]] ~ dnorm(pop_atoms_x[i] * linear_pred_x[i,norm_idx_x[j]],
                                     sd = pop_atoms_x[i] * sqrt(sigma2_x[j]))
          
        }
        
        for(m in 1:(S_x-J_x)) {
          # Normal: Use Cong et al. algorithm
          nat_x[m,j] <- length(xlatent_ind[m,1]:xlatent_ind[m,2])
          
          mu_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), j]<- pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])] * linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]), norm_idx_x[j]]
          
          cov_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), (xlatent_ind[m,1]):(xlatent_ind[m,2]), j] <-sigma2_x[j] * diag(pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])]^2)
          
          # Sample from unconstrained multivariate normal
          w_x[xlatent_ind[m,1]:xlatent_ind[m,2], j] ~ dmnorm(
            mu_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), j],
            cov = cov_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), (xlatent_ind[m,1]):(xlatent_ind[m,2]), j]
          )
          
          # Project onto constraint hyperplane
          x_latent[xlatent_ind[m,1]:xlatent_ind[m,2], norm_idx_x[j]] <-
            w_x[xlatent_ind[m,1]:xlatent_ind[m,2], j] +
            (1/nat_x[m,j]) * (x[J_x+m,norm_idx_x[j]] - sum(w_x[xlatent_ind[m,1]:xlatent_ind[m,2], j])) * rep(1, nat_x[m,j])
          
        }
      }
    }
    
    if (pois_n_x>0){
      for (j in 1:pois_n_x){
        for(i in 1:J_x) {
          x[i,pois_idx_x[j]] ~ dpois(pop_atoms_x[i]*exp(linear_pred_x[i,pois_idx_x[j]]))
        }
        for(m in 1:(S_x-J_x)) {
          # Poisson: Use multinomial
          
          prob_pois_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]),j] <- (pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])] * exp(linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),pois_idx_x[j]])) /
            sum(pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])] * exp(linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),pois_idx_x[j]]))
          
          x_latent[xlatent_ind[m,1]:xlatent_ind[m,2],pois_idx_x[j]] ~
            dmulti(size=x[J_x+m,pois_idx_x[j]],
                   prob=prob_pois_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]),j])
          
        }
      }
    }
    
    if (binom_n_x>0){
      for (j in 1:binom_n_x){
        for(i in 1:J_x) {
          # Binomial
          x[i,binom_idx_x[j]] ~ dbinom(size = pop_atoms_x[i],
                                       prob = exp(linear_pred_x[i,binom_idx_x[j]])/(1+exp(linear_pred_x[i,binom_idx_x[j]])))
        }
        for(m in 1:(S_x-J_x)) {
          # Binomial: Use multivariate non-central hypergeometric
          x_latent[xlatent_ind[m,1]:xlatent_ind[m,2], binom_idx_x[j]] ~
            dmfnchypg(total = x[J_x+m,binom_idx_x[j]],
                      odds = exp(linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),binom_idx_x[j]]),
                      ni = pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])])
        }
      }
    }
    
    for (j in 1:p_x){
      for(i in 1:J_x) {
        temp_x[i,j] <- x[i,j] / pop_atoms_x[i]  # Convert to per-capita
      }
      for(i in 1:(D-J_x)) {
        temp_x[J_x+i,j] <- x_latent[i,j] / pop_atoms_x[J_x+i]  # Convert to per-capita
      }
    }
    
    # Reorder X matrices
    for(j in 1:p_x) {
      for(d in 1:D) {
        x_atomord[d,j] <- temp_x[x_reorder[d],j]
      }
    }
    
    ################################
    ## MODELING Y-GRID COVARIATES ##
    ################################
    
    # Similar structure for Y covariates
    for(j in 1:p_y) {
      for(d in 1:D) {
        # Linear predictor (common to all types)
        linear_pred_yx[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j]
      }
    }
    
    
    if (norm_n_y>0){
      for(j in 1:norm_n_y) {
        for(i in 1:J_y) {
          # Normal
          yx_obs[i,norm_idx_y[j]] ~ dnorm(pop_atoms_y[i] * linear_pred_yx[i,norm_idx_y[j]],
                                          sd = pop_atoms_y[i] * sqrt(sigma2_yx[j]))
        }
        for(m in 1:(S_y-J_y)) {
          # Normal: Use Cong et al. algorithm
          nat_yx[m,j] <- length(ylatent_ind[m,1]:ylatent_ind[m,2])
          
          
          mu_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]), j] <-  pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]), norm_idx_y[j]]
          
          cov_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2]),j]<- sigma2_yx[j] * diag(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]^2)
          
          # Sample from unconstrained multivariate normal
          w_yx[ylatent_ind[m,1]:ylatent_ind[m,2], j] ~ dmnorm(
            mu_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]), j],
            cov = cov_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2]),j]
          )
          
          # Project onto constraint hyperplane
          yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2], norm_idx_y[j]] <-
            w_yx[ylatent_ind[m,1]:ylatent_ind[m,2],j] +
            (1/nat_yx[m,j]) * (yx_obs[J_y+m,norm_idx_y[j]] - sum(w_yx[ylatent_ind[m,1]:ylatent_ind[m,2], j])) * rep(1, nat_yx[m,j])
          
        }
      }
    }
    
    if (pois_n_y>0){
      for(j in 1:pois_n_y) {
        for(i in 1:J_y) {
          # Poisson
          yx_obs[i,pois_idx_y[j]] ~ dpois(pop_atoms_y[i]*exp(linear_pred_yx[i,pois_idx_y[j]]))
        }
        for(m in 1:(S_y-J_y)) {
          # Poisson: Use multinomial
          
          prob_pois_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),j]<-(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * exp(linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),pois_idx_y[j]])) /
            sum(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * exp(linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),pois_idx_y[j]]))
          
          yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],pois_idx_y[j]] ~
            dmulti(size=yx_obs[J_y+m,pois_idx_y[j]],
                   prob=prob_pois_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),j])
        }
      }
    }
    
    if (binom_n_y>0){
      for(j in 1:binom_n_y) {
        for(i in 1:J_y) {
          # Binomial
          yx_obs[i,binom_idx_y[j]] ~ dbinom(size = pop_atoms_y[i],
                                            prob = exp(linear_pred_yx[i,binom_idx_y[j]])/(1+exp(linear_pred_yx[i,binom_idx_y[j]])))
        }
        for(m in 1:(S_y-J_y)) {
          # Binomial: Use multivariate non-central hypergeometric
          yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2], binom_idx_y[j]] ~
            dmfnchypg(total = yx_obs[J_y+m,binom_idx_y[j]],
                      odds = exp(linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),binom_idx_y[j]]),
                      ni = pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])
        }
      }
    }
    
    for (j in 1:p_y){
      for(i in 1:J_y) {
        temp_yx[i,j] <- yx_obs[i,j] / pop_atoms_y[i]  # Convert to per-capita
      }
      for(i in 1:(D-J_y)) {
        temp_yx[J_y+i,j] <- yx_latent[i,j] / pop_atoms_y[J_y+i]  # Convert to per-capita
      }
    }
    
    #######################
    ## MODELING Y OUTCOME ##
    #######################
    
    # Calculate outcome parameters
    for(d in 1:D) {
      # Linear predictor with covariate effects
      if (p_x==1){
        x_effect_sum[d] <- x_atomord[y_to_atom[d], 1] * beta_y[1]
      } else{
        x_effect_sum[d] <- sum(x_atomord[y_to_atom[d], 1:p_x] * beta_y[1:p_x])
      }
      
      if (p_y==1){
        y_effect_sum[d] <- temp_yx[d, 1] * beta_y[(p_x+1)]
      } else{
        y_effect_sum[d] <- sum(temp_yx[d, 1:p_y] * beta_y[(p_x+1):(p_x+p_y)])
      }
      
      linear_pred_y[d] <- beta_0_y + x_effect_sum[d] + y_effect_sum[d] + phi_y[expand_y[d]]
    }
    
    # Observed Y outcome values
    for(i in 1:J_y) {
      if(vartype_y == 1) {
        # Normal
        y_obs[i] ~ dnorm(pop_atoms_y[i] * linear_pred_y[i],
                         sd = pop_atoms_y[i] * sqrt(sigma2_y))
      } else if(vartype_y == 2) {
        # Poisson
        y_obs[i] ~ dpois(pop_atoms_y[i]*exp(linear_pred_y[i]))
      } else if(vartype_y == 3) {
        y_obs[i] ~ dbinom(size = pop_atoms_y[i],
                          prob = exp(linear_pred_y[i])/(1+exp(linear_pred_y[i])))
      }
    }
    
    # Latent Y outcome values
    for(m in 1:(S_y-J_y)) {
      if(vartype_y == 1) {
        # Normal: Use Cong et al. algorithm
        nat_y[m] <- length(ylatent_ind[m,1]:ylatent_ind[m,2])
        
        mu_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])]<- pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]
        
        cov_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2])]<-sigma2_y * diag(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]^2)
        
        w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ dmnorm(
          mu_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])],
          cov = cov_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2])]
        )
        
        y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] <-
          w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] +
          (1/nat_y[m]) * (y_obs[J_y+m] - sum(w_y[ylatent_ind[m,1]:ylatent_ind[m,2]])) * rep(1, nat_y[m])
        
      } else if(vartype_y == 2) {
        # Poisson: Use multinomial
        
        prob_pois_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])]<-(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * exp(linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])) /
          sum(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * exp(linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]))
        
        y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~
          dmulti(size = y_obs[J_y+m],
                 prob =prob_pois_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])])
        
      } else if(vartype_y == 3) {
        # Binomial: Use multivariate non-central hypergeometric
        y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~
          dmfnchypg(total = y_obs[J_y+m],
                    odds = exp(linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]),
                    ni = pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])
      }
    }
  })
  
  return(mixed_abrm)
}

#' Run NIMBLE Model with Diagnostics
#'
#' @param constants List of model constants
#' @param data List of data
#' @param inits List of initial values
#' @param sim_metadata List with simulation metadata (optional)
#' @param model_code NIMBLE code object
#' @param niter Number of MCMC iterations (default: 50000)
#' @param nburnin Number of burn-in iterations (default: 30000)
#' @param nchains Number of MCMC chains (default: 2)
#' @param thin Thinning interval (default: 10)
#' @param save_plots Logical, whether to save diagnostic plots (default: TRUE)
#' @param output_dir Directory for saving plots (default: NULL)
#'
#' @return List containing MCMC samples, summary, and convergence diagnostics
#' @export
#' @importFrom nimble nimbleMCMC
#' @importFrom grDevices pdf dev.off
run_nimble_model <- function(constants, data, inits, sim_metadata = NULL, 
                             model_code, niter = 50000, nburnin = 30000, 
                             nchains = 2, thin = 10, 
                             save_plots = TRUE, output_dir = NULL) {
  
  register_nimble_distributions()
  
  if(!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("\nRunning NIMBLE MCMC...\n")
  
  mcmc.out <- nimble::nimbleMCMC(
    code = model_code,
    constants = constants,
    data = data,
    inits = inits,
    monitors = c('beta_0_y', 'beta_y'),
    thin = thin,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    summary = TRUE
  )
  
  message("\nCalculating convergence diagnostics...\n")
  diagnostics <- check_mcmc_diagnostics(mcmc.out, sim_metadata)
  
  print_convergence_summary(diagnostics)
  
  if(save_plots && !is.null(diagnostics$plots)) {
    if(is.null(sim_metadata)) {
      plot_file <- "mcmc_diagnostics.pdf"
    } else {
      plot_file <- sprintf("mcmc_diagnostics_sim%d_xcor%.1f_ycor%.1f.pdf",
                           sim_metadata$sim_number,
                           sim_metadata$x_correlation,
                           sim_metadata$y_correlation)
    }
    
    if(!is.null(output_dir)) {
      plot_file <- file.path(output_dir, plot_file)
    }
    
    grDevices::pdf(plot_file, width = 12, height = 8)
    print(diagnostics$plots$trace)
    print(diagnostics$plots$density)
    grDevices::dev.off()
    message("\nDiagnostic plots saved to", plot_file, "\n")
  }
  
  mcmc.out$convergence <- diagnostics
  
  return(mcmc.out)
}