#' R Wrapper Function for BiasedUrn Sampling
#'
#' Wraps the BiasedUrn::rMFNCHypergeo function for use in NIMBLE models
#'
#' @param total Integer, total number of items to sample
#' @param odds Numeric vector of odds for each category
#' @param ni Integer vector of population sizes
#' @return Numeric vector of sampled counts
#' @export
#' @keywords internal
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

#' R Wrapper for BiasedUrn Exact Density
#'
#' Computes the exact Fisher non-central hypergeometric density including the
#' normalising constant M(N,n,w). The original hand-rolled nimbleFunction
#' omitted M, which depends on the odds parameters and therefore biases the
#' MCMC posterior whenever binomial variables are present.
#'
#' @param x Vector of counts
#' @param total Total number of items
#' @param odds Numeric vector of odds
#' @param ni Integer vector of population sizes
#' @param log Integer 0/1; return log probability if 1
#' @return Numeric scalar density or log-density
#' @export
#' @keywords internal
biasedUrn_dmfnc <- function(x, total, odds, ni, log = 0L) {
  x   <- as.integer(round(x))
  ni  <- as.integer(round(ni))
  tot <- as.integer(round(total))
  odds[odds <= 0] <- 1e-10
  if (any(is.na(x)) || any(is.na(ni)) || is.na(tot)) {
    return(if (log == 1L) -Inf else 0)
  }
  if (tot > sum(ni)) tot <- sum(ni)
  if (length(ni) == 1L) return(if (log == 1L) 0 else 1)
  p <- BiasedUrn::dMFNCHypergeo(x, ni, tot, odds)
  if (log == 1L) return(log(max(p, .Machine$double.xmin)))
  return(p)
}

#' Density Function for Multivariate Non-Central Hypergeometric Distribution
#'
#' Thin NIMBLE wrapper around \code{biasedUrn_dmfnc}, which delegates to
#' \code{BiasedUrn::dMFNCHypergeo} for a correctly normalised density.
#'
#' @param x Vector of counts
#' @param total Total number of items
#' @param odds Vector of odds
#' @param ni Vector of category sizes
#' @param log Logical, return log probability
#' @return The log-probability (if log=1) or probability (if log=0)
#' @export
#' @keywords internal
#' @import nimble
dmfnchypg <- nimble::nimbleFunction(
  run = function(x = double(1), total = double(0), odds = double(1),
                 ni = double(1), log = integer(0, default = 0)) {
    returnType(double(0))
    return(Rdmfnchypg(x, total, odds, ni, log))
  }
)

#' NIMBLE R Call Wrapper for BiasedUrn Sampler
#'
#' Internal wrapper to call R function from compiled Nimble code.
#'
#' @param total Total number of items
#' @param odds Vector of odds
#' @param ni Vector of category sizes
#' @return Vector of sampled counts
#' @export
#' @keywords internal
#' @rdname Rmfnchypg_rcall
Rmfnchypg <- nimble::nimbleRcall(
  prototype = function(total = double(0), odds = double(1), ni = double(1)) {},
  returnType = double(1),
  Rfun = "biasedUrn_rmfnc"
)

#' NIMBLE R Call Wrapper for BiasedUrn Exact Density
#'
#' Calls \code{biasedUrn_dmfnc} (which uses \code{BiasedUrn::dMFNCHypergeo})
#' from compiled NIMBLE code, providing a correctly normalised log-likelihood.
#'
#' @param x Vector of counts
#' @param total Total number of items
#' @param odds Vector of odds
#' @param ni Vector of category sizes
#' @param log Integer 0/1
#' @return Scalar density or log-density
#' @export
#' @keywords internal
Rdmfnchypg <- nimble::nimbleRcall(
  prototype = function(x = double(1), total = double(0),
                       odds = double(1), ni = double(1),
                       log = integer(0)) {},
  returnType = double(0),
  Rfun = "biasedUrn_dmfnc"
)

#' Random Generation for Multivariate Non-Central Hypergeometric Distribution
#' 
#' @param n number of observations (only n=1 is used)
#' @param total Total number of items
#' @param odds Vector of odds
#' @param ni Vector of category sizes
#' @return Vector of sampled counts
#' @export
#' @keywords internal
#' @rdname rmfnchypg
rmfnchypg <- nimble::nimbleFunction(
  run = function(n = integer(0), total = double(0), odds = double(1), 
                 ni = double(1)) {
    returnType(double(1))
    return(Rmfnchypg(total, odds, ni))
  }
)

#' Register Custom NIMBLE Distributions
#' Registers the custom distributions for use in NIMBLE models.
#' @return Invisible TRUE
#' @export
#' @keywords internal
register_nimble_distributions <- function() {
  if (exists("distributions_registered", envir = .pkg_env) && 
      .pkg_env$distributions_registered) {
    return(invisible(TRUE))
  }
  
  # Only the registration call here - functions are already exported above
  nimble::registerDistributions(list(
    dmfnchypg = list(
      BUGSdist = "dmfnchypg(total, odds, ni)",
      discrete = TRUE,
      types = c('value = double(1)', 'total = double(0)', 
                'odds = double(1)', 'ni = double(1)'),
      mixedSizes = TRUE,
      pqAvail = FALSE,
      range = c(0, Inf)
    )
  ))
  
  .pkg_env$distributions_registered <- TRUE
  invisible(TRUE)
}

#' Get ABRM Model Code for NIMBLE
#'
#' Returns the NIMBLE code for the Atom-Based Regression Model with mixed-type variables.
#' Automatically registers custom distributions if not already registered.
#' @return A \code{nimbleCode} object containing the NIMBLE model specification
#'   for the ABRM. Pass this directly to \code{\link{run_abrm}}.
#' @examples
#' \donttest{
#'   model_code <- get_abrm_model()
#'   print(model_code)
#' }
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
    
    if (p_a>0){
      for (k in 1:p_a){
        for (j in 1:p_x){
          beta_a_x[k,j] ~ dnorm(0, sd = 2)
        }
      }
      
      for (k in 1:p_a){
        for (j in 1:p_y){
          beta_a_yx[k,j] ~ dnorm(0, sd = 2)
        }
      }
      
      for (j in (p_x+p_y+1):(p_x+p_y+p_a)){
        beta_y[j] ~ dnorm(0, sd = 1)
      }
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
    if (p_a>0){
      if (p_a==1){
        for(j in 1:p_x) {
          for(d in 1:D) {
            # Linear predictor (common to all types)
            linear_pred_x[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + covars_atom_x[d, 1] * beta_a_x[1,j]
          }
        }
      } else{
        for(j in 1:p_x) {
          for(d in 1:D) {
            # Linear predictor (common to all types)
            linear_pred_x[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + sum(covars_atom_x[d, 1:p_a] * beta_a_x[1:p_a,j])
          }
        }
      }
    } else{
      for(j in 1:p_x) {
        for(d in 1:D) {
          # Linear predictor (common to all types)
          linear_pred_x[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j]
        }
      }
    }
    
    # Observed X values (atom-equivalent observations)
    if (norm_n_x>0){
      for (j in 1:norm_n_x){
        for(i in 1:J_x) {
          # Normal
          x[i,norm_idx_x[j]] ~ dnorm(pop_atoms_x[i,norm_idx_x[j]] * linear_pred_x[i,norm_idx_x[j]],
                                     sd = pop_atoms_x[i,norm_idx_x[j]] * sqrt(sigma2_x[j]))

        }

        for(m in 1:(S_x-J_x)) {
          # Normal: Use Cong et al. algorithm
          mu_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), j]<- pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),norm_idx_x[j]] * linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]), norm_idx_x[j]]

          cov_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), (xlatent_ind[m,1]):(xlatent_ind[m,2]), j] <-sigma2_x[j] * diag(pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),norm_idx_x[j]]^2)

          # Sample from unconstrained multivariate normal
          w_x[xlatent_ind[m,1]:xlatent_ind[m,2], j] ~ dmnorm(
            mu_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), j],
            cov = cov_norm_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]), (xlatent_ind[m,1]):(xlatent_ind[m,2]), j]
          )

          # Project onto constraint hyperplane (population-weighted residual distribution)
          x_latent[xlatent_ind[m,1]:xlatent_ind[m,2], norm_idx_x[j]] <-
            w_x[xlatent_ind[m,1]:xlatent_ind[m,2], j] +
            (pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),norm_idx_x[j]] /
             sum(pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),norm_idx_x[j]])) *
            (x[J_x+m,norm_idx_x[j]] - sum(w_x[xlatent_ind[m,1]:xlatent_ind[m,2], j]))

        }
      }
    }
    
    if (pois_n_x>0){
      for (j in 1:pois_n_x){
        for(i in 1:J_x) {
          x[i,pois_idx_x[j]] ~ dpois(pop_atoms_x[i,pois_idx_x[j]]*exp(linear_pred_x[i,pois_idx_x[j]]))
        }
        for(m in 1:(S_x-J_x)) {
          # Poisson: Use multinomial

          prob_pois_x[(xlatent_ind[m,1]):(xlatent_ind[m,2]),j] <- (pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),pois_idx_x[j]] * exp(linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),pois_idx_x[j]])) /
            sum(pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),pois_idx_x[j]] * exp(linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),pois_idx_x[j]]))

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
          x[i,binom_idx_x[j]] ~ dbinom(size = pop_atoms_x[i,binom_idx_x[j]],
                                       prob = exp(linear_pred_x[i,binom_idx_x[j]])/(1+exp(linear_pred_x[i,binom_idx_x[j]])))
        }
        for(m in 1:(S_x-J_x)) {
          # Binomial: Use multivariate non-central hypergeometric
          x_latent[xlatent_ind[m,1]:xlatent_ind[m,2], binom_idx_x[j]] ~
            dmfnchypg(total = x[J_x+m,binom_idx_x[j]],
                      odds = exp(linear_pred_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),binom_idx_x[j]]),
                      ni = pop_atoms_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),binom_idx_x[j]])
        }
      }
    }
    
    for (j in 1:p_x){
      for(i in 1:J_x) {
        temp_x[i,j] <- x[i,j] / pop_atoms_x[i,j]
      }
      for(i in 1:(D-J_x)) {
        temp_x[J_x+i,j] <- x_latent[i,j] / pop_atoms_x[J_x+i,j]
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
    if (p_a>0){
      if (p_a==1){
        for(j in 1:p_y) {
          for(d in 1:D) {
            # Linear predictor (common to all types)
            linear_pred_yx[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + covars_atom_y[d, 1] * beta_a_yx[1,j]
          }
        }
      } else{
        for(j in 1:p_y) {
          for(d in 1:D) {
            # Linear predictor (common to all types)
            linear_pred_yx[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + sum(covars_atom_y[d, 1:p_a] * beta_a_yx[1:p_a,j])
          }
        }
      }
    } else{
      for(j in 1:p_y) {
        for(d in 1:D) {
          # Linear predictor (common to all types)
          linear_pred_yx[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j]
        }
      }
    }
    
    
    if (norm_n_y>0){
      for(j in 1:norm_n_y) {
        for(i in 1:J_y) {
          # Normal
          yx_obs[i,norm_idx_y[j]] ~ dnorm(pop_atoms_y[i,norm_idx_y[j]] * linear_pred_yx[i,norm_idx_y[j]],
                                          sd = pop_atoms_y[i,norm_idx_y[j]] * sqrt(sigma2_yx[j]))
        }
        for(m in 1:(S_y-J_y)) {
          # Normal: Use Cong et al. algorithm
          mu_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]), j] <-  pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),norm_idx_y[j]] * linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]), norm_idx_y[j]]

          cov_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2]),j]<- sigma2_yx[j] * diag(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),norm_idx_y[j]]^2)

          # Sample from unconstrained multivariate normal
          w_yx[ylatent_ind[m,1]:ylatent_ind[m,2], j] ~ dmnorm(
            mu_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]), j],
            cov = cov_norm_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2]),j]
          )

          # Project onto constraint hyperplane (population-weighted residual distribution)
          yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2], norm_idx_y[j]] <-
            w_yx[ylatent_ind[m,1]:ylatent_ind[m,2],j] +
            (pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),norm_idx_y[j]] /
             sum(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),norm_idx_y[j]])) *
            (yx_obs[J_y+m,norm_idx_y[j]] - sum(w_yx[ylatent_ind[m,1]:ylatent_ind[m,2], j]))

        }
      }
    }
    
    if (pois_n_y>0){
      for(j in 1:pois_n_y) {
        for(i in 1:J_y) {
          # Poisson
          yx_obs[i,pois_idx_y[j]] ~ dpois(pop_atoms_y[i,pois_idx_y[j]]*exp(linear_pred_yx[i,pois_idx_y[j]]))
        }
        for(m in 1:(S_y-J_y)) {
          # Poisson: Use multinomial

          prob_pois_yx[(ylatent_ind[m,1]):(ylatent_ind[m,2]),j]<-(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),pois_idx_y[j]] * exp(linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),pois_idx_y[j]])) /
            sum(pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),pois_idx_y[j]] * exp(linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),pois_idx_y[j]]))

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
          yx_obs[i,binom_idx_y[j]] ~ dbinom(size = pop_atoms_y[i,binom_idx_y[j]],
                                            prob = exp(linear_pred_yx[i,binom_idx_y[j]])/(1+exp(linear_pred_yx[i,binom_idx_y[j]])))
        }
        for(m in 1:(S_y-J_y)) {
          # Binomial: Use multivariate non-central hypergeometric
          yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2], binom_idx_y[j]] ~
            dmfnchypg(total = yx_obs[J_y+m,binom_idx_y[j]],
                      odds = exp(linear_pred_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),binom_idx_y[j]]),
                      ni = pop_atoms_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),binom_idx_y[j]])
        }
      }
    }
    
    for (j in 1:p_y){
      for(i in 1:J_y) {
        temp_yx[i,j] <- yx_obs[i,j] / pop_atoms_y[i,j]
      }
      for(i in 1:(D-J_y)) {
        temp_yx[J_y+i,j] <- yx_latent[i,j] / pop_atoms_y[J_y+i,j]
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
      
      if (p_a>0){
        if (p_a==1){
          a_effect_sum[d] <- covars_atom_y[d, 1] * beta_y[(p_x+p_y+1)]
        } else{
          a_effect_sum[d] <- sum(covars_atom_y[d, 1:p_a] * beta_y[(p_x+p_y+1):(p_x+p_y+p_a)])
        }
        linear_pred_y[d] <- beta_0_y + x_effect_sum[d] + y_effect_sum[d] + a_effect_sum[d] + phi_y[expand_y[d]]
      } else{
        linear_pred_y[d] <- beta_0_y + x_effect_sum[d] + y_effect_sum[d] + phi_y[expand_y[d]]
      }
    }
    
    # Observed Y outcome values
    for(i in 1:J_y) {
      if(vartype_y == 1) {
        # Normal
        y_obs[i] ~ dnorm(pop_atoms_outcome[i] * linear_pred_y[i],
                         sd = pop_atoms_outcome[i] * sqrt(sigma2_y))
      } else if(vartype_y == 2) {
        # Poisson
        y_obs[i] ~ dpois(pop_atoms_outcome[i]*exp(linear_pred_y[i]))
      } else if(vartype_y == 3) {
        y_obs[i] ~ dbinom(size = pop_atoms_outcome[i],
                          prob = exp(linear_pred_y[i])/(1+exp(linear_pred_y[i])))
      }
    }

    # Latent Y outcome values
    for(m in 1:(S_y-J_y)) {
      if(vartype_y == 1) {
        # Normal: Use Cong et al. algorithm
        mu_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])]<- pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]

        cov_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2])]<-sigma2_y * diag(pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]^2)

        w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ dmnorm(
          mu_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])],
          cov = cov_norm_y[(ylatent_ind[m,1]):(ylatent_ind[m,2]),(ylatent_ind[m,1]):(ylatent_ind[m,2])]
        )

        # Project onto constraint hyperplane (population-weighted residual distribution)
        y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] <-
          w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] +
          (pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] /
           sum(pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])) *
          (y_obs[J_y+m] - sum(w_y[ylatent_ind[m,1]:ylatent_ind[m,2]]))

      } else if(vartype_y == 2) {
        # Poisson: Use multinomial

        prob_pois_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])]<-(pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * exp(linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])) /
          sum(pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])] * exp(linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]))

        y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~
          dmulti(size = y_obs[J_y+m],
                 prob =prob_pois_y[(ylatent_ind[m,1]):(ylatent_ind[m,2])])

      } else if(vartype_y == 3) {
        # Binomial: Use multivariate non-central hypergeometric
        y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~
          dmfnchypg(total = y_obs[J_y+m],
                    odds = exp(linear_pred_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]),
                    ni = pop_atoms_outcome[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])
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
#' @param seed Integer seed for reproducibility. Each chain uses seed+(chain_number-1) (default: NULL)
#' @param save_plots Logical, whether to save trace and density plots to a PDF file (default: FALSE)
#' @param output_dir Directory for saving plots (default: NULL)
#' @param compute_waic Logical; if \code{TRUE}, WAIC is computed by NIMBLE and stored in the result, enabling \code{\link{waic}()} for model comparison. Default \code{FALSE} to keep fitting fast.
#' @return A named list with components \code{samples} (per-chain MCMC sample
#'   matrices, with the outcome coefficients named \code{beta_0_y},
#'   \code{beta_x[]}, \code{beta_y[]} and \code{beta_a[]}), \code{summary}
#'   (posterior summary statistics), \code{convergence} (Gelman-Rubin
#'   statistics and effective sample sizes; see \code{\link{mcmc_diagnostics}})
#'   and, when \code{compute_waic = TRUE}, \code{WAIC}.
#' @examples
#' \donttest{
#'   # The recommended workflow is to call run_abrm(), which calls
#'   # run_nimble_model() internally after preparing all inputs.
#'   sim_data <- simulate_misaligned_data(
#'     seed = 1,
#'     dist_covariates_x = "norm",
#'     dist_covariates_y = "norm",
#'     dist_y = "norm",
#'     x_intercepts = 0,
#'     y_intercepts = 0,
#'     beta0_y = 0,
#'     beta_x = 0.1,
#'     beta_y = -0.1
#'   )
#'   model_code <- get_abrm_model()
#'   results <- run_abrm(
#'     gridx = sim_data$gridx,
#'     gridy = sim_data$gridy,
#'     atoms = sim_data$atoms,
#'     model_code = model_code,
#'     x_types = "norm",
#'     y_types = "norm",
#'     dist_y = "norm",
#'     true_params = sim_data$true_params,
#'     niter = 1000, nburnin = 500, nchains = 2,
#'     seed = 1, save_plots = FALSE
#'   )
#'   summary(results)
#' }
#'
#' @export
#' @import nimble
#' @importFrom grDevices pdf dev.off
run_nimble_model <- function(constants, data, inits, sim_metadata = NULL, 
                             model_code, niter = 50000, nburnin = 30000, 
                             nchains = 2, thin = 10, seed = NULL,
                             save_plots = FALSE, output_dir = NULL,
                             compute_waic = FALSE) {
  
  register_nimble_distributions()
  
  # NIMBLE's internal functions (e.g. getNimbleOption, called from
  # init_isDataEnv) require nimble to be *attached* to the search path, not
  # merely imported into the package namespace via NAMESPACE/Imports.
  # Explicitly attaching it here prevents the
  # "could not find function 'getNimbleOption'" error during nimbleMCMC init.
  if (!"package:nimble" %in% search()) {
    suppressPackageStartupMessages(
      do.call("require", list(package = "nimble", quietly = TRUE))
    )
  }
  
  if(!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("Running NIMBLE MCMC...")
  
  mcmc.out <- nimble::nimbleMCMC(
    code = model_code,
    constants = constants,
    data = data,
    inits = inits,
    monitors = c('beta_0_y', 'beta_y', 'linear_pred_y'),
    thin = thin,
    niter = niter,
    setSeed = if(!is.null(seed)) seed + (0:(nchains-1)) else FALSE,
    nburnin = nburnin,
    nchains = nchains,
    summary = TRUE,
    WAIC = compute_waic
  )
  
  # With a single chain NIMBLE returns bare matrices; normalise to the
  # multi-chain list layout so downstream code has one shape to handle.
  if (is.matrix(mcmc.out$samples)) {
    mcmc.out$samples <- list(chain1 = mcmc.out$samples)
  }
  if (is.matrix(mcmc.out$summary)) {
    mcmc.out$summary <- list(chain1 = mcmc.out$summary, all.chains = mcmc.out$summary)
  }

  # NIMBLE stores all outcome coefficients in one vector beta_y[1:(p_x+p_y+p_a)].
  # Split it into beta_x[], beta_y[] and beta_a[] once, here, so samples,
  # summaries and diagnostics all use the same names.
  p_a <- if (is.null(constants$p_a)) 0L else constants$p_a
  rename <- function(nms) .rename_beta_params(nms, constants$p_x, constants$p_y, p_a)
  for (ch in names(mcmc.out$samples)) {
    colnames(mcmc.out$samples[[ch]]) <- rename(colnames(mcmc.out$samples[[ch]]))
  }
  for (nm in names(mcmc.out$summary)) {
    rownames(mcmc.out$summary[[nm]]) <- rename(rownames(mcmc.out$summary[[nm]]))
  }

  # Strip linear_pred_y from MCMC samples before diagnostics.
  # The summary table keeps these rows so run_abrm() can compute fitted values,
  # but they must not appear in trace/density plots, vcov(), or confint().
  for (ch in names(mcmc.out$samples)) {
    keep <- !grepl("^linear_pred_y", colnames(mcmc.out$samples[[ch]]))
    mcmc.out$samples[[ch]] <- mcmc.out$samples[[ch]][, keep, drop = FALSE]
  }

  message("Calculating convergence diagnostics...")
  diagnostics <- check_mcmc_diagnostics(mcmc.out)

  print_convergence_summary(diagnostics)

  if (save_plots) {
    if (is.null(sim_metadata)) {
      plot_file <- "mcmc_diagnostics.pdf"
      plot_title <- "MCMC diagnostics"
    } else {
      plot_file <- sprintf("mcmc_diagnostics_sim%d_xcor%.1f_ycor%.1f.pdf",
                           as.integer(sim_metadata$sim_number),
                           sim_metadata$x_correlation,
                           sim_metadata$y_correlation)
      plot_title <- sprintf("MCMC diagnostics (x_cor=%.1f, y_cor=%.1f, sim=%d)",
                            sim_metadata$x_correlation, sim_metadata$y_correlation,
                            as.integer(sim_metadata$sim_number))
    }
    if (!is.null(output_dir)) {
      plot_file <- file.path(output_dir, plot_file)
    }
    grDevices::pdf(plot_file, width = 12, height = 8)
    print(create_trace_plot(mcmc.out$samples, title = plot_title))
    print(create_density_plot(mcmc.out$samples, title = plot_title))
    grDevices::dev.off()
    message("Diagnostic plots saved to ", plot_file)
  }

  mcmc.out$convergence <- diagnostics
  
  return(mcmc.out)
}