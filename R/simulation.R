#' Generate Correlated Spatial Effects
#'
#' @param W Spatial adjacency matrix
#' @param n_vars Number of variables
#' @param rho Spatial correlation parameter (default = 0.6)
#' @param var_spat Spatial variance (default = 1)
#' @param correlation Correlation between variables (default = 0.5)
#' @param verify Logical for verification (default = FALSE)
#'
#' @return Matrix of spatial effects
#' @importFrom MASS mvrnorm
#' @importFrom stats rpois rbinom
#' @importFrom methods as
#' @export
gen_correlated_spat <- function(W, n_vars, rho = 0.6, var_spat = 1, correlation = 0.5, verify = FALSE) {
  n <- nrow(W)
  
  # Create precision matrix for variables (Lambda)
  Sigma_vars <- matrix(correlation, n_vars, n_vars)
  diag(Sigma_vars) <- 1
  Lambda <- solve(Sigma_vars)
  
  # Spatial precision matrix
  precision_matrix <- diag(colSums(W)) - rho * W
  
  # Create full precision matrix using Kronecker product
  full_prec <- kronecker(Lambda, precision_matrix)
  
  # Generate from MCAR distribution
  all_effects <- MASS::mvrnorm(1, mu = rep(0, n * n_vars), Sigma = solve(full_prec))
  spatial_effects <- matrix(all_effects, n, n_vars)
  
  return(spatial_effects)
}

#' Simulate Misaligned Spatial Data
#'
#' @param seed Random seed (default = 2)
#' @param dist_covariates_x Vector specifying distribution type for each synthetic X-grid covariate ('poisson', 'binomial', or 'normal')
#' @param dist_covariates_y Vector specifying distribution type for each synthetic Y-grid covariate ('poisson', 'binomial', or 'normal')
#' @param dist_y Distribution type for synthetic outcome variable (one of 'poisson', 'binomial', or 'normal')
#' @param x_intercepts Intercepts for X covariates
#' @param y_intercepts Intercepts for Y covariates
#' @param rho_x Spatial correlation parameter for X-grid covariates (0 to 1 with higher values yielding more spatial correlation, default = 0.6)
#' @param rho_y Spatial correlation parameter for Y-grid covariates and outcome (0 to 1 with higher values yielding more spatial correlation, default = 0.6)
#' @param x_correlation Between-variable correlation for all pairs of X-grid covariates (default = 0.5)
#' @param y_correlation Between-variable correlation for all pairs of Y-grid covariates (default = 0.5)
#' @param beta0_y Intercept for outcome model
#' @param beta_x Outcome model coefficients for X-grid covariates
#' @param beta_y Outcome model coefficients for Y-grid covariates
#' @param diff_pops Logical, indicating whether the atoms should be generated with different population sizes (diff_pops = TRUE) or a common population size (diff_pops = FALSE)
#' @param xy_cov_cor Logical, indicating whether the atom-level spatial random effects for X-grid and Y-grid covariates should be correlated (xy_cov_cor = TRUE) or not. When set to TRUE, the x_correlation and rho_x parameters are used to generate all covariates (separate correlation parameters are not allowed for X-grid and Y-grid covariates).
#'
#' @return List containing gridy, gridx, atoms, and true_params
#' @importFrom sp GridTopology SpatialGrid
#' @importFrom sf st_as_sf st_union st_sf st_contains
#' @importFrom spdep poly2nb nb2mat
#' @importFrom raster intersect
#' @importFrom stats plogis
#' @export
simulate_misaligned_data <- function(seed = 2,
                                     dist_covariates_x = c('normal','poisson','binomial'),
                                     dist_covariates_y = c('normal','poisson','binomial'),
                                     dist_y = 'poisson',
                                     x_intercepts = rep(0, 3),
                                     y_intercepts = rep(0, 3),
                                     rho_x = 0.6,
                                     rho_y = 0.6,
                                     x_correlation = 0.5, 
                                     y_correlation = 0.5,
                                     beta0_y = NULL,
                                     beta_x = NULL, 
                                     beta_y = NULL,
                                     diff_pops = TRUE,
                                     xy_cov_cor = FALSE) {
  set.seed(seed)
  
  n_covariates_x <- length(dist_covariates_x)
  n_covariates_y <- length(dist_covariates_y)
  
  res1 <- c(5, 5)    # Y grid resolution (coarser)
  res2 <- c(10, 10)  # X grid resolution (finer)
  
  # Create spatial grids
  grid1 <- sp::GridTopology(cellcentre.offset = c(0.5, 0.5),
                            cellsize = c(2, 2),
                            cells.dim = res1)
  sp_grid1 <- sp::SpatialGrid(grid1)
  
  grid2 <- sp::GridTopology(cellcentre.offset = c(0, 0),
                            cellsize = c(1, 1),
                            cells.dim = res2)
  sp_grid2 <- sp::SpatialGrid(grid2)
  
  # Convert grids to sf objects
  sp_grid2_poly <- sf::st_as_sf(as(sp_grid2, 'SpatialPolygons'))
  sp_grid1_poly <- sf::st_as_sf(as(sp_grid1, 'SpatialPolygons'))
  
  # Create non-nested misalignment
  pair_list <- list(c(2,3), c(4,5), c(6,7), c(8,9))
  union_ids <- matrix(0, res2[1], res2[2])
  id_num <- 1
  for (i in seq(1, res2[1], 2)) {
    union_ids[i, pair_list[[sample(1:4, size = 1)]]] <- id_num
    id_num <- id_num + 1
  }
  
  sp_grid2_poly$union_ids <- c(t(union_ids))
  
  # Union cells
  store_xunions <- NULL
  for (i in 1:max(sp_grid2_poly$union_ids)) {
    temp <- sp_grid2_poly[sp_grid2_poly$union_ids == i, ] %>%
      sf::st_union() %>%
      sf::st_sf()
    store_xunions <- rbind(store_xunions, temp)
  }
  
  # Merge remaining cells
  grid2_nounion <- subset(sp_grid2_poly, union_ids == 0)
  store_iunions <- NULL
  for (i in 1:nrow(sp_grid1_poly)) {
    temp <- grid2_nounion[c(sf::st_contains(sp_grid1_poly[i,], grid2_nounion, sparse = FALSE)) == TRUE,] %>%
      sf::st_union() %>%
      sf::st_sf()
    store_iunions <- rbind(store_iunions, temp)
  }
  
  grid2_final <- rbind(store_iunions, store_xunions)
  
  # Rename grids
  gridx <- grid2_final
  gridy <- sp_grid1_poly
  gridy$ID <- 1:nrow(gridy)
  gridx$ID <- 1:nrow(gridx)
  
  # Create atoms
  atoms <- sf::st_as_sf(raster::intersect(as(gridy, 'Spatial'), as(gridx, 'Spatial')))
  names(atoms)[which(names(atoms) == 'ID')] <- c("ID_y")
  names(atoms)[which(names(atoms) == 'ID.1')] <- c("ID_x")
  atoms$ID_atomorder <- 1:nrow(atoms)
  
  # Generate atom-level population
  if (diff_pops==TRUE){
    atoms$population <- rpois(nrow(atoms), lambda = 50) + 10
  } else{
    atoms$population <- 1
  }
  
  # Create spatial adjacency matrix for atoms
  atoms_sp <- as(atoms, "Spatial")
  neighbors_atoms <- spdep::poly2nb(atoms_sp, queen = TRUE)
  W_atoms <- spdep::nb2mat(neighbors_atoms, style = "B", zero.policy = TRUE)
  
  ## generate atom level spatial random effects for covariates
  if (xy_cov_cor == FALSE){
    x_atom_results <- gen_correlated_spat(W_atoms, n_covariates_x, correlation = x_correlation, rho=rho_x)
    y_atom_results <- gen_correlated_spat(W_atoms, n_covariates_y , correlation = y_correlation, rho = rho_y)
  } else{
    xy_atom_results <- gen_correlated_spat(W_atoms, n_covariates_x + n_covariates_y, correlation = x_correlation, rho=rho_x)
    x_atom_results <-  xy_atom_results[,1:n_covariates_x]
    y_atom_results <- xy_atom_results[,(n_covariates_x+1):ncol(xy_atom_results)]
  }
  
  # Generate X covariates at atom level
  for(i in 1:n_covariates_x) {
    if (dist_covariates_x[i] == 'poisson') {
      log_rate_percapita <- x_intercepts[i] + x_atom_results[,i]
      lambda_total <- atoms$population * exp(log_rate_percapita)
      atoms[[paste0("covariate_x_", i)]] <- rpois(nrow(atoms), lambda = lambda_total)
      
    } else if (dist_covariates_x[i] == 'binomial') {
      probabilities <- plogis(x_intercepts[i] + x_atom_results[,i])
      probabilities <- pmin(pmax(probabilities, 0.001), 0.999)
      atoms[[paste0("covariate_x_", i)]] <- rbinom(nrow(atoms), size = atoms$population, prob = probabilities)
      
    } else if (dist_covariates_x[i] == 'normal') {
      atoms[[paste0("covariate_x_", i)]] <- atoms$population * rnorm(nrow(atoms), mean = x_intercepts[i] + x_atom_results[,i], sd = 1)
      
    } else {
      stop(paste0('Distribution of X-grid covariate ', i, ' not recognized'))
    }
  }
  
  # Generate Y covariates at atom level
  
  for(i in 1:n_covariates_y) {
    if (dist_covariates_y[i] == 'poisson') {
      log_rate_percapita <- y_intercepts[i] + y_atom_results[,i]
      lambda_total <- atoms$population * exp(log_rate_percapita)
      atoms[[paste0("covariate_y_", i)]] <- rpois(nrow(atoms), lambda = lambda_total)
      
    } else if (dist_covariates_y[i] == 'binomial') {
      probabilities <- plogis(y_intercepts[i] + y_atom_results[,i])
      probabilities <- pmin(pmax(probabilities, 0.001), 0.999)
      atoms[[paste0("covariate_y_", i)]] <- rbinom(nrow(atoms), size = atoms$population, prob = probabilities)
      
    } else if (dist_covariates_y[i] == 'normal') {
      atoms[[paste0("covariate_y_", i)]] <- atoms$population * rnorm(nrow(atoms), mean = y_intercepts[i] + y_atom_results[,i], sd = 1)
      
    } else {
      stop(paste0('Distribution of Y-grid covariate ', i, ' not recognized'))
    }
  }
  
  # Generate atom-level outcome
  outcome_re <- gen_correlated_spat(W_atoms, 1, correlation = 1, rho = rho_y)
  linear_pred_y_atom <- rep(beta0_y, nrow(atoms))
  
  # Add X covariate effects
  for(i in 1:n_covariates_x) {
    percapita_rate_x <- atoms[[paste0("covariate_x_", i)]] / atoms$population
    linear_pred_y_atom <- linear_pred_y_atom + beta_x[i] * percapita_rate_x
  }
  
  # Add Y covariate effects
  for(i in 1:n_covariates_y) {
    percapita_rate_y <- atoms[[paste0("covariate_y_", i)]] / atoms$population
    linear_pred_y_atom <- linear_pred_y_atom + beta_y[i] * percapita_rate_y
  }
  
  # Generate final outcome
  if (dist_y == 'poisson') {
    rate_percapita_y <- exp(linear_pred_y_atom + outcome_re[,1])
    atoms$y <- rpois(nrow(atoms), lambda = atoms$population * rate_percapita_y)
  } else if (dist_y == 'binomial') {
    outcome_prob_atom <- plogis(linear_pred_y_atom + outcome_re[,1])
    atoms$y <- rbinom(nrow(atoms), size = atoms$population, prob = outcome_prob_atom)
  } else if (dist_y == 'normal') {
    atoms$y <- rnorm(nrow(atoms), 
                     mean = atoms$population * (linear_pred_y_atom + outcome_re[,1]), 
                     sd = atoms$population * 1)
  } else {
    stop('Distribution of outcome not recognized.')
  }
  
  # Aggregate to grid level
  for(i in 1:n_covariates_x) {
    gridx[[paste0("covariate_x_", i)]] <- sapply(1:nrow(gridx), function(j) {
      atom_indices <- which(atoms$ID_x == j)
      sum(atoms[[paste0("covariate_x_", i)]][atom_indices])
    })
  }
  
  for(i in 1:n_covariates_y) {
    gridy[[paste0("covariate_y_", i)]] <- sapply(1:nrow(gridy), function(j) {
      atom_indices <- which(atoms$ID_y == j)
      sum(atoms[[paste0("covariate_y_", i)]][atom_indices])
    })
  }
  
  gridy$y <- sapply(1:nrow(gridy), function(j) {
    atom_indices <- which(atoms$ID_y == j)
    sum(atoms$y[atom_indices])
  })
  
  gridy$population <- sapply(1:nrow(gridy), function(j) {
    atom_indices <- which(atoms$ID_y == j)
    sum(atoms$population[atom_indices])
  })
  
  gridx$population <- sapply(1:nrow(gridx), function(j) {
    atom_indices <- which(atoms$ID_x == j)
    sum(atoms$population[atom_indices])
  })
  
  # Store true parameters
  true_params <- list(
    beta_x = beta_x,
    beta_y = beta_y,
    x_correlation = x_correlation,
    y_correlation = y_correlation
  )
  
  result <- list(
    gridy = gridy,
    gridx = gridx,
    atoms = atoms,
    true_params = true_params
  )
  class(result) <- "misaligned_data"
  return(result)
}