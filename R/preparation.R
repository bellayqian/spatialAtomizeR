#' Prepare Spatial Bookkeeping
#'
#' @param data List containing gridy, gridx, and atoms from simulate_misaligned_data
#'
#' @return List of bookkeeping objects for NIMBLE model
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom nimble as.carAdjacency
#' @importFrom stats rmultinom
#' @importFrom methods as
prepare_spatial_bookkeeping <- function(data) {
  gridy <- data$gridy
  gridx <- data$gridx
  atoms <- data$atoms
  S_x <- nrow(data$gridx)
  S_y <- nrow(data$gridy)
  
  ## add an atom-order ID
  atoms$ID_atomorder<-1:nrow(atoms)
  
  # Ensure IDs are present
  # gridx$ID <- 1:nrow(gridx)
  # gridy$ID <- 1:nrow(gridy)
  
  # Identify non-atom grids
  ## x grids that aren't atoms
  # x_vars <- c("x", grep("covariate_x_", names(gridx), value = TRUE))
  x_vars <- c(grep("covariate_x_", names(gridx), value = TRUE))
  x_nonatoms <- unique(atoms$ID_x[which(duplicated(atoms$ID_x)==T)])
  x_nonatoms <- x_nonatoms[order(x_nonatoms)]
  J_x <- nrow(gridx) - length(x_nonatoms)
  
  ## y grids that aren't atoms
  y_vars <- c("y", grep("covariate_y_", names(gridy), value = TRUE))
  y_nonatoms <- unique(atoms$ID_y[which(duplicated(atoms$ID_y)==T)])
  y_nonatoms <- y_nonatoms[order(y_nonatoms)]
  J_y <- nrow(gridy) - length(y_nonatoms)
  
  # Re-order x and y so that first J_x/J_y units are atom-equivalents
  y_nonatom_id <- sapply(y_nonatoms,FUN=function(j) which(gridy$ID==j))
  gridy_yorder <- gridy[c((1:S_y)[-y_nonatom_id], y_nonatom_id),]
  names(gridy_yorder)[which(names(gridy_yorder)=='ID')] <- 'ID_y'
  gridy_yorder$ID_yorder <- 1:nrow(gridy_yorder)
  
  x_nonatom_id <- sapply(x_nonatoms,FUN=function(j) which(gridx$ID==j))
  gridx_xorder <- gridx[c((1:S_x)[-x_nonatom_id], x_nonatom_id),]
  names(gridx_xorder)[which(names(gridx_xorder)=='ID')] <- 'ID_x'
  gridx_xorder$ID_xorder <- 1:nrow(gridx_xorder)
  
  # Prepare ID crosswalk. Link the new IDs into the atoms dataset so we can find
  # correspondence between new Y grid ids and atoms
  IDxwalk <- merge(st_drop_geometry(atoms), st_drop_geometry(gridy_yorder), by='ID_y')
  IDxwalk <- merge(IDxwalk, st_drop_geometry(gridx_xorder), by='ID_x')
  
  # Prepare mapping vectors, map atoms to grids
  ## get x/y_to_atom
  x_to_atom <- IDxwalk$ID_atomorder[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder)]
  y_to_atom <- IDxwalk$ID_atomorder[order(IDxwalk$ID_yorder, IDxwalk$ID_atomorder)]
  
  ## get expand_x/y
  expand_x <- IDxwalk$ID_xorder[order(IDxwalk$ID_xorder)]
  expand_y <- IDxwalk$ID_yorder[order(IDxwalk$ID_yorder)]
  
  # Prepare latent indices
  ## get x/y_latentind
  IDxwalk <- IDxwalk[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder),]
  x_nonatoms <- unique(IDxwalk$ID_xorder[which(duplicated(IDxwalk$ID_xorder)==T)])
  x_latentind <- matrix(0, length(x_nonatoms), 2)
  for (i in 1:length(x_nonatoms)) {
    x_latentind[i,] <- c(min(which(IDxwalk$ID_xorder==x_nonatoms[i])), max(which(IDxwalk$ID_xorder==x_nonatoms[i]))) - J_x
  }
  
  IDxwalk <- IDxwalk[order(IDxwalk$ID_yorder, IDxwalk$ID_atomorder),]
  y_nonatoms <- unique(IDxwalk$ID_yorder[which(duplicated(IDxwalk$ID_yorder)==T)])
  y_latentind <- matrix(0, length(y_nonatoms), 2)
  for (i in 1:length(y_nonatoms)) {
    y_latentind[i,] <- c(min(which(IDxwalk$ID_yorder==y_nonatoms[i])), max(which(IDxwalk$ID_yorder==y_nonatoms[i]))) - J_y
  }
  
  # Prepare x_reorder
  IDxwalk <- IDxwalk[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder),]
  IDxwalk$xro <- 1:nrow(IDxwalk)
  IDxwalk <- IDxwalk[order(IDxwalk$ID_atomorder),]
  x_reorder <- IDxwalk$xro
  
  return(list(J_x = J_x, J_y = J_y, 
              x_to_atom = x_to_atom, y_to_atom = y_to_atom,
              expand_x = expand_x, expand_y = expand_y, 
              x_latentind = x_latentind, y_latentind = y_latentind, 
              x_reorder = x_reorder,
              gridy_yorder = gridy_yorder, gridx_xorder = gridx_xorder,
              x_vars = x_vars, y_vars = y_vars))
}

#' Prepare Adjacency Matrices
#'
#' @param gridy_yorder Reordered Y grid
#' @param gridx_xorder Reordered X grid
#'
#' @return List containing W_x and W_y adjacency matrices
#' @export
#' @importFrom spdep poly2nb nb2mat
#' @importFrom dplyr select
prepare_adjacency_matrices <- function(gridy_yorder, gridx_xorder) {
  neighbors_x <- poly2nb(as(gridx_xorder, "Spatial"), queen = TRUE)
  W_x <- nb2mat(neighbors_x, style = "B", zero.policy = TRUE)
  
  neighbors_y <- poly2nb(as(gridy_yorder, "Spatial"), queen = TRUE)
  W_y <- nb2mat(neighbors_y, style = "B", zero.policy = TRUE)
  
  return(list(W_x = W_x, W_y = W_y))
}

#' Prepare NIMBLE Model Inputs
#'
#' @param bookkeeping Output from prepare_spatial_bookkeeping
#' @param adjacency Output from prepare_adjacency_matrices
#' @param data Original simulated data
#' @param norm_idx_x Indices of normal-distributed X covariates
#' @param pois_idx_x Indices of Poisson-distributed X covariates
#' @param binom_idx_x Indices of binomial-distributed X covariates
#' @param norm_idx_y Indices of normal-distributed Y covariates
#' @param pois_idx_y Indices of Poisson-distributed Y covariates
#' @param binom_idx_y Indices of binomial-distributed Y covariates
#' @param dist_y Distribution type for outcome (1=normal, 2=poisson, 3=binomial)
#'
#' @return List containing constants, data, and inits for NIMBLE
#' @export
#' @importFrom BiasedUrn rMFNCHypergeo
#' @importFrom stats rWishart rnorm
#' @importFrom nimble as.carAdjacency
prepare_nimble_inputs <- function(bookkeeping, adjacency, data,
                                  norm_idx_x = NULL, pois_idx_x = NULL, binom_idx_x = NULL,
                                  norm_idx_y = NULL, pois_idx_y = NULL, binom_idx_y = NULL,
                                  dist_y = 2) {
  # Basic dimensions
  gridy_yorder <- bookkeeping$gridy_yorder
  gridx_xorder <- bookkeeping$gridx_xorder
  atoms <- data$atoms
  n_atoms <- if(is.null(nrow(atoms))) length(atoms) else nrow(atoms)
  
  # Extract population data
  pop_atoms <-atoms$population
  
  # Identify variables
  x_vars <- c(grep("covariate_x_", names(gridx_xorder), value = TRUE))
  y_vars <- c("y", grep("covariate_y_", names(gridy_yorder), value = TRUE))
  
  # Create covariate matrices
  create_covariate_matrix <- function(grid_data, vars, grid_size) {
    if("ID_x" %in% names(grid_data)) {
      # For X grid
      covar_matrix <- matrix(NA, nrow = grid_size, ncol = length(vars))
      colnames(covar_matrix) <- vars
      
      # Fill in values
      for(i in 1:length(vars)) {
        covar_matrix[,i] <- as.numeric(grid_data[[vars[i]]])
      }
    } else {
      # For Y grid - p_y covariates first, then y
      p_y <- length(vars) - 1  # subtract 1 because 'y' is included in vars
      covar_matrix <- matrix(NA, nrow = grid_size, ncol = p_y + 1)
      
      # Get covariate names (excluding 'y')
      covar_names <- vars[vars != "y"]
      
      # Fill in covariates first
      for(i in 1:p_y) {
        covar_matrix[,i] <- as.numeric(grid_data[[covar_names[i]]])
      }
      
      # Fill in y as the last column
      covar_matrix[,p_y + 1] <- as.numeric(grid_data[["y"]])
      colnames(covar_matrix) <- c(covar_names, "y")
    }
    
    return(covar_matrix)
  }
  
  print("\nCreating and standardizing covariate matrices...")
  covar_x <- create_covariate_matrix(grid_data=gridx_xorder[order(gridx_xorder$ID_xorder), ],
                                     vars=bookkeeping$x_vars,
                                     grid_size=nrow(gridx_xorder))
  
  covar_y <- create_covariate_matrix(grid_data=gridy_yorder[order(gridy_yorder$ID_yorder), ],
                                     vars=bookkeeping$y_vars,
                                     grid_size=nrow(gridy_yorder))
  
  # Create index matrices
  p_x <- ncol(covar_x)
  p_y <- ncol(covar_y) - 1  # subtract y column
  
  # Scale matrices for Wishart priors
  R_x <- diag(p_x)  # For X grid covariates
  R_yx <- diag(p_y) # For Y grid covariates
  
  D <- n_atoms
  S_x <- nrow(gridx_xorder)
  S_y <- nrow(gridy_yorder)
  
  # Create initial precision matrices that are guaranteed positive definite
  initialize_precision <- function(dim) {
    # Create a well-conditioned positive definite matrix
    scale <- diag(dim)
    df <- dim + 2  # Slightly more degrees of freedom for stability
    W <- rWishart(1, df, scale)[,,1]
    # Scale the matrix to have reasonable magnitudes
    W <- W / mean(diag(W))
    return(W)
  }
  
  constants <- list(
    # Dimensions
    p_x = p_x,                  # Number of X covariates
    p_y = p_y,                  # Number of Y covariates
    D = n_atoms,
    S_x = nrow(gridx_xorder),   # X grid size
    S_y = nrow(gridy_yorder),   # Y grid size
    J_x = bookkeeping$J_x,      # Number of atom-equivalent X grids
    J_y = bookkeeping$J_y,      # Number of atom-equivalent Y grids
    
    # Population data
    pop_atoms_x = atoms$population[bookkeeping$x_to_atom],      # Population for each atom
    pop_atoms_y = atoms$population[bookkeeping$y_to_atom],      # Population for each atom
    
    # Mapping between spaces
    # x_to_atom = bookkeeping$x_to_atom,  # Maps X grids to atoms
    y_to_atom = bookkeeping$y_to_atom,  # Maps Y grids to atoms
    expand_x = bookkeeping$expand_x,     # Expands atoms to X grid
    expand_y = bookkeeping$expand_y,     # Expands atoms to Y grid
    x_reorder = bookkeeping$x_reorder,
    
    # Latent indices
    xlatent_ind = bookkeeping$x_latentind,
    ylatent_ind = bookkeeping$y_latentind,
    
    # CAR structures
    num_x = as.numeric(as.carAdjacency(adjacency$W_x)$num),
    weights_x = as.numeric(as.carAdjacency(adjacency$W_x)$weights),
    adj_x = as.numeric(as.carAdjacency(adjacency$W_x)$adj),
    num_y = as.numeric(as.carAdjacency(adjacency$W_y)$num),
    weights_y = as.numeric(as.carAdjacency(adjacency$W_y)$weights),
    adj_y = as.numeric(as.carAdjacency(adjacency$W_y)$adj),
    
    # Wishart prior matrices
    R_x = R_x,         # Add scale matrix for X covariates
    R_yx = R_yx,       # Add scale matrix for Y covariates
    df_x = p_x + 2,    # Degrees of freedom for Wishart
    df_yx = p_y + 2,    # Degrees of freedom for Wishart
    
    ## hack: you have to add an extra 0 onto the end of the variable type id's to account for the dimension 1 edge case
    norm_n_x=length(norm_idx_x),
    norm_idx_x=c(norm_idx_x,0),
    pois_n_x=length(pois_idx_x),
    pois_idx_x=c(pois_idx_x,0),
    binom_n_x=length(binom_idx_x),
    binom_idx_x=c(binom_idx_x,0),
    
    norm_n_y=length(norm_idx_y),
    norm_idx_y=c(norm_idx_y,0),
    pois_n_y=length(pois_idx_y),
    pois_idx_y=c(pois_idx_y,0),
    binom_n_y=length(binom_idx_y),
    binom_idx_y=c(binom_idx_y,0),
    
    vartype_y = dist_y
  )
  
  data <- list(
    x = covar_x,     # X covariates
    yx_obs = covar_y[, 1:p_y, drop=F],  # Add explicit variables
    y_obs = covar_y[, p_y + 1]  # Y outcome
  )
  
  inits <- list(
    # Beta parameters with small non-zero values to avoid NA warnings
    beta_0_x = rep(0.01, p_x),
    beta_0_y = 0.01,
    beta_0_yx = rep(0.01, p_y),
    beta_y = rep(0.01, p_x + p_y),
    
    # Spatial random effects
    tau_x = rep(1, p_x),
    tau_y = 1,
    tau_yx = rep(1, p_y),
    
    # Matrices for spatial correlation
    # Prec_x = initialize_precision(p_x),
    # Prec_yx = initialize_precision(p_y),
    
    # Spatial random effects
    # psi_x = matrix(rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x),
    # psi_yx = matrix(rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y),
    phi_y = rnorm(S_y, 0, 0.1),
    phi_yx = matrix(rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y),
    phi_x = matrix(rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x),
    
    # Latent variables
    x_latent = matrix(1, nrow = D - constants$J_x, ncol = p_x),
    y_latent = rep(1, D - constants$J_y),
    yx_latent = matrix(1, nrow = D - constants$J_y, ncol = p_y)
    
  )
  
  if (p_x>1){
    inits$psi_x <- matrix(rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x)
    inits$Prec_x <- initialize_precision(p_x)
  }
  
  if (p_y>1){
    inits$psi_yx <- matrix(rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y)
    inits$Prec_yx <- initialize_precision(p_y)
  }
  
  if (length(norm_idx_x)>0){
    inits$sigma2_x<-rep(1,length(norm_idx_x))
    inits$w_x<-matrix(rnorm(length(norm_idx_x)*(n_atoms-bookkeeping$J_x)),nrow=(n_atoms-bookkeeping$J_x),ncol=length(norm_idx_x))
  }
  
  if (length(norm_idx_y)>0){
    inits$sigma2_yx<-rep(1,length(norm_idx_y))
    inits$w_yx<-matrix(rnorm(length(norm_idx_y)*(n_atoms-bookkeeping$J_y)),nrow=(n_atoms-bookkeeping$J_y),ncol=length(norm_idx_y))
  }
  
  if (dist_y==1){
    inits$sigma2_y<-1
    inits$w_y<-rnorm((n_atoms-bookkeeping$J_y))
  }
  
  # Add this diagnostic function to check the indexing
  check_nimble_indexing <- function(constants, data) {
    message("=== NIMBLE Indexing Diagnostics ===\n")
    
    D <- constants$D
    J_y <- constants$J_y
    
    # Check if y_to_atom exists
    if(!"y_to_atom" %in% names(constants)) {
      message("ERROR: y_to_atom missing from constants!\n")
      return(FALSE)
    }
    
    # Check dimensions
    message("D (atoms):", D, "\n")
    message("J_y (atom-equivalent Y grids):", J_y, "\n")
    message("Length of y_to_atom:", length(constants$y_to_atom), "\n")
    message("Length of expand_y:", length(constants$expand_y), "\n")
    message("Length of pop_atoms:", length(constants$pop_atoms), "\n")
    
    # Check population alignment
    # Check y_to_atom mapping
    if(any(constants$y_to_atom > D) || any(constants$y_to_atom < 1)) {
      message("Warning: y_to_atom contains invalid indices!\n")
      return(FALSE)
    }
    
    # Check temp_yx dimensions
    if(nrow(data$yx_obs) != constants$S_y) {
      message("Warning: yx_obs rows don't match S_y!\n")
      return(FALSE)
    }
    
    message("Basic indexing checks passed\n")
    return(TRUE)
  }
  # check_nimble_indexing(constants, data)
  
  # Initialize latent values with proper population-scaled splits
  # For x_latent
  for(j in 1:p_x) {
    for(m in 1:(S_x-constants$J_x)) {
      total <- data$x[constants$J_x + m, j]
      n_atoms <- diff(constants$xlatent_ind[m,]) + 1
      pops_atoms<-constants$pop_atoms_x[constants$J_x+(constants$xlatent_ind[m,1]:constants$xlatent_ind[m,2])]
      
      if (j %in% constants$norm_idx_x){
        ## if variable is normal, divide grid total evenly between atoms
        inits$x_latent[constants$xlatent_ind[m,1]:constants$xlatent_ind[m,2], j] <- total/n_atoms
      } else{
        ## if variable is count, divide grid total using multinomial
        inits$x_latent[constants$xlatent_ind[m,1]:constants$xlatent_ind[m,2], j] <- rmultinom(n=1,size=total,prob=pops_atoms/sum(pops_atoms))
      }
    }
  }
  
  # For y_latent
  for(m in 1:(S_y-constants$J_y)) {
    total<- data$y_obs[constants$J_y + m]
    n_atoms <- diff(constants$ylatent_ind[m,]) + 1
    pops_atoms<-constants$pop_atoms_y[constants$J_y+(constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2])]
    
    if (dist_y==1){
      inits$y_latent[constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2]] <- total/n_atoms
    } else{
      inits$y_latent[constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2]] <- rmultinom(n=1,size=total,prob=pops_atoms/sum(pops_atoms))
    }
  }
  
  # For yx_latent
  for(j in 1:p_y) {
    for(m in 1:(S_y-constants$J_y)) {
      total <- data$yx_obs[constants$J_y + m, j]
      n_atoms <- diff(constants$ylatent_ind[m,]) + 1
      pops_atoms<-constants$pop_atoms_y[constants$J_y+(constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2])]
      
      if (j %in% constants$norm_idx_y){
        inits$yx_latent[constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2], j] <- total/n_atoms
      } else{
        inits$yx_latent[constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2], j] <- rmultinom(n=1,size=total,prob=pops_atoms/sum(pops_atoms))
      }
    }
  }
  
  
  return(list(
    constants = constants,
    data = data,
    inits = inits
  ))
}




