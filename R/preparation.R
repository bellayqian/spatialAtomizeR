# The functions in this file operate on *standardised* inputs built by
# run_abrm(): gridx has an integer `ID` (1..S_x) and covariates
# `covariate_x_1..p_x`; gridy has an integer `ID` (1..S_y), covariates
# `covariate_y_1..p_y` and outcome `y`; atoms has integer `ID_x`, `ID_y`,
# population column(s) and optional atom covariates `covariate_a_1..p_a`.
# Users never need these names: run_abrm() maps their own column names onto
# them and maps results back.

#' Prepare Spatial Bookkeeping
#'
#' @param data List containing standardised gridy, gridx, and atoms
#'
#' @return List of bookkeeping objects for NIMBLE model
#' @keywords internal
#' @noRd
#' @importFrom sf st_drop_geometry
prepare_spatial_bookkeeping <- function(data) {
  gridy <- data$gridy
  gridx <- data$gridx
  atoms <- data$atoms
  S_x <- nrow(gridx)
  S_y <- nrow(gridy)

  atoms$ID_atomorder <- seq_len(nrow(atoms))

  x_vars <- grep("^covariate_x_", names(gridx), value = TRUE)
  y_vars <- c("y", grep("^covariate_y_", names(gridy), value = TRUE))

  # Grid cells containing more than one atom are "non-atoms": their values
  # must be split across atoms as latent variables.
  x_nonatoms <- sort(unique(atoms$ID_x[duplicated(atoms$ID_x)]))
  y_nonatoms <- sort(unique(atoms$ID_y[duplicated(atoms$ID_y)]))
  J_x <- S_x - length(x_nonatoms)
  J_y <- S_y - length(y_nonatoms)

  # Re-order grids so that the first J_x / J_y cells are atom-equivalents
  y_nonatom_id <- match(y_nonatoms, gridy$ID)
  gridy_yorder <- gridy[c(setdiff(seq_len(S_y), y_nonatom_id), y_nonatom_id), ]
  names(gridy_yorder)[names(gridy_yorder) == "ID"] <- "ID_y"
  gridy_yorder$ID_yorder <- seq_len(nrow(gridy_yorder))

  x_nonatom_id <- match(x_nonatoms, gridx$ID)
  gridx_xorder <- gridx[c(setdiff(seq_len(S_x), x_nonatom_id), x_nonatom_id), ]
  names(gridx_xorder)[names(gridx_xorder) == "ID"] <- "ID_x"
  gridx_xorder$ID_xorder <- seq_len(nrow(gridx_xorder))

  # ID crosswalk: link the re-ordered grid IDs to every atom
  id_cols_atoms <- sf::st_drop_geometry(atoms)[, c("ID_x", "ID_y", "ID_atomorder")]
  IDxwalk <- merge(id_cols_atoms,
                   sf::st_drop_geometry(gridy_yorder)[, c("ID_y", "ID_yorder")],
                   by = "ID_y")
  IDxwalk <- merge(IDxwalk,
                   sf::st_drop_geometry(gridx_xorder)[, c("ID_x", "ID_xorder")],
                   by = "ID_x")

  # Map atoms to grids
  x_to_atom <- IDxwalk$ID_atomorder[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder)]
  y_to_atom <- IDxwalk$ID_atomorder[order(IDxwalk$ID_yorder, IDxwalk$ID_atomorder)]

  expand_x <- IDxwalk$ID_xorder[order(IDxwalk$ID_xorder)]
  expand_y <- IDxwalk$ID_yorder[order(IDxwalk$ID_yorder)]

  # Latent index ranges: for each non-atom cell, the [first, last] positions
  # of its atoms among the latent (post-J) atoms.
  latent_ranges <- function(order_id, J) {
    nonatoms <- unique(order_id[duplicated(order_id)])
    out <- matrix(0L, length(nonatoms), 2L)
    for (i in seq_along(nonatoms)) {
      pos <- which(order_id == nonatoms[i])
      out[i, ] <- c(min(pos), max(pos)) - J
    }
    storage.mode(out) <- "integer"
    out
  }

  IDxwalk <- IDxwalk[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder), ]
  x_latentind <- latent_ranges(IDxwalk$ID_xorder, J_x)
  IDxwalk$xro <- seq_len(nrow(IDxwalk))

  IDxwalk <- IDxwalk[order(IDxwalk$ID_yorder, IDxwalk$ID_atomorder), ]
  y_latentind <- latent_ranges(IDxwalk$ID_yorder, J_y)

  # x_reorder: position (in X order) of each atom, listed in atom order
  x_reorder <- IDxwalk$xro[order(IDxwalk$ID_atomorder)]

  list(J_x = as.integer(J_x), J_y = as.integer(J_y),
       x_to_atom = as.integer(x_to_atom), y_to_atom = as.integer(y_to_atom),
       expand_x = as.integer(expand_x), expand_y = as.integer(expand_y),
       x_latentind = x_latentind, y_latentind = y_latentind,
       x_reorder = as.integer(x_reorder),
       gridy_yorder = gridy_yorder, gridx_xorder = gridx_xorder,
       x_vars = x_vars, y_vars = y_vars)
}

#' Prepare Adjacency Matrices
#'
#' @param gridy_yorder Reordered Y grid (sf)
#' @param gridx_xorder Reordered X grid (sf)
#'
#' @return List containing W_x and W_y binary adjacency matrices
#' @keywords internal
#' @noRd
#' @importFrom spdep poly2nb nb2mat
prepare_adjacency_matrices <- function(gridy_yorder, gridx_xorder) {
  neighbors_x <- spdep::poly2nb(gridx_xorder, queen = TRUE)
  W_x <- spdep::nb2mat(neighbors_x, style = "B", zero.policy = TRUE)

  neighbors_y <- spdep::poly2nb(gridy_yorder, queen = TRUE)
  W_y <- spdep::nb2mat(neighbors_y, style = "B", zero.policy = TRUE)

  list(W_x = W_x, W_y = W_y)
}

#' Prepare NIMBLE Model Inputs
#'
#' @param bookkeeping Output from prepare_spatial_bookkeeping
#' @param adjacency Output from prepare_adjacency_matrices
#' @param data Standardised data list (gridx, gridy, atoms)
#' @param x_types Character vector (length p_x) with elements in
#'   "norm", "pois", "binom"
#' @param y_types Character vector (length p_y) with elements in
#'   "norm", "pois", "binom"
#' @param dist_y Outcome distribution: "norm", "pois" or "binom"
#'
#' @return List containing constants, data, and inits for NIMBLE
#' @keywords internal
#' @noRd
#' @importFrom stats rWishart rnorm rmultinom
#' @importFrom nimble as.carAdjacency
prepare_nimble_inputs <- function(bookkeeping, adjacency, data,
                                  x_types, y_types, dist_y = "pois") {
  gridy_yorder <- bookkeeping$gridy_yorder
  gridx_xorder <- bookkeeping$gridx_xorder
  atoms <- data$atoms
  n_atoms <- nrow(atoms)

  a_vars <- grep("^covariate_a_", names(atoms), value = TRUE)

  # Covariate matrices in model (re-ordered) grid order
  gx <- gridx_xorder[order(gridx_xorder$ID_xorder), ]
  gy <- gridy_yorder[order(gridy_yorder$ID_yorder), ]
  covar_x <- vapply(bookkeeping$x_vars, function(v) as.numeric(gx[[v]]),
                    numeric(nrow(gx)))
  covar_x <- matrix(covar_x, nrow = nrow(gx), dimnames = list(NULL, bookkeeping$x_vars))
  y_covar_names <- setdiff(bookkeeping$y_vars, "y")
  covar_yx <- vapply(y_covar_names, function(v) as.numeric(gy[[v]]), numeric(nrow(gy)))
  covar_yx <- matrix(covar_yx, nrow = nrow(gy), dimnames = list(NULL, y_covar_names))
  y_obs <- as.numeric(gy[["y"]])

  p_x <- ncol(covar_x)
  p_y <- ncol(covar_yx)
  p_a <- length(a_vars)

  atoms_df <- sf::st_drop_geometry(atoms)
  if (p_a > 0L) {
    covar_matrix_a <- as.matrix(atoms_df[, a_vars, drop = FALSE])
    storage.mode(covar_matrix_a) <- "double"
  }

  # Population matrices (D x p_x, D x p_y) and outcome vector (D)
  if ("population" %in% names(atoms_df)) {
    pop_vec <- atoms_df$population
    pop_matrix_x <- matrix(pop_vec, nrow = n_atoms, ncol = p_x)
    pop_matrix_y <- matrix(pop_vec, nrow = n_atoms, ncol = p_y)
    pop_outcome  <- pop_vec
  } else {
    pop_matrix_x <- as.matrix(atoms_df[, paste0("population_covariate_x_", seq_len(p_x)),
                                       drop = FALSE])
    pop_matrix_y <- as.matrix(atoms_df[, paste0("population_covariate_y_", seq_len(p_y)),
                                       drop = FALSE])
    pop_outcome  <- atoms_df$population_y
  }

  S_x <- nrow(gridx_xorder)
  S_y <- nrow(gridy_yorder)
  J_x <- bookkeeping$J_x
  J_y <- bookkeeping$J_y

  # Well-conditioned positive-definite initial precision matrix
  initialize_precision <- function(dim) {
    W <- stats::rWishart(1, dim + 2, diag(dim))[, , 1]
    W / mean(diag(W))
  }

  # Integer indices of each variable type. A trailing 0 is appended so that
  # NIMBLE always sees a vector (avoids the length-1 edge case).
  type_idx <- function(types, t) c(which(types == t), 0L)
  car_x <- nimble::as.carAdjacency(adjacency$W_x)
  car_y <- nimble::as.carAdjacency(adjacency$W_y)

  constants <- list(
    # Dimensions
    p_x = p_x, p_y = p_y, p_a = p_a,
    D = n_atoms,
    S_x = S_x, S_y = S_y,
    J_x = J_x, J_y = J_y,

    # Populations, re-ordered to model atom order
    pop_atoms_x = pop_matrix_x[bookkeeping$x_to_atom, , drop = FALSE],
    pop_atoms_y = pop_matrix_y[bookkeeping$y_to_atom, , drop = FALSE],
    pop_atoms_outcome = pop_outcome[bookkeeping$y_to_atom],

    # Mapping between spaces
    y_to_atom = bookkeeping$y_to_atom,
    expand_x = bookkeeping$expand_x,
    expand_y = bookkeeping$expand_y,
    x_reorder = bookkeeping$x_reorder,

    # Latent index ranges
    xlatent_ind = bookkeeping$x_latentind,
    ylatent_ind = bookkeeping$y_latentind,

    # CAR structures (num/adj are integer indices, weights are numeric)
    num_x = as.integer(car_x$num),
    weights_x = as.numeric(car_x$weights),
    adj_x = as.integer(car_x$adj),
    num_y = as.integer(car_y$num),
    weights_y = as.numeric(car_y$weights),
    adj_y = as.integer(car_y$adj),

    # Wishart priors
    R_x = diag(p_x), R_yx = diag(p_y),
    df_x = p_x + 2L, df_yx = p_y + 2L,

    # Variable-type indices
    norm_n_x = sum(x_types == "norm"),   norm_idx_x = type_idx(x_types, "norm"),
    pois_n_x = sum(x_types == "pois"),   pois_idx_x = type_idx(x_types, "pois"),
    binom_n_x = sum(x_types == "binom"), binom_idx_x = type_idx(x_types, "binom"),
    norm_n_y = sum(y_types == "norm"),   norm_idx_y = type_idx(y_types, "norm"),
    pois_n_y = sum(y_types == "pois"),   pois_idx_y = type_idx(y_types, "pois"),
    binom_n_y = sum(y_types == "binom"), binom_idx_y = type_idx(y_types, "binom"),

    vartype_y = unname(.dist_code[dist_y])
  )

  if (p_a > 0L) {
    constants$covars_atom_x <- covar_matrix_a[bookkeeping$x_to_atom, , drop = FALSE]
    constants$covars_atom_y <- covar_matrix_a[bookkeeping$y_to_atom, , drop = FALSE]
  }

  nimble_data <- list(
    x = covar_x,
    yx_obs = covar_yx,
    y_obs = y_obs
  )

  n_norm_x <- constants$norm_n_x
  n_norm_y <- constants$norm_n_y

  inits <- list(
    beta_0_x = rep(0.01, p_x),
    beta_0_y = 0.01,
    beta_0_yx = rep(0.01, p_y),
    tau_x = rep(1, p_x),
    tau_y = 1,
    tau_yx = rep(1, p_y),
    phi_y = stats::rnorm(S_y, 0, 0.1),
    phi_yx = matrix(stats::rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y),
    phi_x = matrix(stats::rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x),
    x_latent = matrix(1, nrow = n_atoms - J_x, ncol = p_x),
    y_latent = rep(1, n_atoms - J_y),
    yx_latent = matrix(1, nrow = n_atoms - J_y, ncol = p_y),
    beta_y = rep(0.01, p_x + p_y + p_a)
  )

  if (p_x > 1L) {
    inits$psi_x <- matrix(stats::rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x)
    inits$Prec_x <- initialize_precision(p_x)
  }
  if (p_y > 1L) {
    inits$psi_yx <- matrix(stats::rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y)
    inits$Prec_yx <- initialize_precision(p_y)
  }
  if (n_norm_x > 0L) {
    inits$sigma2_x <- rep(1, n_norm_x)
    inits$w_x <- matrix(stats::rnorm(n_norm_x * (n_atoms - J_x)),
                        nrow = n_atoms - J_x, ncol = n_norm_x)
  }
  if (n_norm_y > 0L) {
    inits$sigma2_yx <- rep(1, n_norm_y)
    inits$w_yx <- matrix(stats::rnorm(n_norm_y * (n_atoms - J_y)),
                         nrow = n_atoms - J_y, ncol = n_norm_y)
  }
  if (dist_y == "norm") {
    inits$sigma2_y <- 1
    inits$w_y <- stats::rnorm(n_atoms - J_y)
  }
  if (p_a > 0L) {
    inits$beta_a_x <- matrix(0.01, nrow = p_a, ncol = p_x)
    inits$beta_a_yx <- matrix(0.01, nrow = p_a, ncol = p_y)
  }

  # Initialise latent atom values with population-proportional splits of the
  # observed cell totals (equal splits for normal variables).
  split_total <- function(total, pops, is_normal) {
    if (is_normal) rep(total / length(pops), length(pops))
    else stats::rmultinom(n = 1, size = total, prob = pops / sum(pops))
  }

  for (j in seq_len(p_x)) {
    for (m in seq_len(S_x - J_x)) {
      rng <- constants$xlatent_ind[m, 1]:constants$xlatent_ind[m, 2]
      inits$x_latent[rng, j] <- split_total(nimble_data$x[J_x + m, j],
                                            constants$pop_atoms_x[J_x + rng, j],
                                            x_types[j] == "norm")
    }
  }
  for (m in seq_len(S_y - J_y)) {
    rng <- constants$ylatent_ind[m, 1]:constants$ylatent_ind[m, 2]
    inits$y_latent[rng] <- split_total(nimble_data$y_obs[J_y + m],
                                       constants$pop_atoms_outcome[J_y + rng],
                                       dist_y == "norm")
  }
  for (j in seq_len(p_y)) {
    for (m in seq_len(S_y - J_y)) {
      rng <- constants$ylatent_ind[m, 1]:constants$ylatent_ind[m, 2]
      inits$yx_latent[rng, j] <- split_total(nimble_data$yx_obs[J_y + m, j],
                                             constants$pop_atoms_y[J_y + rng, j],
                                             y_types[j] == "norm")
    }
  }

  list(constants = constants, data = nimble_data, inits = inits)
}
