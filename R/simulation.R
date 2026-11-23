#' Generate Correlated Spatial Effects
#'
#' Draws one realisation of a multivariate proper conditional autoregressive
#' (MCAR) random field over the areal units described by the adjacency matrix
#' \code{W}. Each column of the result is one spatially structured variable,
#' and the columns are correlated with one another according to
#' \code{correlation}.
#'
#' @details
#' Let \eqn{Q = D_W - \rho W}, where \eqn{D_W} is the diagonal matrix of
#' neighbour counts, and let \eqn{C} be the \code{n_vars} x \code{n_vars}
#' between-variable correlation matrix. The generated matrix \eqn{\Phi} has
#' \deqn{\mathrm{Cov}(\Phi_{\cdot j}, \Phi_{\cdot k}) = \sigma^2 C_{jk} Q^{-1},}
#' with \eqn{\sigma^2 = } \code{var_spat}; that is,
#' \eqn{\mathrm{vec}(\Phi) \sim N(0, \sigma^2 \, C \otimes Q^{-1})}.
#'
#' The draw is computed as \eqn{\Phi = \sigma R_Q^{-1} Z U_C}, where
#' \eqn{Q = R_Q^\top R_Q} and \eqn{C = U_C^\top U_C} are Cholesky
#' factorisations and \eqn{Z} is an \eqn{n \times p} matrix of independent
#' standard normals. This avoids forming and inverting the \eqn{np \times np}
#' Kronecker-product precision matrix, so the cost is \eqn{O(n^3 + p^3)}
#' rather than \eqn{O(n^3 p^3)}, and no dense joint covariance is ever stored.
#'
#' \eqn{Q} is positive definite when \eqn{|\rho| < 1} and every unit has at
#' least one neighbour; both conditions are checked, with informative errors.
#'
#' @param W Symmetric, non-negative \eqn{n \times n} adjacency matrix (for
#'   example from \code{spdep::nb2mat(..., style = "B")}). Every unit must
#'   have at least one neighbour.
#' @param n_vars Number of variables (columns) to generate.
#' @param rho Spatial dependence parameter, \eqn{-1 < \rho < 1}. Values closer
#'   to 1 give stronger spatial smoothing (default 0.6).
#' @param var_spat Positive scale (variance multiplier) applied to the whole
#'   field (default 1).
#' @param correlation Between-variable correlation. Either a single number,
#'   used for every pair of variables (exchangeable correlation), or a full
#'   \code{n_vars} x \code{n_vars} correlation matrix. Ignored when
#'   \code{n_vars = 1}. Default 0.5.
#' @param verify Logical. If \code{TRUE}, a message compares the target
#'   between-variable correlation matrix with the empirical correlation of the
#'   generated columns (default \code{FALSE}).
#' @return Numeric matrix with \code{nrow(W)} rows and \code{n_vars} columns,
#'   where each column is one spatially correlated variable.
#'
#' @examples
#' # A simple 4-node chain graph adjacency matrix
#' W <- matrix(c(0, 1, 0, 0,
#'               1, 0, 1, 0,
#'               0, 1, 0, 1,
#'               0, 0, 1, 0), nrow = 4, byrow = TRUE)
#'
#' # Generate 2 correlated spatial effects over the 4 areas
#' set.seed(1)
#' effects <- gen_correlated_spat(W, n_vars = 2, rho = 0.5, correlation = 0.4)
#' dim(effects)        # 4 rows (areas) x 2 columns (variables)
#'
#' # A full correlation matrix can be given instead of a single value
#' C <- matrix(c(1, 0.3, -0.2,
#'               0.3, 1, 0.5,
#'               -0.2, 0.5, 1), 3, 3)
#' effects3 <- gen_correlated_spat(W, n_vars = 3, correlation = C)
#'
#' @importFrom stats rnorm cor
#' @export
#' @keywords internal
gen_correlated_spat <- function(W, n_vars, rho = 0.6, var_spat = 1,
                                correlation = 0.5, verify = FALSE) {
  # ---- validate inputs ----------------------------------------------------
  W <- as.matrix(W)
  n <- nrow(W)
  if (!is.numeric(W) || n < 2L || ncol(W) != n) {
    stop("'W' must be a numeric square adjacency matrix with at least 2 rows.",
         call. = FALSE)
  }
  if (anyNA(W) || any(W < 0)) {
    stop("'W' must not contain missing or negative values.", call. = FALSE)
  }
  if (!isSymmetric(unname(W))) {
    stop("'W' must be symmetric.", call. = FALSE)
  }
  if (length(n_vars) != 1L || is.na(n_vars) || n_vars < 1 ||
      n_vars != round(n_vars)) {
    stop("'n_vars' must be a single positive integer.", call. = FALSE)
  }
  n_vars <- as.integer(n_vars)
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || abs(rho) >= 1) {
    stop("'rho' must be a single number strictly between -1 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(var_spat) || length(var_spat) != 1L || is.na(var_spat) ||
      var_spat <= 0) {
    stop("'var_spat' must be a single positive number.", call. = FALSE)
  }
  n_nb <- rowSums(W)
  if (any(n_nb == 0)) {
    iso <- which(n_nb == 0)
    stop(sprintf(paste0(
      "%d unit(s) in 'W' have no neighbours (first: row %d), which makes the ",
      "CAR precision matrix singular. Remove isolated units or give them ",
      "at least one neighbour."), length(iso), iso[1]), call. = FALSE)
  }

  # ---- between-variable correlation --------------------------------------
  if (n_vars == 1L) {
    C <- matrix(1, 1L, 1L)
  } else if (length(correlation) == 1L) {
    if (!is.numeric(correlation) || is.na(correlation)) {
      stop("'correlation' must be numeric.", call. = FALSE)
    }
    C <- matrix(correlation, n_vars, n_vars)
    diag(C) <- 1
  } else {
    C <- as.matrix(correlation)
    if (!is.numeric(C) || !identical(dim(C), c(n_vars, n_vars))) {
      stop(sprintf("'correlation' must be a single number or a %d x %d matrix.",
                   n_vars, n_vars), call. = FALSE)
    }
    if (!isSymmetric(unname(C)) || any(abs(diag(C) - 1) > 1e-8)) {
      stop("A 'correlation' matrix must be symmetric with a unit diagonal.",
           call. = FALSE)
    }
  }
  U_C <- tryCatch(chol(C), error = function(e) {
    msg <- "'correlation' does not define a positive-definite correlation matrix."
    if (length(correlation) == 1L && n_vars > 1L) {
      msg <- sprintf("%s A common correlation for %d variables must lie in (%s, 1).",
                     msg, n_vars, format(round(-1 / (n_vars - 1), 4)))
    }
    stop(msg, call. = FALSE)
  })

  # ---- draw the field -----------------------------------------------------
  Q   <- diag(n_nb, n) - rho * W
  R_Q <- chol(Q)                                   # Q = t(R_Q) %*% R_Q
  Z   <- matrix(stats::rnorm(n * n_vars), n, n_vars)
  # Each column of backsolve(R_Q, Z) ~ N(0, Q^{-1}); right-multiplying by U_C
  # induces between-column covariance t(U_C) %*% U_C = C.
  spatial_effects <- sqrt(var_spat) * (backsolve(R_Q, Z) %*% U_C)
  dimnames(spatial_effects) <- NULL

  if (isTRUE(verify) && n_vars > 1L) {
    emp <- stats::cor(spatial_effects)
    message("gen_correlated_spat(): target vs empirical between-variable correlation\n",
            paste(utils::capture.output(print(round(C, 3))), collapse = "\n"), "\n",
            paste(utils::capture.output(print(round(emp, 3))), collapse = "\n"),
            "\n(A single realisation of a spatially smooth field is noisy; ",
            "agreement improves with more areal units.)")
  }

  spatial_effects
}


#' Create Atoms by Intersecting Two Grids
#'
#' Computes the atoms, i.e. the non-empty polygonal pieces of the spatial
#' intersection of the Y-grid and X-grid, using
#' \code{\link[sf]{st_intersection}}. Lower-dimensional pieces where cell
#' boundaries merely touch (points and line segments) are dropped, as are
#' pieces whose area does not exceed \code{min_area}.
#'
#' @param gridy Y-grid \code{sf} data frame.
#' @param gridx X-grid \code{sf} data frame.
#' @param id_y,id_x Names of the ID columns in \code{gridy} and \code{gridx}.
#'   IDs may be integer, numeric or character, but must be unique within each
#'   grid. Whole-number numeric IDs (for example GEOIDs converted with
#'   \code{as.numeric()}) are stored as integers in the result when they fit;
#'   character IDs are kept as they are. Default \code{"ID"} for both.
#' @param min_area Atoms with area less than or equal to this value (in the
#'   units of \code{sf::st_area()}) are dropped. Useful for removing slivers
#'   caused by slightly mismatched boundaries in real shapefiles. Default 0.
#' @return An \code{sf} data frame with one row per atom and columns
#'   \code{ID_y} and \code{ID_x} holding the IDs (as in \code{gridy} and
#'   \code{gridx}) of the cells each atom belongs to.
#'
#' @examples
#' sim_data <- simulate_misaligned_data(
#'   seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'   dist_y = "norm", beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#' )
#' atoms <- create_atoms(sim_data$gridy, sim_data$gridx)
#' nrow(atoms)
#'
#' @importFrom sf st_sf st_geometry st_intersection st_dimension st_is
#'   st_collection_extract st_cast st_area st_crs
#' @export
create_atoms <- function(gridy, gridx, id_y = "ID", id_x = "ID", min_area = 0) {
  if (!inherits(gridy, "sf") || !inherits(gridx, "sf")) {
    stop("'gridy' and 'gridx' must both be sf data frames.", call. = FALSE)
  }
  if (!id_y %in% names(gridy)) {
    stop(sprintf("'gridy' has no ID column named '%s'.", id_y), call. = FALSE)
  }
  if (!id_x %in% names(gridx)) {
    stop(sprintf("'gridx' has no ID column named '%s'.", id_x), call. = FALSE)
  }
  if (anyDuplicated(gridy[[id_y]]) || anyNA(gridy[[id_y]])) {
    stop(sprintf("IDs in gridy$%s must be unique and non-missing.", id_y), call. = FALSE)
  }
  if (anyDuplicated(gridx[[id_x]]) || anyNA(gridx[[id_x]])) {
    stop(sprintf("IDs in gridx$%s must be unique and non-missing.", id_x), call. = FALSE)
  }
  if (sf::st_crs(gridy) != sf::st_crs(gridx)) {
    stop("'gridy' and 'gridx' must have the same CRS; use sf::st_transform() first.",
         call. = FALSE)
  }

  # Keep only the ID column of each grid, named so the intersection carries
  # unambiguous ID_y / ID_x columns (avoids sf's automatic "ID.1" suffix).
  gy <- sf::st_sf(ID_y = .as_integer_id(gridy[[id_y]]), geometry = sf::st_geometry(gridy))
  gx <- sf::st_sf(ID_x = .as_integer_id(gridx[[id_x]]), geometry = sf::st_geometry(gridx))
  sf::st_agr(gy) <- "constant"
  sf::st_agr(gx) <- "constant"

  atoms <- sf::st_intersection(gy, gx)

  # Cells that share an edge or a corner intersect in a LINESTRING or POINT;
  # those are not atoms. Keep polygonal pieces only.
  dims  <- sf::st_dimension(atoms)
  atoms <- atoms[!is.na(dims) & dims == 2L, ]
  is_gc <- sf::st_is(atoms, "GEOMETRYCOLLECTION")
  if (any(is_gc)) {
    poly_part <- sf::st_collection_extract(atoms[is_gc, ], "POLYGON")
    atoms <- rbind(atoms[!is_gc, ], poly_part)
  }
  if (!all(sf::st_is(atoms, "POLYGON"))) {
    atoms <- sf::st_cast(atoms, "MULTIPOLYGON")
  }

  area  <- as.numeric(sf::st_area(atoms))
  atoms <- atoms[area > min_area, ]

  atoms <- atoms[order(match(atoms$ID_y, gy$ID_y),
                       match(atoms$ID_x, gx$ID_x)), ]
  rownames(atoms) <- NULL
  atoms
}


#' Simulate Misaligned Spatial Data
#'
#' @param seed Random seed (default = 2)
#' @param dist_covariates_x Character vector giving the distribution of each
#'   synthetic X-grid covariate: \code{"norm"}, \code{"pois"} or
#'   \code{"binom"} (the long forms \code{"normal"}, \code{"poisson"},
#'   \code{"binomial"} are also accepted).
#' @param dist_covariates_y As \code{dist_covariates_x}, for the Y-grid covariates.
#' @param dist_y Distribution of the synthetic outcome: one of \code{"norm"},
#'   \code{"pois"} or \code{"binom"} (long forms also accepted).
#' @param x_intercepts Intercepts for X covariates (one per covariate).
#' @param y_intercepts Intercepts for Y covariates (one per covariate).
#' @param rho_x Spatial correlation parameter for X-grid covariates (0 to 1 with higher values yielding more spatial correlation, default = 0.6)
#' @param rho_y Spatial correlation parameter for Y-grid covariates and outcome (0 to 1 with higher values yielding more spatial correlation, default = 0.6)
#' @param x_correlation Between-variable correlation for all pairs of X-grid covariates (default = 0.5)
#' @param y_correlation Between-variable correlation for all pairs of Y-grid covariates (default = 0.5)
#' @param beta0_y Intercept for outcome model (required).
#' @param beta_x Outcome model coefficients for X-grid covariates (required; one per X covariate).
#' @param beta_y Outcome model coefficients for Y-grid covariates (required; one per Y covariate).
#' @param diff_pops Logical, indicating whether the atoms should be generated with different population sizes (diff_pops = TRUE) or a common population size (diff_pops = FALSE). Ignored when \code{diff_pops_per_var = TRUE}.
#' @param diff_pops_per_var Logical. When \code{TRUE}, each covariate and the outcome receive their own population denominator column in the atoms object (named \code{population_covariate_x_1}, \code{population_covariate_y_1}, ..., and \code{population_y}). When \code{FALSE} (default), a single \code{population} column is used for all variables.
#' @param xy_cov_cor Logical, indicating whether the atom-level spatial random effects for X-grid and Y-grid covariates should be correlated (xy_cov_cor = TRUE) or not. When set to TRUE, the x_correlation and rho_x parameters are used to generate all covariates (separate correlation parameters are not allowed for X-grid and Y-grid covariates).
#'
#' @return An object of class \code{"misaligned_data"}: a named list with
#'   components \code{gridy} (Y-grid sf data frame with integer \code{ID},
#'   outcome \code{y} and Y covariates \code{covariate_y_1}, ...),
#'   \code{gridx} (X-grid sf data frame with integer \code{ID} and X
#'   covariates \code{covariate_x_1}, ...), \code{atoms} (atom-level sf data
#'   frame, the spatial intersection of the two grids, with integer
#'   \code{ID_y}, \code{ID_x} and \code{ID_atomorder}), and
#'   \code{true_params} (an \code{\link{abrm_params}} object holding the true
#'   regression coefficients used for data generation, ready to pass to
#'   \code{\link{run_abrm}}).
#'
#' @examples
#' sim_data <- simulate_misaligned_data(
#'   seed = 1,
#'   dist_covariates_x = c("norm", "pois", "binom"),
#'   dist_covariates_y = c("norm", "pois", "binom"),
#'   dist_y = "pois",
#'   x_intercepts = c(4, -1, -1),
#'   y_intercepts = c(4, -1, -1),
#'   beta0_y = -1,
#'   x_correlation = 0.5,
#'   y_correlation = 0.5,
#'   beta_x = c(-0.03, 0.1, -0.2),
#'   beta_y = c(0.03, -0.1, 0.2)
#' )
#' class(sim_data)              # "misaligned_data"
#' print(sim_data)              # grid dimensions and atom count
#' summary(sim_data)            # includes true beta values
#' names(sim_data$atoms)        # atom-level variables
#' sim_data$true_params$beta_x  # true X-grid coefficients
#'
#' @importFrom sp GridTopology SpatialGrid
#' @importFrom sf st_as_sf st_union st_sf st_contains
#' @importFrom spdep poly2nb nb2mat
#' @importFrom stats plogis rpois rbinom rnorm
#' @importFrom methods as
#' @export
simulate_misaligned_data <- function(seed = 2,
                                     dist_covariates_x = c("norm", "pois", "binom"),
                                     dist_covariates_y = c("norm", "pois", "binom"),
                                     dist_y = "pois",
                                     x_intercepts = rep(0, length(dist_covariates_x)),
                                     y_intercepts = rep(0, length(dist_covariates_y)),
                                     rho_x = 0.6,
                                     rho_y = 0.6,
                                     x_correlation = 0.5,
                                     y_correlation = 0.5,
                                     beta0_y = NULL,
                                     beta_x = NULL,
                                     beta_y = NULL,
                                     diff_pops = TRUE,
                                     diff_pops_per_var = FALSE,
                                     xy_cov_cor = FALSE) {
  # ---- validate inputs ----------------------------------------------------
  dist_covariates_x <- .normalize_dist(dist_covariates_x, "dist_covariates_x")
  dist_covariates_y <- .normalize_dist(dist_covariates_y, "dist_covariates_y")
  dist_y <- .normalize_dist(dist_y, "dist_y")
  if (length(dist_y) != 1L) stop("'dist_y' must be a single value.", call. = FALSE)

  n_covariates_x <- length(dist_covariates_x)
  n_covariates_y <- length(dist_covariates_y)

  check_len <- function(v, n, name, what) {
    if (is.null(v)) stop(sprintf("'%s' must be supplied.", name), call. = FALSE)
    if (!is.numeric(v) || length(v) != n) {
      stop(sprintf("'%s' must be a numeric vector of length %d (%s).", name, n, what),
           call. = FALSE)
    }
  }
  check_len(x_intercepts, n_covariates_x, "x_intercepts", "one per X-grid covariate")
  check_len(y_intercepts, n_covariates_y, "y_intercepts", "one per Y-grid covariate")
  check_len(beta0_y, 1L, "beta0_y", "the outcome intercept")
  check_len(beta_x, n_covariates_x, "beta_x", "one per X-grid covariate")
  check_len(beta_y, n_covariates_y, "beta_y", "one per Y-grid covariate")

  set.seed(seed)

  res1 <- c(5L, 5L)    # Y grid resolution (coarser)
  res2 <- c(10L, 10L)  # X grid resolution (finer)

  # ---- build the two grids ------------------------------------------------
  grid1 <- sp::GridTopology(cellcentre.offset = c(0.5, 0.5),
                            cellsize = c(2, 2),
                            cells.dim = res1)
  sp_grid1 <- sp::SpatialGrid(grid1)

  grid2 <- sp::GridTopology(cellcentre.offset = c(0, 0),
                            cellsize = c(1, 1),
                            cells.dim = res2)
  sp_grid2 <- sp::SpatialGrid(grid2)

  sp_grid2_poly <- sf::st_as_sf(as(sp_grid2, "SpatialPolygons"))
  sp_grid1_poly <- sf::st_as_sf(as(sp_grid1, "SpatialPolygons"))

  # Create non-nested misalignment: in every other row of the fine grid, merge
  # a randomly chosen horizontal pair of cells that straddles a coarse-cell
  # boundary.
  pair_list <- list(c(2L, 3L), c(4L, 5L), c(6L, 7L), c(8L, 9L))
  union_ids <- matrix(0L, res2[1], res2[2])
  id_num <- 1L
  for (i in seq(1L, res2[1], by = 2L)) {
    union_ids[i, pair_list[[sample(1:4, size = 1)]]] <- id_num
    id_num <- id_num + 1L
  }
  sp_grid2_poly$union_ids <- c(t(union_ids))

  store_xunions <- NULL
  for (i in seq_len(max(sp_grid2_poly$union_ids))) {
    temp <- sf::st_sf(geometry = sp_grid2_poly[sp_grid2_poly$union_ids == i, ] |>
                        sf::st_union())
    store_xunions <- rbind(store_xunions, temp)
  }

  # Merge the remaining fine cells within each coarse cell
  grid2_nounion <- sp_grid2_poly[sp_grid2_poly$union_ids == 0L, ]
  store_iunions <- NULL
  for (i in seq_len(nrow(sp_grid1_poly))) {
    inside <- sf::st_contains(sp_grid1_poly[i, ], grid2_nounion, sparse = FALSE)[1, ]
    temp <- sf::st_sf(geometry = grid2_nounion[inside, ] |> sf::st_union())
    store_iunions <- rbind(store_iunions, temp)
  }

  gridx <- rbind(store_iunions, store_xunions)
  gridy <- sp_grid1_poly
  gridy$ID <- seq_len(nrow(gridy))
  gridx$ID <- seq_len(nrow(gridx))

  # ---- atoms: spatial intersection of the two grids -----------------------
  # create_atoms() calls sf::st_intersection(gridy, gridx) and drops the
  # zero-area point/line pieces produced where cell edges touch.
  atoms <- create_atoms(gridy, gridx, id_y = "ID", id_x = "ID")
  atoms$ID_atomorder <- seq_len(nrow(atoms))
  n_atoms <- nrow(atoms)

  # ---- atom-level populations --------------------------------------------
  draw_pop <- function() stats::rpois(n_atoms, lambda = 50) + 10L
  if (diff_pops_per_var) {
    for (i in seq_len(n_covariates_x)) {
      atoms[[paste0("population_covariate_x_", i)]] <- draw_pop()
    }
    for (i in seq_len(n_covariates_y)) {
      atoms[[paste0("population_covariate_y_", i)]] <- draw_pop()
    }
    atoms$population_y <- draw_pop()
  } else if (diff_pops) {
    atoms$population <- draw_pop()
  } else {
    atoms$population <- rep(1L, n_atoms)
  }
  pop_of <- function(col) if (diff_pops_per_var) atoms[[col]] else atoms$population

  # ---- atom-level spatial random effects ---------------------------------
  neighbors_atoms <- spdep::poly2nb(atoms, queen = TRUE)
  W_atoms <- spdep::nb2mat(neighbors_atoms, style = "B", zero.policy = TRUE)

  if (!xy_cov_cor) {
    x_atom_results <- gen_correlated_spat(W_atoms, n_covariates_x,
                                          correlation = x_correlation, rho = rho_x)
    y_atom_results <- gen_correlated_spat(W_atoms, n_covariates_y,
                                          correlation = y_correlation, rho = rho_y)
  } else {
    xy_atom_results <- gen_correlated_spat(W_atoms, n_covariates_x + n_covariates_y,
                                           correlation = x_correlation, rho = rho_x)
    x_atom_results <- xy_atom_results[, seq_len(n_covariates_x), drop = FALSE]
    y_atom_results <- xy_atom_results[, n_covariates_x + seq_len(n_covariates_y),
                                      drop = FALSE]
  }

  # Draw one atom-level covariate given its distribution, intercept, random
  # effect and population denominator.
  draw_covariate <- function(dist, intercept, re, pop) {
    switch(dist,
      pois  = stats::rpois(n_atoms, lambda = pop * exp(intercept + re)),
      binom = stats::rbinom(n_atoms, size = pop,
                            prob = pmin(pmax(stats::plogis(intercept + re), 0.001), 0.999)),
      norm  = pop * stats::rnorm(n_atoms, mean = intercept + re, sd = 1)
    )
  }

  for (i in seq_len(n_covariates_x)) {
    atoms[[paste0("covariate_x_", i)]] <- draw_covariate(
      dist_covariates_x[i], x_intercepts[i], x_atom_results[, i],
      pop_of(paste0("population_covariate_x_", i)))
  }
  for (i in seq_len(n_covariates_y)) {
    atoms[[paste0("covariate_y_", i)]] <- draw_covariate(
      dist_covariates_y[i], y_intercepts[i], y_atom_results[, i],
      pop_of(paste0("population_covariate_y_", i)))
  }

  # ---- atom-level outcome --------------------------------------------------
  outcome_re  <- gen_correlated_spat(W_atoms, 1L, rho = rho_y)[, 1]
  pop_outcome <- pop_of("population_y")

  linear_pred_y_atom <- rep(beta0_y, n_atoms)
  for (i in seq_len(n_covariates_x)) {
    rate <- atoms[[paste0("covariate_x_", i)]] / pop_of(paste0("population_covariate_x_", i))
    linear_pred_y_atom <- linear_pred_y_atom + beta_x[i] * rate
  }
  for (i in seq_len(n_covariates_y)) {
    rate <- atoms[[paste0("covariate_y_", i)]] / pop_of(paste0("population_covariate_y_", i))
    linear_pred_y_atom <- linear_pred_y_atom + beta_y[i] * rate
  }
  eta <- linear_pred_y_atom + outcome_re

  atoms$y <- switch(dist_y,
    pois  = stats::rpois(n_atoms, lambda = pop_outcome * exp(eta)),
    binom = stats::rbinom(n_atoms, size = pop_outcome, prob = stats::plogis(eta)),
    norm  = stats::rnorm(n_atoms, mean = pop_outcome * eta, sd = pop_outcome)
  )

  # ---- aggregate to grid level ----------------------------------------------
  S_x <- nrow(gridx)
  S_y <- nrow(gridy)
  for (i in seq_len(n_covariates_x)) {
    v <- paste0("covariate_x_", i)
    gridx[[v]] <- .sum_by_id(atoms[[v]], atoms$ID_x, S_x)
  }
  for (i in seq_len(n_covariates_y)) {
    v <- paste0("covariate_y_", i)
    gridy[[v]] <- .sum_by_id(atoms[[v]], atoms$ID_y, S_y)
  }
  gridy$y <- .sum_by_id(atoms$y, atoms$ID_y, S_y)

  if (diff_pops_per_var) {
    for (i in seq_len(n_covariates_x)) {
      v <- paste0("population_covariate_x_", i)
      gridx[[v]] <- .sum_by_id(atoms[[v]], atoms$ID_x, S_x)
    }
    for (i in seq_len(n_covariates_y)) {
      v <- paste0("population_covariate_y_", i)
      gridy[[v]] <- .sum_by_id(atoms[[v]], atoms$ID_y, S_y)
    }
    gridy$population_y <- .sum_by_id(atoms$population_y, atoms$ID_y, S_y)
  } else {
    gridy$population <- .sum_by_id(atoms$population, atoms$ID_y, S_y)
    gridx$population <- .sum_by_id(atoms$population, atoms$ID_x, S_x)
  }

  true_params <- abrm_params(
    beta0_y = beta0_y,
    beta_x = beta_x,
    beta_y = beta_y,
    x_correlation = x_correlation,
    y_correlation = y_correlation
  )

  structure(
    list(gridy = gridy, gridx = gridx, atoms = atoms, true_params = true_params),
    class = "misaligned_data"
  )
}
