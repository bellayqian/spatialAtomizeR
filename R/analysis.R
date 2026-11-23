#' Run ABRM Analysis
#'
#' Fits an Atom-Based Regression Model (ABRM) to spatially misaligned data
#' using Bayesian MCMC via NIMBLE. Works with both simulated and real-world
#' data. Variables can have any column names: tell \code{run_abrm()} which
#' columns to use through \code{x_vars}, \code{y_vars}, \code{outcome},
#' \code{atom_vars}, the ID arguments and \code{population}.
#'
#' @param gridx The X-grid \code{sf} data frame, containing an area ID column
#'   (see \code{id_x}) and the X-grid covariates (see \code{x_vars}).
#' @param gridy The Y-grid \code{sf} data frame, containing an area ID column
#'   (see \code{id_y}), the Y-grid covariates (see \code{y_vars}) and the
#'   outcome (see \code{outcome}).
#' @param atoms The atom \code{sf} data frame (for example from
#'   \code{\link{create_atoms}}). Must contain the X-grid and Y-grid ID of
#'   each atom (see \code{atom_id_x}, \code{atom_id_y}), the population
#'   denominator column(s) (see \code{population}) and any atom-level
#'   covariates (see \code{atom_vars}).
#' @param model_code NIMBLE model code from \code{\link{get_abrm_model}}.
#' @param x_types Character vector with one element per X-grid covariate (in
#'   the order of \code{x_vars}) giving its distribution: \code{"norm"},
#'   \code{"pois"} or \code{"binom"}. The long forms \code{"normal"},
#'   \code{"poisson"} and \code{"binomial"} are also accepted.
#' @param y_types As \code{x_types}, for the Y-grid covariates (in the order of
#'   \code{y_vars}).
#' @param dist_y Distribution of the outcome: \code{"norm"}, \code{"pois"}
#'   (default) or \code{"binom"} (long forms also accepted).
#' @param x_vars Character vector of X-grid covariate column names in
#'   \code{gridx}. If \code{NULL} (default), all columns named
#'   \code{covariate_x_1}, \code{covariate_x_2}, ... are used.
#' @param y_vars Character vector of Y-grid covariate column names in
#'   \code{gridy}. If \code{NULL} (default), all columns named
#'   \code{covariate_y_1}, \code{covariate_y_2}, ... are used.
#' @param outcome Name of the outcome column in \code{gridy}. Default \code{"y"}.
#' @param atom_vars Optional character vector of atom-level covariate column
#'   names in \code{atoms}. If \code{NULL} (default), any columns named
#'   \code{covariate_a_1}, \code{covariate_a_2}, ... are used; use
#'   \code{character(0)} to fit without atom-level covariates.
#' @param id_x,id_y Names of the area ID columns in \code{gridx} and
#'   \code{gridy}. IDs may be integer, numeric or character but must be unique
#'   within each grid. Default \code{"ID"}.
#' @param atom_id_x,atom_id_y Names of the columns in \code{atoms} holding each
#'   atom's X-grid and Y-grid ID. Defaults \code{"ID_x"} and \code{"ID_y"}
#'   (the names produced by \code{\link{create_atoms}}).
#' @param population Population denominator(s) in \code{atoms}. Either
#'   \itemize{
#'     \item a single column name used for every variable, or
#'     \item a named list \code{list(x = , y = , outcome = )} giving one column
#'       name per X-grid covariate, one per Y-grid covariate (each may be a
#'       single name to recycle) and one for the outcome, or
#'     \item \code{NULL} (default): use a column named \code{population} if
#'       present; otherwise per-variable columns named
#'       \code{population_covariate_x_1}, ..., \code{population_covariate_y_1},
#'       ..., and \code{population_y} (as produced by
#'       \code{\link{simulate_misaligned_data}}).
#'   }
#' @param true_params Optional true values of the outcome-model coefficients,
#'   used to report bias and credible-interval coverage. Build it with
#'   \code{\link{abrm_params}} (recommended), or pass the
#'   \code{true_params} component returned by
#'   \code{\link{simulate_misaligned_data}}, which has the same format. A plain
#'   named list with the same elements is also accepted. See
#'   \strong{Format of \code{true_params}} below.
#' @param niter Number of MCMC iterations (default: 50000)
#' @param nburnin Number of burn-in iterations (default: 30000)
#' @param nchains Number of MCMC chains (default: 2)
#' @param thin Thinning interval (default: 10)
#' @param seed Integer seed for reproducibility. Each chain uses seed+(chain_number-1) (default: NULL)
#' @param sim_metadata Optional named list of simulation metadata (e.g.,
#'   \code{sim_number}, \code{x_correlation}, \code{y_correlation}) passed
#'   through to \code{\link{run_nimble_model}} and used to label diagnostic
#'   plot file names when \code{save_plots = TRUE}. For non-simulation use,
#'   leave as \code{NULL}.
#' @param save_plots Logical, whether to save diagnostic plots to a PDF (default: FALSE)
#' @param output_dir Directory for saving outputs (default: NULL)
#' @param compute_waic Logical; if \code{TRUE}, NIMBLE computes the Widely
#'   Applicable Information Criterion (WAIC) during MCMC sampling and stores it
#'   in the returned object. Retrieve it afterwards with \code{\link{waic}}.
#'   Default \code{FALSE}.
#' @param ... Only used to support the deprecated arguments
#'   \code{norm_idx_x}, \code{pois_idx_x}, \code{binom_idx_x},
#'   \code{norm_idx_y}, \code{pois_idx_y} and \code{binom_idx_y}, which are
#'   converted to \code{x_types} / \code{y_types} with a warning. Any other
#'   argument is an error.
#'
#' @section Format of \code{true_params}:
#' \code{true_params} is a list, best created with \code{\link{abrm_params}()},
#' with the following elements. Any of them may be omitted; coefficients
#' without a true value get \code{NA} in the results and are left out of the
#' bias and coverage summaries.
#' \describe{
#'   \item{\code{beta0_y}}{Numeric scalar: the outcome-model intercept.
#'     (\code{beta_0_y} is accepted as an alias.)}
#'   \item{\code{beta_x}}{Numeric vector of length \code{length(x_vars)}: the
#'     coefficients of the X-grid covariates. Unnamed vectors are matched by
#'     position (the order of \code{x_vars}); named vectors are matched by
#'     covariate name.}
#'   \item{\code{beta_y}}{Numeric vector of length \code{length(y_vars)}: the
#'     coefficients of the Y-grid covariates, matched as for \code{beta_x}.}
#'   \item{\code{beta_a}}{Numeric vector of length \code{length(atom_vars)}:
#'     the coefficients of the atom-level covariates, matched as for
#'     \code{beta_x}.}
#' }
#' Other elements (such as \code{x_correlation} from
#' \code{simulate_misaligned_data()}) are ignored. For example:
#' \preformatted{
#' true_params = abrm_params(beta0_y = -1,
#'                           beta_x  = c(0.10, -0.05),
#'                           beta_y  = c(income = -0.10, pct_old = 0.05))
#' }
#'
#' @return An object of class \code{"abrm"}: a named list with components
#'   \describe{
#'     \item{\code{mcmc_results}}{Output of \code{\link{run_nimble_model}}:
#'       MCMC samples, posterior summary, convergence diagnostics (see
#'       \code{\link{mcmc_diagnostics}}) and, optionally, WAIC.}
#'     \item{\code{parameter_estimates}}{Data frame with one row per
#'       regression coefficient: \code{variable} (the covariate name, or
#'       \code{"intercept"}), \code{parameter} (the internal parameter name,
#'       e.g. \code{"beta_x[1]"}), \code{estimated_beta}, \code{std_error},
#'       \code{ci_lower}, \code{ci_upper}, and, when \code{true_params} is
#'       supplied, \code{true_beta}, \code{bias}, \code{relative_bias} and
#'       \code{within_ci}.}
#'     \item{\code{all_parameters}}{Posterior summaries of all monitored
#'       parameters.}
#'     \item{\code{fitted_values}}{Numeric vector of Y-grid-level predicted
#'       outcome values on the original outcome scale.}
#'     \item{\code{y_grid_ids}}{Y-grid cell IDs (as in \code{gridy[[id_y]]})
#'       in model order.}
#'     \item{\code{y_observed}}{Observed outcome values in model order.}
#'     \item{\code{param_labels}}{Named character vector mapping internal
#'       parameter names to covariate names.}
#'     \item{\code{model_info}}{List recording \code{x_vars}, \code{y_vars},
#'       \code{atom_vars}, \code{outcome}, \code{x_types}, \code{y_types} and
#'       \code{dist_y}.}
#'   }
#'   Use \code{fitted()} for a table of observed vs. fitted values,
#'   \code{\link{plot_trace}} / \code{\link{plot_density}} for MCMC plots,
#'   \code{\link{mcmc_diagnostics}} for R-hat and effective sample sizes, and
#'   \code{waic()} for model fit criteria when \code{compute_waic = TRUE}.
#'
#' @examples
#' \donttest{
#'   # Simulate misaligned spatial data with one normal covariate per grid
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
#'
#'   # Retrieve the NIMBLE model code
#'   model_code <- get_abrm_model()
#'
#'   # Fit the ABRM (use small niter/nburnin for illustration only)
#'   results <- run_abrm(
#'     gridx       = sim_data$gridx,
#'     gridy       = sim_data$gridy,
#'     atoms       = sim_data$atoms,
#'     model_code  = model_code,
#'     x_types     = "norm",
#'     y_types     = "norm",
#'     dist_y      = "norm",
#'     true_params = abrm_params(beta0_y = 0, beta_x = 0.1, beta_y = -0.1),
#'     niter       = 1000,
#'     nburnin     = 500,
#'     nchains     = 2,
#'     seed        = 1
#'   )
#'
#'   print(results)             # concise model summary
#'   summary(results)           # full parameter table
#'   plot_trace(results)        # MCMC trace plots (a ggplot object)
#'   plot_density(results)      # posterior density plots (a ggplot object)
#'   mcmc_diagnostics(results)  # R-hat and effective sample sizes
#'
#'   # Covariates can have any names
#'   gx <- sim_data$gridx; names(gx)[names(gx) == "covariate_x_1"] <- "income"
#'   gy <- sim_data$gridy; names(gy)[names(gy) == "covariate_y_1"] <- "age"
#'   results2 <- run_abrm(
#'     gridx = gx, gridy = gy, atoms = sim_data$atoms, model_code = model_code,
#'     x_vars = "income", y_vars = "age", outcome = "y",
#'     x_types = "norm", y_types = "norm", dist_y = "norm",
#'     niter = 1000, nburnin = 500, seed = 1
#'   )
#'   coef(results2)
#' }
#'
#' @export
run_abrm <- function(gridx, gridy, atoms, model_code,
                     x_types, y_types, dist_y = "pois",
                     x_vars = NULL, y_vars = NULL, outcome = "y",
                     atom_vars = NULL,
                     id_x = "ID", id_y = "ID",
                     atom_id_x = "ID_x", atom_id_y = "ID_y",
                     population = NULL,
                     true_params = NULL,
                     niter = 50000, nburnin = 30000,
                     nchains = 2, thin = 10, seed = NULL,
                     sim_metadata = NULL, save_plots = FALSE, output_dir = NULL,
                     compute_waic = FALSE, ...) {

  # ---- deprecated / unknown arguments --------------------------------------
  dots <- list(...)
  old_args <- c("norm_idx_x", "pois_idx_x", "binom_idx_x",
                "norm_idx_y", "pois_idx_y", "binom_idx_y")
  unknown <- setdiff(names(dots), old_args)
  if (length(dots) > 0L && (is.null(names(dots)) || any(names(dots) == "") ||
                            length(unknown) > 0L)) {
    stop("Unused argument(s) in run_abrm(): ",
         paste(if (length(unknown)) unknown else "<unnamed>", collapse = ", "),
         call. = FALSE)
  }

  # ---- input checks: grids and atoms ---------------------------------------
  for (nm in c("gridx", "gridy", "atoms")) {
    if (!inherits(get(nm), "sf")) stop(sprintf("'%s' must be an sf data frame.", nm),
                                       call. = FALSE)
  }
  need_col <- function(df, col, df_name, what) {
    if (!is.character(col) || length(col) != 1L || !col %in% names(df)) {
      stop(sprintf("%s: column '%s' not found in %s.", what, col, df_name), call. = FALSE)
    }
  }
  need_col(gridx, id_x, "gridx", "id_x")
  need_col(gridy, id_y, "gridy", "id_y")
  need_col(atoms, atom_id_x, "atoms", "atom_id_x")
  need_col(atoms, atom_id_y, "atoms", "atom_id_y")
  need_col(gridy, outcome, "gridy", "outcome")

  # ---- variable names --------------------------------------------------------
  legacy_vars <- function(df, prefix) {
    v <- grep(paste0("^", prefix, "\\d+$"), names(df), value = TRUE)
    v[order(as.integer(sub(prefix, "", v, fixed = TRUE)))]
  }
  if (is.null(x_vars)) x_vars <- legacy_vars(gridx, "covariate_x_")
  if (is.null(y_vars)) y_vars <- legacy_vars(gridy, "covariate_y_")
  if (is.null(atom_vars)) atom_vars <- legacy_vars(atoms, "covariate_a_")
  if (length(x_vars) == 0L) {
    stop("No X-grid covariates: supply their column names in gridx via 'x_vars'.",
         call. = FALSE)
  }
  if (length(y_vars) == 0L) {
    stop("No Y-grid covariates: supply their column names in gridy via 'y_vars'.",
         call. = FALSE)
  }
  check_vars <- function(df, vars, df_name, arg) {
    missing_v <- setdiff(vars, names(df))
    if (length(missing_v)) {
      stop(sprintf("Column(s) named in '%s' not found in %s: %s", arg, df_name,
                   paste(missing_v, collapse = ", ")), call. = FALSE)
    }
    if (anyDuplicated(vars)) stop(sprintf("'%s' contains duplicates.", arg), call. = FALSE)
    for (v in vars) {
      if (!is.numeric(df[[v]])) {
        stop(sprintf("Column '%s' in %s must be numeric.", v, df_name), call. = FALSE)
      }
      if (anyNA(df[[v]])) {
        stop(sprintf("Column '%s' in %s contains missing values, which are not supported.",
                     v, df_name), call. = FALSE)
      }
    }
  }
  check_vars(gridx, x_vars, "gridx", "x_vars")
  check_vars(gridy, c(y_vars, outcome), "gridy", "y_vars/outcome")
  check_vars(atoms, atom_vars, "atoms", "atom_vars")
  p_x <- length(x_vars)
  p_y <- length(y_vars)
  p_a <- length(atom_vars)

  # ---- variable types --------------------------------------------------------
  types_from_idx <- function(side, p) {
    idx <- dots[paste0(c("norm", "pois", "binom"), "_idx_", side)]
    if (all(vapply(idx, is.null, logical(1)))) return(NULL)
    out <- rep(NA_character_, p)
    out[idx[[1]]] <- "norm"; out[idx[[2]]] <- "pois"; out[idx[[3]]] <- "binom"
    if (length(out) != p) {
      stop(sprintf("Deprecated *_idx_%s arguments refer to covariates beyond the %d supplied.",
                   side, p), call. = FALSE)
    }
    out
  }
  resolve_types <- function(types, side, p, vars) {
    arg <- paste0(side, "_types")
    from_idx <- types_from_idx(side, p)
    if (is.null(types)) {
      if (is.null(from_idx)) {
        stop(sprintf(paste0("'%s' must be supplied: a character vector of length %d ",
                            "with elements in \"norm\", \"pois\", \"binom\" (one per ",
                            "covariate in %s_vars: %s)."),
                     arg, p, side, paste(vars, collapse = ", ")), call. = FALSE)
      }
      warning(sprintf(paste0("The norm_idx_%1$s/pois_idx_%1$s/binom_idx_%1$s arguments ",
                             "are deprecated; use %2$s = c(%3$s) instead."),
                      side, arg, paste0("\"", from_idx, "\"", collapse = ", ")),
              call. = FALSE)
      types <- from_idx
    } else if (!is.null(from_idx)) {
      stop(sprintf("Supply either '%s' or the deprecated *_idx_%s arguments, not both.",
                   arg, side), call. = FALSE)
    }
    if (anyNA(types)) {
      stop(sprintf("No distribution given for %s covariate(s): %s", side,
                   paste(vars[is.na(types)], collapse = ", ")), call. = FALSE)
    }
    types <- .normalize_dist(types, arg)
    if (length(types) != p) {
      stop(sprintf("'%s' has length %d but there are %d %s-grid covariates (%s).",
                   arg, length(types), p, toupper(side), paste(vars, collapse = ", ")),
           call. = FALSE)
    }
    types
  }
  x_types <- resolve_types(if (missing(x_types)) NULL else x_types, "x", p_x, x_vars)
  y_types <- resolve_types(if (missing(y_types)) NULL else y_types, "y", p_y, y_vars)

  if (is.numeric(dist_y)) {
    if (length(dist_y) != 1L || !dist_y %in% 1:3) {
      stop("'dist_y' must be \"norm\", \"pois\" or \"binom\".", call. = FALSE)
    }
    new_dist <- names(.dist_code)[dist_y]
    warning(sprintf("Numeric 'dist_y' is deprecated; use dist_y = \"%s\" instead.", new_dist),
            call. = FALSE)
    dist_y <- new_dist
  }
  dist_y <- .normalize_dist(dist_y, "dist_y")
  if (length(dist_y) != 1L) stop("'dist_y' must be a single value.", call. = FALSE)

  # ---- IDs: map user IDs to internal integer indices ------------------------
  ids_x <- gridx[[id_x]]
  ids_y <- gridy[[id_y]]
  if (anyDuplicated(ids_x) || anyNA(ids_x)) {
    stop(sprintf("gridx$%s must be unique and non-missing.", id_x), call. = FALSE)
  }
  if (anyDuplicated(ids_y) || anyNA(ids_y)) {
    stop(sprintf("gridy$%s must be unique and non-missing.", id_y), call. = FALSE)
  }
  atom_ix <- match(atoms[[atom_id_x]], ids_x)
  atom_iy <- match(atoms[[atom_id_y]], ids_y)
  if (anyNA(atom_ix)) {
    stop(sprintf("%d atom(s) have an %s value that does not match any gridx$%s.",
                 sum(is.na(atom_ix)), atom_id_x, id_x), call. = FALSE)
  }
  if (anyNA(atom_iy)) {
    stop(sprintf("%d atom(s) have an %s value that does not match any gridy$%s.",
                 sum(is.na(atom_iy)), atom_id_y, id_y), call. = FALSE)
  }
  if (length(setdiff(seq_along(ids_x), atom_ix))) {
    stop("Some gridx cells contain no atoms; every grid cell must contain at least one atom.",
         call. = FALSE)
  }
  if (length(setdiff(seq_along(ids_y), atom_iy))) {
    stop("Some gridy cells contain no atoms; every grid cell must contain at least one atom.",
         call. = FALSE)
  }

  # ---- population columns ----------------------------------------------------
  pop_cols <- .resolve_population(population, atoms, p_x, p_y)

  # ---- standardised inputs for the internal pipeline -------------------------
  build_sf <- function(cols, geom) {
    df <- as.data.frame(cols, optional = TRUE)
    sf::st_sf(df, geometry = geom)
  }
  std_gridx <- build_sf(
    c(list(ID = seq_along(ids_x)),
      stats::setNames(lapply(x_vars, function(v) as.numeric(gridx[[v]])),
                      paste0("covariate_x_", seq_len(p_x)))),
    sf::st_geometry(gridx))
  std_gridy <- build_sf(
    c(list(ID = seq_along(ids_y)),
      stats::setNames(lapply(y_vars, function(v) as.numeric(gridy[[v]])),
                      paste0("covariate_y_", seq_len(p_y))),
      list(y = as.numeric(gridy[[outcome]]))),
    sf::st_geometry(gridy))
  atom_cols <- list(ID_x = atom_ix, ID_y = atom_iy)
  if (pop_cols$single) {
    atom_cols$population <- atoms[[pop_cols$x[1]]]
  } else {
    atom_cols <- c(atom_cols,
      stats::setNames(lapply(pop_cols$x, function(v) atoms[[v]]),
                      paste0("population_covariate_x_", seq_len(p_x))),
      stats::setNames(lapply(pop_cols$y, function(v) atoms[[v]]),
                      paste0("population_covariate_y_", seq_len(p_y))),
      list(population_y = atoms[[pop_cols$outcome]]))
  }
  if (p_a > 0L) {
    atom_cols <- c(atom_cols,
      stats::setNames(lapply(atom_vars, function(v) as.numeric(atoms[[v]])),
                      paste0("covariate_a_", seq_len(p_a))))
  }
  std_atoms <- build_sf(atom_cols, sf::st_geometry(atoms))
  std_data <- list(gridx = std_gridx, gridy = std_gridy, atoms = std_atoms)

  message("Preparing ABRM inputs...")
  bookkeeping   <- prepare_spatial_bookkeeping(std_data)
  adjacency     <- prepare_adjacency_matrices(bookkeeping$gridy_yorder,
                                              bookkeeping$gridx_xorder)
  nimble_inputs <- prepare_nimble_inputs(bookkeeping, adjacency, std_data,
                                         x_types = x_types, y_types = y_types,
                                         dist_y = dist_y)

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
    output_dir = output_dir,
    compute_waic = compute_waic
  )
  # run_nimble_model() has already split NIMBLE's stacked beta_y vector into
  # beta_x[], beta_y[] and beta_a[] in the samples, summary and diagnostics.

  # ---- labels -----------------------------------------------------------------
  raw_names <- c("beta_0_y",
                 sprintf("beta_x[%d]", seq_len(p_x)),
                 sprintf("beta_y[%d]", seq_len(p_y)),
                 sprintf("beta_a[%d]", seq_len(p_a)))
  labels <- c("intercept", x_vars, y_vars, atom_vars)
  if (anyDuplicated(labels)) {
    # The same column name can occur on both grids; disambiguate.
    labels <- c("intercept", paste0("x:", x_vars), paste0("y:", y_vars),
                paste0("a:", atom_vars))
  }
  param_labels <- stats::setNames(labels, raw_names)

  # ---- parameter tables ---------------------------------------------------------
  summ <- abrm_results$summary$all.chains
  all_parameters <- data.frame(
    variable = rownames(summ),
    estimated_beta = summ[, "Mean"],
    std_error = summ[, "St.Dev."],
    ci_lower = summ[, "95%CI_low"],
    ci_upper = summ[, "95%CI_upp"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Posterior mean of the atom-level linear predictor, used for fitted values
  linpred <- all_parameters$estimated_beta[grepl("^linear_pred_y", all_parameters$variable)]
  all_parameters <- all_parameters[!grepl("^linear_pred_y", all_parameters$variable), ,
                                   drop = FALSE]
  rownames(all_parameters) <- NULL

  beta_rows <- match(raw_names, all_parameters$variable)
  parameter_estimates <- data.frame(
    variable = labels,
    parameter = raw_names,
    all_parameters[beta_rows, c("estimated_beta", "std_error", "ci_lower", "ci_upper")],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (!is.null(true_params)) {
    truth <- .resolve_true_params(true_params, x_vars, y_vars, atom_vars)
    pe <- parameter_estimates
    pe$true_beta <- truth
    pe$bias <- pe$estimated_beta - pe$true_beta
    pe$relative_bias <- (pe$estimated_beta - pe$true_beta) / pe$true_beta * 100
    pe$within_ci <- pe$true_beta >= pe$ci_lower & pe$true_beta <= pe$ci_upper
    parameter_estimates <- pe
  }

  # ---- fitted values -------------------------------------------------------------
  pop_out <- nimble_inputs$constants$pop_atoms_outcome
  fitted_atom <- switch(dist_y,
    norm  = pop_out * linpred,
    pois  = pop_out * exp(linpred),
    binom = pop_out * stats::plogis(linpred)
  )
  J_y      <- nimble_inputs$constants$J_y
  ylat_ind <- nimble_inputs$constants$ylatent_ind
  fitted_values <- c(
    fitted_atom[seq_len(J_y)],
    apply(ylat_ind, 1, function(i) sum(fitted_atom[(i[1] + J_y):(i[2] + J_y)]))
  )

  structure(
    list(
      mcmc_results = abrm_results,
      parameter_estimates = parameter_estimates,
      all_parameters = all_parameters,
      fitted_values = fitted_values,
      y_grid_ids = ids_y[bookkeeping$gridy_yorder$ID_y],
      y_observed = nimble_inputs$data$y_obs,
      param_labels = param_labels,
      model_info = list(x_vars = x_vars, y_vars = y_vars, atom_vars = atom_vars,
                        outcome = outcome, x_types = x_types, y_types = y_types,
                        dist_y = dist_y)
    ),
    class = "abrm"
  )
}


#' Work out which atom columns hold the population denominators
#' @return list(single = TRUE/FALSE, x = <p_x names>, y = <p_y names>,
#'   outcome = <name>); when single, all three point to the same column.
#' @keywords internal
#' @noRd
.resolve_population <- function(population, atoms, p_x, p_y) {
  nms <- names(atoms)
  check_pop <- function(cols) {
    missing_c <- setdiff(cols, nms)
    if (length(missing_c)) {
      stop("Population column(s) not found in atoms: ",
           paste(unique(missing_c), collapse = ", "), call. = FALSE)
    }
    for (v in unique(cols)) {
      pv <- atoms[[v]]
      if (!is.numeric(pv) || anyNA(pv) || any(pv <= 0)) {
        stop(sprintf("Population column '%s' must be numeric, positive and non-missing.", v),
             call. = FALSE)
      }
    }
  }

  if (is.null(population)) {
    legacy_x <- paste0("population_covariate_x_", seq_len(p_x))
    legacy_y <- paste0("population_covariate_y_", seq_len(p_y))
    has_single <- "population" %in% nms
    has_per_var <- "population_covariate_x_1" %in% nms
    if (has_single && has_per_var) {
      stop(paste0("atoms has both 'population' and 'population_covariate_x_1'. ",
                  "Say which to use with the 'population' argument."), call. = FALSE)
    }
    if (has_single) {
      population <- "population"
    } else if (has_per_var) {
      population <- list(x = legacy_x, y = legacy_y, outcome = "population_y")
    } else {
      stop(paste0("No population column found in atoms. Name it with the 'population' ",
                  "argument, e.g. population = \"pop\", or ",
                  "population = list(x = ..., y = ..., outcome = ...)."), call. = FALSE)
    }
  }

  if (is.character(population)) {
    if (length(population) != 1L) {
      stop("A character 'population' must be a single column name.", call. = FALSE)
    }
    check_pop(population)
    return(list(single = TRUE, x = rep(population, p_x), y = rep(population, p_y),
                outcome = population))
  }

  if (!is.list(population) || !all(c("x", "y", "outcome") %in% names(population))) {
    stop("'population' must be a column name or a list with elements x, y and outcome.",
         call. = FALSE)
  }
  recycle <- function(v, p, what) {
    if (length(v) == 1L) v <- rep(v, p)
    if (length(v) != p) {
      stop(sprintf("population$%s must have length 1 or %d.", what, p), call. = FALSE)
    }
    as.character(v)
  }
  px <- recycle(population$x, p_x, "x")
  py <- recycle(population$y, p_y, "y")
  po <- as.character(population$outcome)
  if (length(po) != 1L) stop("population$outcome must be a single column name.", call. = FALSE)
  check_pop(c(px, py, po))
  list(single = FALSE, x = px, y = py, outcome = po)
}


#' Align user-supplied true parameters with the model's coefficients
#' @return Numeric vector: intercept, X, Y, atom-level (NA where unknown)
#' @keywords internal
#' @noRd
.resolve_true_params <- function(true_params, x_vars, y_vars, atom_vars) {
  if (!is.list(true_params)) {
    stop("'true_params' must be an abrm_params() object (or a list with the same ",
         "elements); see ?abrm_params.", call. = FALSE)
  }
  known <- c("beta0_y", "beta_0_y", "beta_x", "beta_y", "beta_a")
  if (!inherits(true_params, "abrm_params") &&
      !any(names(true_params) %in% known)) {
    stop("'true_params' has none of the expected elements (beta0_y, beta_x, ",
         "beta_y, beta_a); see ?abrm_params.", call. = FALSE)
  }
  align <- function(v, vars, el) {
    n <- length(vars)
    if (is.null(v)) return(rep(NA_real_, n))
    if (!is.numeric(v) || length(v) != n) {
      stop(sprintf("true_params$%s must be a numeric vector of length %d (one per covariate: %s).",
                   el, n, paste(vars, collapse = ", ")), call. = FALSE)
    }
    if (!is.null(names(v)) && all(names(v) != "")) {
      if (!setequal(names(v), vars)) {
        stop(sprintf("Names of true_params$%s (%s) do not match the covariates (%s).",
                     el, paste(names(v), collapse = ", "), paste(vars, collapse = ", ")),
             call. = FALSE)
      }
      v <- v[vars]
    }
    unname(as.numeric(v))
  }
  b0 <- if (!is.null(true_params$beta0_y)) true_params$beta0_y else true_params$beta_0_y
  if (is.null(b0)) b0 <- NA_real_
  if (length(b0) != 1L || !(is.numeric(b0) || is.na(b0))) {
    stop("true_params$beta0_y must be a single number.", call. = FALSE)
  }
  c(as.numeric(b0),
    align(true_params$beta_x, x_vars, "beta_x"),
    align(true_params$beta_y, y_vars, "beta_y"),
    align(true_params$beta_a, atom_vars, "beta_a"))
}
