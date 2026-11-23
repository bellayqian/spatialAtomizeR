#' Construct True Parameter Values for \code{run_abrm()}
#'
#' Builds the \code{true_params} argument of \code{\link{run_abrm}}: the true
#' regression coefficients of the outcome model. When \code{true_params} is
#' supplied, \code{run_abrm()} reports the bias, relative bias and
#' credible-interval coverage of each estimated coefficient. This is mainly
#' useful with simulated data, where the truth is known;
#' \code{\link{simulate_misaligned_data}} already returns its truth in this
#' format.
#'
#' @details
#' The outcome model is
#' \deqn{g(\mu) = \beta_0 + \sum_j \beta_{x,j} x_j + \sum_k \beta_{y,k} y_k +
#'   \sum_l \beta_{a,l} a_l + \ldots}
#' for X-grid covariates \eqn{x_j}, Y-grid covariates \eqn{y_k} and (optional)
#' atom-level covariates \eqn{a_l}. Each of \code{beta_x}, \code{beta_y} and
#' \code{beta_a} is either
#' \itemize{
#'   \item an \strong{unnamed} vector, matched to the covariates by
#'     \emph{position} (the order of \code{x_vars}, \code{y_vars} and
#'     \code{atom_vars} in \code{run_abrm()}), or
#'   \item a \strong{named} vector, matched by covariate \emph{name}, in any
#'     order. The names must be exactly the covariate names.
#' }
#' Components left as \code{NULL} are treated as unknown: the corresponding
#' rows of the results have \code{NA} true values and are excluded from the
#' bias and coverage summaries.
#'
#' Lengths and names are checked against the data when the object is passed to
#' \code{run_abrm()}.
#'
#' @param beta0_y Numeric scalar, the intercept of the outcome model.
#' @param beta_x Numeric vector of coefficients for the X-grid covariates.
#' @param beta_y Numeric vector of coefficients for the Y-grid covariates.
#' @param beta_a Numeric vector of coefficients for atom-level covariates
#'   (only needed if the model uses \code{atom_vars}).
#' @param ... Optional further named values to carry along, such as
#'   \code{x_correlation} (used by \code{simulate_misaligned_data()}). They
#'   are stored but ignored by \code{run_abrm()}.
#'
#' @return An object of class \code{"abrm_params"}: a named list with the
#'   components \code{beta0_y}, \code{beta_x}, \code{beta_y} and \code{beta_a}
#'   (\code{NULL} when not supplied) plus anything passed through \code{...}.
#'   Elements are accessed as for any list, e.g. \code{p$beta_x}.
#'
#' @examples
#' # Matched by position
#' p <- abrm_params(beta0_y = -1, beta_x = c(0.1, -0.05), beta_y = 0.2)
#' p
#'
#' # Matched by covariate name (any order)
#' abrm_params(beta0_y = -1,
#'             beta_x = c(clinics = 0.1, income = -0.05),
#'             beta_y = c(age = 0.2))
#'
#' \donttest{
#'   sim <- simulate_misaligned_data(
#'     seed = 1, dist_covariates_x = "norm", dist_covariates_y = "norm",
#'     dist_y = "norm", beta0_y = 0, beta_x = 0.1, beta_y = -0.1
#'   )
#'   fit <- run_abrm(
#'     sim$gridx, sim$gridy, sim$atoms, get_abrm_model(),
#'     x_types = "norm", y_types = "norm", dist_y = "norm",
#'     true_params = abrm_params(beta0_y = 0, beta_x = 0.1, beta_y = -0.1),
#'     niter = 1000, nburnin = 500, seed = 1
#'   )
#'   fit$parameter_estimates
#' }
#'
#' @export
abrm_params <- function(beta0_y = NULL, beta_x = NULL, beta_y = NULL,
                        beta_a = NULL, ...) {
  check_scalar <- function(v, nm) {
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || is.na(v))) {
      stop(sprintf("'%s' must be a single non-missing number.", nm), call. = FALSE)
    }
  }
  check_vec <- function(v, nm) {
    if (is.null(v)) return(invisible())
    if (!is.numeric(v) || length(v) == 0L || anyNA(v)) {
      stop(sprintf("'%s' must be a non-empty numeric vector without missing values.", nm),
           call. = FALSE)
    }
    nms <- names(v)
    if (!is.null(nms) && (any(nms == "") || anyDuplicated(nms))) {
      stop(sprintf(paste0("If '%s' is named, every element needs a unique name ",
                          "(covariate names); otherwise leave it unnamed."), nm),
           call. = FALSE)
    }
  }
  check_scalar(beta0_y, "beta0_y")
  check_vec(beta_x, "beta_x")
  check_vec(beta_y, "beta_y")
  check_vec(beta_a, "beta_a")

  extra <- list(...)
  if (length(extra) && (is.null(names(extra)) || any(names(extra) == ""))) {
    stop("Additional values passed through '...' must be named.", call. = FALSE)
  }

  out <- c(list(beta0_y = beta0_y, beta_x = beta_x, beta_y = beta_y, beta_a = beta_a),
           extra)
  structure(out, class = c("abrm_params", "list"))
}

#' @param x An \code{"abrm_params"} object.
#' @param ... Unused.
#' @rdname abrm_params
#' @return \code{print()} returns \code{x} invisibly.
#' @export
print.abrm_params <- function(x, ...) {
  cat("ABRM true parameters\n")
  show <- function(label, v) {
    if (is.null(v)) {
      cat(sprintf("  %-8s (not supplied)\n", label))
    } else if (is.null(names(v))) {
      cat(sprintf("  %-8s %s  [by position]\n", label, paste(format(v), collapse = ", ")))
    } else {
      cat(sprintf("  %-8s %s  [by name]\n", label,
                  paste0(names(v), " = ", format(v), collapse = ", ")))
    }
  }
  show("beta0_y", x$beta0_y)
  show("beta_x", x$beta_x)
  show("beta_y", x$beta_y)
  show("beta_a", x$beta_a)
  invisible(x)
}
