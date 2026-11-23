# Internal helpers shared across the package. None of these are exported.

# ---------------------------------------------------------------------------
# Distribution names
# ---------------------------------------------------------------------------

# Accepted spellings for variable distributions, mapped to the canonical
# short form used throughout the package ("norm", "pois", "binom").
.dist_lookup <- c(
  norm = "norm", normal = "norm", gaussian = "norm",
  pois = "pois", poisson = "pois",
  binom = "binom", binomial = "binom"
)

# Integer codes used inside the NIMBLE model for the outcome (vartype_y).
.dist_code <- c(norm = 1L, pois = 2L, binom = 3L)

#' Normalise distribution names to "norm" / "pois" / "binom"
#' @param x Character vector of distribution names.
#' @param arg Argument name, used in error messages.
#' @return Character vector of canonical names, same length as \code{x}.
#' @keywords internal
#' @noRd
.normalize_dist <- function(x, arg) {
  if (is.null(x) || length(x) == 0L) {
    stop(sprintf("'%s' must be supplied.", arg), call. = FALSE)
  }
  if (!is.character(x)) {
    stop(sprintf(
      "'%s' must be a character vector with elements in \"norm\", \"pois\", \"binom\".",
      arg), call. = FALSE)
  }
  out <- .dist_lookup[tolower(x)]
  bad <- is.na(out)
  if (any(bad)) {
    stop(sprintf(
      "'%s' contains unrecognised distribution(s): %s. Use \"norm\", \"pois\" or \"binom\".",
      arg, paste0("\"", x[bad], "\"", collapse = ", ")), call. = FALSE)
  }
  unname(out)
}

# ---------------------------------------------------------------------------
# Parameter naming
# ---------------------------------------------------------------------------

#' Split NIMBLE's stacked beta_y vector into beta_x / beta_y / beta_a
#'
#' Inside the model the outcome coefficients are a single vector
#' \code{beta_y[1:(p_x + p_y + p_a)]}. This maps those names to
#' \code{beta_x[1..p_x]}, \code{beta_y[1..p_y]} and \code{beta_a[1..p_a]}.
#' Any other names are returned unchanged.
#' @keywords internal
#' @noRd
.rename_beta_params <- function(nms, p_x, p_y, p_a = 0L) {
  if (is.null(nms)) return(nms)
  m <- regmatches(nms, regexec("^beta_y\\[(\\d+)\\]$", nms))
  k <- vapply(m, function(z) if (length(z) == 2L) as.integer(z[2]) else NA_integer_,
              integer(1))
  hit <- !is.na(k)
  kk <- k[hit]
  nms[hit] <- ifelse(
    kk <= p_x, sprintf("beta_x[%d]", kk),
    ifelse(kk <= p_x + p_y, sprintf("beta_y[%d]", kk - p_x),
           sprintf("beta_a[%d]", kk - p_x - p_y))
  )
  nms
}

#' Map raw parameter names (beta_x[1], ...) to user-facing labels
#'
#' Uses the \code{param_labels} stored on the fitted object. Objects created
#' by older versions of the package do not carry labels, so a fallback based
#' on the legacy \code{covariate_x_<i>} naming is used for them.
#' @keywords internal
#' @noRd
.abrm_label <- function(object, raw) {
  lab <- object$param_labels
  if (!is.null(lab)) {
    out <- unname(lab[raw])
    out[is.na(out)] <- raw[is.na(out)]
    return(out)
  }
  out <- gsub("^beta_x\\[(\\d+)\\]$", "covariate_x_\\1", raw)
  out <- gsub("^beta_y\\[(\\d+)\\]$", "covariate_y_\\1", out)
  out <- gsub("^beta_a\\[(\\d+)\\]$", "covariate_a_\\1", out)
  out[out %in% c("beta_0_y", "beta0_y")] <- "intercept"
  out
}

#' Pool MCMC chains of an abrm object into one matrix
#' @keywords internal
#' @noRd
.abrm_sample_matrix <- function(object) {
  samples <- object$mcmc_results$samples
  if (is.null(samples)) {
    stop("MCMC samples not found in the abrm object.", call. = FALSE)
  }
  if (is.matrix(samples)) return(samples)
  if (is.list(samples) && all(vapply(samples, is.matrix, logical(1)))) {
    return(do.call(rbind, samples))
  }
  as.matrix(samples)
}

#' Ordered raw names of the regression coefficients in an abrm object
#'
#' Returns intercept, then X-grid, Y-grid and atom-level coefficients, each in
#' index order. If \code{parm} is given, the result is restricted to the
#' parameters whose raw name (e.g. \code{"beta_x[1]"}) or label (e.g. the
#' covariate name) is in \code{parm}; unknown entries raise an error.
#' @keywords internal
#' @noRd
.abrm_beta_names <- function(object, available, parm = NULL) {
  by_index <- function(pattern) {
    v <- grep(pattern, available, value = TRUE)
    v[order(as.integer(gsub("\\D", "", v)))]
  }
  raw <- c(
    grep("^beta_0_y$", available, value = TRUE),
    by_index("^beta_x\\["),
    by_index("^beta_y\\["),
    by_index("^beta_a\\[")
  )
  if (is.null(parm)) return(raw)
  labs <- .abrm_label(object, raw)
  unknown <- setdiff(parm, c(raw, labs))
  if (length(unknown) > 0L) {
    stop("Unknown parameter(s) in 'parm': ", paste(unknown, collapse = ", "),
         ". Available: ", paste(labs, collapse = ", "), call. = FALSE)
  }
  raw[raw %in% parm | labs %in% parm]
}

# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------

#' Sum an atom-level variable within each grid cell (cells numbered 1..n)
#' @keywords internal
#' @noRd
.sum_by_id <- function(values, ids, n) {
  vapply(seq_len(n), function(j) sum(values[ids == j]), numeric(1))
}

#' Store whole-number numeric IDs as integers
#'
#' Doubles that hold whole numbers within the integer range are converted to
#' integer so that IDs compare exactly and match the package's other index
#' variables. Character IDs, non-whole values and values too large for an
#' integer (e.g. 11-digit census-tract GEOIDs) are returned unchanged.
#' @keywords internal
#' @noRd
.as_integer_id <- function(x) {
  if (is.double(x) && !anyNA(x) && all(x == round(x)) &&
      all(abs(x) <= .Machine$integer.max)) {
    return(as.integer(x))
  }
  x
}
