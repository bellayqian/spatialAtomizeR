# ─────────────────────────────────────────────────────────────────────────────
# test_0.3.0.R
# Tests the updated spatialAtomizeR 0.3.0 with population denominator changes.
# Runs ABRM under both population formats (single and per-variable).
#
# Usage:
#   Rscript test_0.3.0.R
#   # or interactively: source("test_0.3.0.R")
# ─────────────────────────────────────────────────────────────────────────────

PKG_PATH <- "~/Desktop/有事/Harvard/Research/1_RachelNethery_HealthEquity/spatialAtomizeR_0703RCN_1118"

cat("Loading package from source...\n")
if ("package:spatialAtomizeR" %in% search()) {
  detach("package:spatialAtomizeR", unload = TRUE, character.only = TRUE)
}
if ("spatialAtomizeR" %in% loadedNamespaces()) {
  unloadNamespace("spatialAtomizeR")
}
devtools::load_all(PKG_PATH, quiet = TRUE)

# ── 1. Option 1: Single population (backward compatibility) ─────────────────
cat("\n=== TEST 1: Single population column (Option 1) ===\n")

sim_data <- simulate_misaligned_data(
  seed               = 42,
  dist_covariates_x  = "poisson",
  dist_covariates_y  = "normal",
  dist_y             = "poisson",
  x_intercepts       = -1,
  y_intercepts       =  3,
  beta0_y            = -2,
  beta_x             =  0.5,
  beta_y             = -0.2
)
print(sim_data)

cat("  Atoms population columns: ")
pop_cols <- grep("^population", names(sim_data$atoms), value = TRUE)
cat(paste(pop_cols, collapse = ", "), "\n")
stopifnot("Option 1 should have single 'population' column" = identical(pop_cols, "population"))

cat(sprintf("  True beta_x : %.2f\n", sim_data$true_params$beta_x))
cat(sprintf("  True beta_y : %.2f\n", sim_data$true_params$beta_y))

model_code <- get_abrm_model()

cat("  Running ABRM (niter=10000, nburnin=5000)...\n")
t1 <- system.time({
  fit1 <- run_abrm(
    gridx       = sim_data$gridx,
    gridy       = sim_data$gridy,
    atoms       = sim_data$atoms,
    model_code  = model_code,
    true_params = sim_data$true_params,
    x_types  = "pois",
    y_types     = "norm",
    dist_y      = "pois",
    niter       = 10000,
    nburnin     = 5000,
    nchains     = 2,
    thin        = 5,
    seed        = 1,
    compute_waic = FALSE
  )
})

cat(sprintf("  Fit time: %.1f s\n", t1["elapsed"]))
print(fit1)

cat("\n  Parameter estimates:\n")
print(fit1$parameter_estimates[, c("variable", "estimated_beta", "ci_lower", "ci_upper", "true_beta", "bias")])

# ── 2. Option 2: Per-variable populations ────────────────────────────────────
cat("\n=== TEST 2: Per-variable population columns (Option 2) ===\n")

sim_data_pv <- simulate_misaligned_data(
  seed               = 42,
  dist_covariates_x  = "poisson",
  dist_covariates_y  = "normal",
  dist_y             = "poisson",
  x_intercepts       = -1,
  y_intercepts       =  3,
  beta0_y            = -2,
  beta_x             =  0.5,
  beta_y             = -0.2,
  diff_pops_per_var  = TRUE
)
print(sim_data_pv)

cat("  Atoms population columns: ")
pop_cols_pv <- grep("^population", names(sim_data_pv$atoms), value = TRUE)
cat(paste(pop_cols_pv, collapse = ", "), "\n")
stopifnot("Option 2 should NOT have single 'population' column" = !"population" %in% pop_cols_pv)
stopifnot("Option 2 should have population_covariate_x_1" = "population_covariate_x_1" %in% pop_cols_pv)
stopifnot("Option 2 should have population_y" = "population_y" %in% pop_cols_pv)

cat("  Running ABRM with per-variable populations...\n")
t2 <- system.time({
  fit2 <- run_abrm(
    gridx       = sim_data_pv$gridx,
    gridy       = sim_data_pv$gridy,
    atoms       = sim_data_pv$atoms,
    model_code  = model_code,
    true_params = sim_data_pv$true_params,
    x_types  = "pois",
    y_types     = "norm",
    dist_y      = "pois",
    niter       = 10000,
    nburnin     = 5000,
    nchains     = 2,
    thin        = 5,
    seed        = 1,
    compute_waic = FALSE
  )
})

cat(sprintf("  Fit time: %.1f s\n", t2["elapsed"]))
print(fit2)

cat("\n  Parameter estimates:\n")
print(fit2$parameter_estimates[, c("variable", "estimated_beta", "ci_lower", "ci_upper", "true_beta", "bias")])

# ── 3. Error handling: both formats present should fail ──────────────────────
cat("\n=== TEST 3: Conflict detection (both formats present) ===\n")

sim_conflict <- sim_data
sim_conflict$atoms$population_covariate_x_1 <- sim_conflict$atoms$population
err_caught <- tryCatch({
  run_abrm(
    gridx = sim_conflict$gridx, gridy = sim_conflict$gridy,
    atoms = sim_conflict$atoms, model_code = model_code,
    x_types = "pois", y_types = "norm", dist_y = "pois",
    niter = 100, nburnin = 50, nchains = 1
  )
  FALSE
}, error = function(e) {
  cat("  Expected error caught:", conditionMessage(e), "\n")
  TRUE
})
stopifnot("Conflict detection should trigger an error" = err_caught)
cat("  PASSED\n")

# ── 4. Multiple covariates with per-variable populations ─────────────────────
cat("\n=== TEST 4: Multiple covariates, per-variable populations ===\n")

sim_multi <- simulate_misaligned_data(
  seed               = 99,
  dist_covariates_x  = c("poisson", "normal"),
  dist_covariates_y  = c("normal", "poisson"),
  dist_y             = "poisson",
  x_intercepts       = c(-1, 2),
  y_intercepts       = c(3, -1),
  beta0_y            = -1,
  beta_x             = c(0.3, -0.1),
  beta_y             = c(-0.2, 0.4),
  diff_pops_per_var  = TRUE
)

cat("  Atoms population columns: ")
pop_cols_multi <- grep("^population", names(sim_multi$atoms), value = TRUE)
cat(paste(pop_cols_multi, collapse = ", "), "\n")
stopifnot(length(pop_cols_multi) == 5)  # x1, x2, y1, y2, population_y

cat("  Running ABRM with multiple covariates + per-variable pops...\n")
t4 <- system.time({
  fit4 <- run_abrm(
    gridx       = sim_multi$gridx,
    gridy       = sim_multi$gridy,
    atoms       = sim_multi$atoms,
    model_code  = model_code,
    true_params = sim_multi$true_params,
    x_types  = c("pois", "norm"),
    y_types     = c("norm", "pois"),
    dist_y      = "pois",
    niter       = 10000,
    nburnin     = 5000,
    nchains     = 2,
    thin        = 5,
    seed        = 1,
    compute_waic = FALSE
  )
})

cat(sprintf("  Fit time: %.1f s\n", t4["elapsed"]))
print(fit4)

cat("\n  Parameter estimates:\n")
print(fit4$parameter_estimates[, c("variable", "estimated_beta", "ci_lower", "ci_upper", "true_beta", "bias")])

# ── 5. Atom-level covariates (covariate_a_*) ─────────────────────────────────
cat("\n=== TEST 5: Atom-level covariates (covariate_a_*) ===\n")

sim_data_a <- simulate_misaligned_data(
  seed               = 42,
  dist_covariates_x  = "poisson",
  dist_covariates_y  = "normal",
  dist_y             = "poisson",
  x_intercepts       = -1,
  y_intercepts       =  3,
  beta0_y            = -2,
  beta_x             =  0.5,
  beta_y             = -0.2
)

set.seed(42)
n_atoms <- nrow(sim_data_a$atoms)
sim_data_a$atoms$covariate_a_1 <- rnorm(n_atoms, 0, 1)
sim_data_a$atoms$covariate_a_2 <- rnorm(n_atoms, 0, 1)

cat("  Atom-level covariate columns: ")
a_cols <- grep("^covariate_a_", names(sim_data_a$atoms), value = TRUE)
cat(paste(a_cols, collapse = ", "), "\n")
stopifnot("Should have 2 atom-level covariates" = length(a_cols) == 2)

cat("  Running ABRM with atom-level covariates...\n")
t5 <- system.time({
  fit5 <- run_abrm(
    gridx       = sim_data_a$gridx,
    gridy       = sim_data_a$gridy,
    atoms       = sim_data_a$atoms,
    model_code  = model_code,
    true_params = sim_data_a$true_params,
    x_types  = "pois",
    y_types     = "norm",
    dist_y      = "pois",
    niter       = 10000,
    nburnin     = 5000,
    nchains     = 2,
    thin        = 5,
    seed        = 1,
    compute_waic = FALSE
  )
})

cat(sprintf("  Fit time: %.1f s\n", t5["elapsed"]))
print(fit5)

cat("\n  Parameter estimates:\n")
print(fit5$parameter_estimates[, c("variable", "estimated_beta", "ci_lower", "ci_upper")])

# ── 6. Both features: atom-level covariates + per-variable populations ───────
cat("\n=== TEST 6: Atom-level covariates + per-variable populations ===\n")

sim_data_both <- simulate_misaligned_data(
  seed               = 42,
  dist_covariates_x  = "poisson",
  dist_covariates_y  = "normal",
  dist_y             = "poisson",
  x_intercepts       = -1,
  y_intercepts       =  3,
  beta0_y            = -2,
  beta_x             =  0.5,
  beta_y             = -0.2,
  diff_pops_per_var  = TRUE
)

set.seed(42)
n_atoms_both <- nrow(sim_data_both$atoms)
sim_data_both$atoms$covariate_a_1 <- rnorm(n_atoms_both, 0, 1)

cat("  Population columns: ")
cat(paste(grep("^population", names(sim_data_both$atoms), value = TRUE), collapse = ", "), "\n")
cat("  Atom-level covariate columns: ")
cat(paste(grep("^covariate_a_", names(sim_data_both$atoms), value = TRUE), collapse = ", "), "\n")

cat("  Running ABRM with both features...\n")
t6 <- system.time({
  fit6 <- run_abrm(
    gridx       = sim_data_both$gridx,
    gridy       = sim_data_both$gridy,
    atoms       = sim_data_both$atoms,
    model_code  = model_code,
    true_params = sim_data_both$true_params,
    x_types  = "pois",
    y_types     = "norm",
    dist_y      = "pois",
    niter       = 10000,
    nburnin     = 5000,
    nchains     = 2,
    thin        = 5,
    seed        = 1,
    compute_waic = FALSE
  )
})

cat(sprintf("  Fit time: %.1f s\n", t6["elapsed"]))
print(fit6)

cat("\n  Parameter estimates:\n")
print(fit6$parameter_estimates[, c("variable", "estimated_beta", "ci_lower", "ci_upper")])

# ── 7. Flexible column names (no fixed 'covariate_x_1' / 'y' / 'ID' naming) ──
cat("\n=== TEST 7: Custom column names, character IDs, named true_params ===\n")

# Same data as Test 1, but every column the model uses gets a user-chosen name.
gx <- sim_data$gridx; gy <- sim_data$gridy; at <- sim_data$atoms
names(gx)[names(gx) == "covariate_x_1"] <- "income"
names(gy)[names(gy) == "covariate_y_1"] <- "pct_over65"
names(gy)[names(gy) == "y"]             <- "deaths"
names(at)[names(at) == "population"]    <- "pop"

# Character IDs under non-default column names
gx$tract_id  <- paste0("T", gx$ID)
gy$county_id <- paste0("C", gy$ID)
at$tract_id  <- paste0("T", at$ID_x)
at$county_id <- paste0("C", at$ID_y)

# (the grids no longer contain any default-named model column; atoms keep their
#  atom-level simulated columns, which run_abrm() does not read)
stopifnot("grids must not contain the default model column names" =
            !any(c("covariate_x_1", "covariate_y_1", "y") %in% c(names(gx), names(gy))),
          "atoms must not contain the default population column" = !"population" %in% names(at))

cat("  Running ABRM with custom names (same data and seeds as a default-named fit)...\n")
set.seed(7)
fit7 <- run_abrm(
  gridx = gx, gridy = gy, atoms = at, model_code = model_code,
  x_vars = "income", y_vars = "pct_over65", outcome = "deaths",
  id_x = "tract_id", id_y = "county_id",
  atom_id_x = "tract_id", atom_id_y = "county_id",
  population = "pop",
  x_types = "pois", y_types = "norm", dist_y = "pois",
  true_params = abrm_params(beta0_y = -2, beta_x = c(income = 0.5),
                            beta_y = c(pct_over65 = -0.2)),
  niter = 3000, nburnin = 1500, nchains = 2, thin = 5, seed = 1
)

set.seed(7)
fit7_default <- run_abrm(
  gridx = sim_data$gridx, gridy = sim_data$gridy, atoms = sim_data$atoms,
  model_code = model_code, x_types = "pois", y_types = "norm", dist_y = "pois",
  true_params = sim_data$true_params,
  niter = 3000, nburnin = 1500, nchains = 2, thin = 5, seed = 1
)

# (a) results are labelled with the user's names
stopifnot("labels should be the user's covariate names" =
            identical(names(coef(fit7)), c("intercept", "income", "pct_over65")))
stopifnot(identical(rownames(confint(fit7)), c("intercept", "income", "pct_over65")))
stopifnot(identical(mcmc_diagnostics(fit7)$variable, c("intercept", "income", "pct_over65")))
stopifnot("model_info should record the names" =
            identical(fit7$model_info$x_vars, "income") &&
            identical(fit7$model_info$outcome, "deaths"))
# (b) character IDs come back as the user's IDs
stopifnot("y_grid_ids should be the user's county IDs" =
            all(fit7$y_grid_ids %in% gy$county_id))
# (c) renaming columns must not change the inference
stopifnot("renaming columns changed the posterior means" =
            isTRUE(all.equal(unname(coef(fit7)), unname(coef(fit7_default)))))
# (d) true_params given by name were matched to the right covariates
stopifnot(isTRUE(all.equal(fit7$parameter_estimates$true_beta, c(-2, 0.5, -0.2))))
# (e) bad names give clear errors
msg <- tryCatch(run_abrm(gx, gy, at, model_code, x_vars = "nope", y_vars = "pct_over65",
                         outcome = "deaths", id_x = "tract_id", id_y = "county_id",
                         atom_id_x = "tract_id", atom_id_y = "county_id", population = "pop",
                         x_types = "pois", y_types = "norm", dist_y = "pois"),
                error = conditionMessage)
cat("  Missing-column error:", msg, "\n")
stopifnot(grepl("nope", msg))
cat("  PASSED (names, character IDs, identical posterior means to the default-named fit)\n")
print(fit7$parameter_estimates[, c("variable", "parameter", "estimated_beta", "true_beta")])

# ── 8. New helpers: plots, diagnostics, atoms (re-use earlier fits; no refit) ─
cat("\n=== TEST 8: ggplot returns, diagnostics accessor, create_atoms ===\n")
stopifnot(inherits(plot_trace(fit1),   "ggplot"))
stopifnot(inherits(plot_density(fit1), "ggplot"))
stopifnot(inherits(plot(fit1, type = "density"), "ggplot"))
stopifnot(inherits(plot(sim_data), "ggplot"), inherits(plot(sim_data, which = "atoms"), "ggplot"))
d1 <- mcmc_diagnostics(fit1)
stopifnot(all(c("variable", "parameter", "rhat", "ess", "converged") %in% names(d1)),
          nrow(d1) == 3)
no_atoms <- sim_data; no_atoms$atoms <- NULL
stopifnot("plot(which='atoms') must error without atoms" =
            inherits(try(plot(no_atoms, which = "atoms"), silent = TRUE), "try-error"))
at_new <- create_atoms(sim_data$gridy, sim_data$gridx)
stopifnot(nrow(at_new) == nrow(sim_data$atoms),
          is.integer(at_new$ID_x), is.integer(at_new$ID_y),
          is.integer(sim_data$atoms$ID_x), is.integer(sim_data$gridx$ID))
cat("  PASSED\n")

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n=== ALL TESTS COMPLETE ===\n")
cat("Test 1: Single population (backward compat) — ")
cat(ifelse(all(fit1$parameter_estimates$within_ci), "all CIs cover truth\n", "some CIs miss\n"))
cat("Test 2: Per-variable population — ")
cat(ifelse(all(fit2$parameter_estimates$within_ci), "all CIs cover truth\n", "some CIs miss\n"))
cat("Test 3: Conflict detection — PASSED\n")
cat("Test 4: Multi-covariate per-variable — ")
cat(ifelse(all(fit4$parameter_estimates$within_ci), "all CIs cover truth\n", "some CIs miss\n"))
cat("Test 5: Atom-level covariates — ran successfully\n")
cat("Test 6: Atom-level covariates + per-variable pops — ran successfully\n")
cat("Test 7: Custom column names + character IDs — PASSED\n")
cat("Test 8: ggplot returns, diagnostics accessor, create_atoms — PASSED\n")
