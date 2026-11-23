# spatialAtomizeR 0.3.0

## Breaking changes

* `run_abrm()`: the `norm_idx_x`, `pois_idx_x`, `binom_idx_x`, `norm_idx_y`,
  `pois_idx_y`, `binom_idx_y` arguments are replaced by `x_types` and
  `y_types`, character vectors with one element per covariate taking values
  `"norm"`, `"pois"` or `"binom"`. The old arguments still work but are
  deprecated (with a warning).
* `run_abrm()`: `dist_y` is now a string (`"norm"`, `"pois"` or `"binom"`).
  Numeric codes 1/2/3 still work but are deprecated.
* `plot.misaligned_data()` and `plot.abrm()` now return `ggplot` objects
  instead of printing and invisibly returning the input.
* `plot.abrm()` gains a `type` argument (`"trace"` or `"density"`) and draws
  one kind of plot; see the new `plot_trace()` and `plot_density()`.
  Diagnostic ggplot objects are no longer stored inside the fitted object.
* `coef()`, `vcov()` and `confint()` now consistently label coefficients with
  covariate names, and the intercept as `"intercept"`.
* Requires R >= 4.1.0 (native pipe). `raster`, `MASS` and `reshape2` are no
  longer dependencies.
* `gen_correlated_spat()` uses a different (Cholesky-based) sampler, so
  simulated data for a given seed differ from earlier versions.

## New features

* Covariates no longer need fixed names: `run_abrm()` gains `x_vars`,
  `y_vars`, `outcome`, `atom_vars`, `id_x`, `id_y`, `atom_id_x`, `atom_id_y`
  and `population` arguments. IDs may be character (e.g. GEOIDs). The old
  `covariate_x_1`, ... naming remains the default.
* New `create_atoms()` builds atoms with `sf::st_intersection()`;
  `simulate_misaligned_data()` uses it instead of `raster::intersect()`.
* New `mcmc_diagnostics()` accessor for R-hat and effective sample sizes.
* `gen_correlated_spat()` now uses `var_spat` and `verify`, accepts a full
  correlation matrix, validates its inputs, and avoids inverting the
  Kronecker-product precision matrix.
* New `abrm_params()` constructor for the `true_params` argument of
  `run_abrm()`. It validates its inputs, vectors can be unnamed (matched by
  position) or named (matched by covariate name), and atom-level coefficients
  (`beta_a`) are supported. `simulate_misaligned_data()` now returns its
  `true_params` as an `abrm_params` object. The format is documented in
  `?abrm_params` and `?run_abrm`.

## Bug fixes and internals

* `plot(x, which = "atoms")` errors informatively when atoms are missing.
* ID and index variables (grid IDs, atom IDs, NIMBLE index constants) are
  integers. `create_atoms()` stores whole-number numeric IDs as integers
  (IDs too large for an integer, such as 11-digit GEOIDs, and character IDs
  are left unchanged), and `run_abrm()` no longer requires numeric IDs: any
  unique IDs are matched and mapped to integer indices internally.
* README, vignette, `data-raw/vignette_cache.R` and `test_0.2.9.R` use the new
  interface.
* `biasedUrn_dmfnc()` is exported so models with binomial variables work
  from an installed package (previously only under `devtools::load_all()`).
* Single-chain fits no longer fail; R-hat is reported as `NA`.
* Atom-level coefficients were dropped from `parameter_estimates` when
  `true_params` was supplied; they are now always included.
* Native pipe `|>` replaces `%>%` throughout.

# spatialAtomizeR 0.2.8

This is a resubmission (v0.2.7) addressing reviewer comments on package quality,
documentation, and replication material. Changes include:

* Fixed atom ID column renaming to handle both older sf (which produces 
  `ID`/`ID.1` after `raster::intersect()`) and newer sf >= 1.0 (which 
  produces `ID_1`/`ID_2`). Affected `simulate_misaligned_data()` and 
  the case example vignette (#issue number if you have one).
* Added `\examples{}` blocks to all exported functions; slow MCMC-based functions
  use `\dontrun{}` while lightweight functions (e.g., `simulate_misaligned_data()`,
  `get_abrm_model()`, `gen_correlated_spat()`, and all S3 methods for
  `misaligned_data`) now have fully executable examples
* Changed `vcov.abrm` example from `\dontrun{}` to `\donttest{}` per CRAN policy for examples that are valid but too slow for routine checking
* Fixed `run_abrm()` silently writing a PDF to the user's working directory:
  changed default of `save_plots` from `TRUE` to `FALSE`; diagnostic plots are
  now obtained explicitly via `plot()` on the returned object
* Added `plot.misaligned_data()` S3 method for spatial data exploration
* Added `coef.abrm()`, `confint.abrm()`, `fitted.abrm()`, `waic.abrm()`, `print.waic_abrm()` S3 methods to improve model result inspection without accessing internal object structure directly
* Updated vignette to set `eval = TRUE` for all simulated-data (Example 1)
  code chunks so the vignette executes and displays output during build
* Removed redundant `library(nimble)` call from vignette (nimble is already
  loaded as a dependency of spatialAtomizeR)
* Fixed man page titles to use title case throughout
  (e.g., "Print Method for ABRM Objects" instead of "Print method for abrm objects")

## Test environments

* local: macOS/Windows/Linux (R 4.4.0)
* GitHub Actions:
  - ubuntu-latest (R-release, R-devel)
  - windows-latest (R-release)
  - macOS-latest (R-release)
* win-builder (R-devel)
* R-hub (ubuntu-gcc-release)

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Previous submission (v0.2.6)

* Fixed Editor's comments
* Fixed `print()` statement causing debug output in console
* Fixed unable to use S3 plot() function bug
* Marked 14 internal helper functions as `@keywords internal` (no longer exported)
* Fixed vignette errors: removed non-existent function parameters and corrected signatures
* Improved S3 method output formatting
* Enhanced documentation with working code examples
* Added vcov() method for variance-covariance matrices
* Fixed coefficient naming for beta_x and beta_y parameters
