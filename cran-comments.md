## Resubmission

This is a patch resubmission fixing a bug reported by a CRAN reviewer.

* Fixed: atom ID column names produced by `raster::intersect()` + `sf::st_as_sf()` 
  vary by sf version ("ID"/"ID.1" in older sf, "ID_1"/"ID_2" in sf >= 1.0). 
  The renaming step in `simulate_misaligned_data()` and the case example vignette 
  now handles both conventions using `grep("^ID", ...)`.
  
## Previous submission (v0.2.7)

Addressing reviewer comments on package quality,
documentation, and replication material. Changes include:

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
