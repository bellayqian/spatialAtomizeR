## Resubmission

This is a resubmission improving package quality. Changes include:

* Fixed Editor's comments
* Fixed `print()` statement causing debug output in console
* Fixed unable to use S3 plot() function bug
* Marked 14 internal helper functions as `@keywords internal` (no longer exported)
* Fixed vignette errors: removed non-existent function parameters and corrected signatures
* Improved S3 method output formatting
* Enhanced documentation with working code examples
* Added vcov() method for variance-covariance matrices
* Fixed coefficient naming for beta_x and beta_y parameters

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

## Previous submission (v0.2.4)

Addressed package environment issue with NIMBLE model compilation by refactoring custom distributions and moving nimble from Imports to Depends.