# spatialAtomizeR 0.2.5

## Major Changes

* Fixed `print()` statement in `prepare_nimble_inputs()` that caused debug output (`[1] "\nCreating..."`) to appear in console
* Marked 14 internal helper functions as `@keywords internal` to prevent export:
  - NIMBLE distribution functions: `biasedUrn_rmfnc`, `dmfnchypg`, `Rmfnchypg`, `rmfnchypg`, `register_nimble_distributions`
  - Diagnostic functions: `check_mcmc_diagnostics`, `create_diagnostic_plots`, `print_convergence_summary`, `create_comparison_plots`, `create_summary_statistics`, `create_sensitivity_summary_plots`
  - Preparation functions: `prepare_spatial_bookkeeping`, `prepare_adjacency_matrices`, `prepare_nimble_inputs`

## Documentation

* **Major vignette overhaul** - Fixed critical errors and improved usability:
  - Removed non-existent parameters (`res1`, `res2`) from `simulate_misaligned_data()` examples
  - Corrected `run_abrm()` function signature (now uses `gridx`, `gridy`, `atoms` instead of `sim_data`)
  - Removed "Advanced" section showing internal functions
  - Added Example 1: Complete simulated data workflow with S3 method demonstrations
  - Added Example 2: Real-world Utah counties and school districts analysis (matches manuscript)
  - Enhanced troubleshooting section with common errors and solutions
* Improved `summary.abrm()` output formatting (cleaner column names, removed redundancy)
* All vignette code examples now tested and verified to run without errors

# spatialAtomizeR 0.2.4

## Major Changes

* **Refactored NIMBLE model architecture**: Custom distributions now defined as top-level exported functions within package namespace (addressing package environment issue from December 4, 2024 submission)
* **Moved `nimble` from Imports to Depends**: Ensures internal compilation tools are correctly available on search path without global assignments

## Bug Fixes

* Fixed package environment issue that prevented proper NIMBLE model compilation

## Notes

* NIMBLE custom distributions (`dmfnchypg`, `rmfnchypg`, `Rmfnchypg`) now properly registered and accessible
