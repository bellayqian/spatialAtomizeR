## R CMD check results

0 errors | 0 warnings | 0 note

## Resubmission

This is a resubmission addressing the package environment issue found in the previous submission (December 4, 2025).

## Changes made:

REFACTORED the Nimble model architecture: custom distributions are now defined as top-level exported functions within the package namespace.

MOVED nimble from Imports to Depends to ensure internal compilation tools are correctly available on the search path without needing global assignments.

## Notes

If there is 1 note about timestamps: The note regarding "future file timestamps" is due to local system clock differences and does not affect package functionality.