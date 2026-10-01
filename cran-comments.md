## R CMD check results

0 errors | 0 warnings | 1 note

* checking HTML version of manual ... NOTE
  Skipping HTML validation because the local HTML Tidy installation is not recent
  enough. This is a local toolchain issue.

## Test environments

* local macOS, R 4.5.1
* R CMD check --as-cran
* GitHub Actions:
  * macOS-latest, R release
  * ubuntu-latest, R release
  * windows-latest, R release
* GitHub Actions pkgdown workflow
* GitHub Actions test coverage workflow

## Release summary

This is a major update from gipsDA 0.1.2 to 1.0.0.

Main changes include:

* fixed prediction edge cases for formula and matrix fits,
* improved handling of one-predictor models and unused factor levels,
* added validation for small-sample debiased prediction,
* clarified interpretation of selected permutation structures,
* updated mathematical vignettes to match the implementation,
* documented fitting metadata stored in `fit_info` and summary objects,
* updated package metadata and release files.