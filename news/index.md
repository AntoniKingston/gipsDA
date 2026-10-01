# Changelog

## gipsDA 1.0.0

### Bug fixes

- Fixed prediction from formula fits when fitting-only arguments such as
  `show_progress_bar` and `store_probabilities` are supplied.
- Fixed prediction from matrix fits using `na.action`.
- Improved handling of edge cases involving one-predictor models and
  unused factor levels.
- Added validation for `method = "debiased"` in small-sample settings.

### Documentation

- Clarified the interpretation of selected permutation structures in
  [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
  versus QDA models.
- Updated the manual reconstruction vignette to match the current
  covariance scaling and class-specific QDA sample-size handling.
- Documented fitting metadata stored in `fit_info` and summary objects.

### Maintenance

- Updated the minimum supported R version.
- Improved build ignores for local coverage files, plots, and generated
  artifacts.
