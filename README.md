# gipsDA

[![Test coverage](https://github.com/AntoniKingston/gipsDA/actions/workflows/test-coverage.yaml/badge.svg?branch=nomult_from_scratch)](https://github.com/AntoniKingston/gipsDA/actions/workflows/test-coverage.yaml)
[![Codecov](https://codecov.io/gh/AntoniKingston/gipsDA/branch/nomult_from_scratch/graph/badge.svg)](https://app.codecov.io/gh/AntoniKingston/gipsDA)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE.md)

`gipsDA` provides linear and quadratic discriminant analysis with covariance
matrices projected onto permutation-invariant structures estimated by
[`gips`](https://github.com/PrzeChoj/gips). The package follows the familiar
`MASS::lda()` and `MASS::qda()` interfaces while replacing their classical
covariance estimates with structured projections.

The package is under active development and is not currently available on
CRAN.

## Models

- `gipslda()` fits LDA using one projected within-class covariance matrix.
- `gipsqda()` fits QDA by projecting each class covariance independently.
- `gipsmultqda()` jointly projects all class covariance matrices, allowing
  them to share an estimated permutation symmetry.

All three functions support formula, matrix, and data-frame interfaces.
Prediction methods provide posterior probabilities and class assignments.
LDA additionally provides discriminant coordinates, coefficients, and
visualization methods.

## Installation

Install the development version from GitHub with
[`pak`](https://pak.r-lib.org/):

```r
install.packages("pak")
pak::pkg_install("AntoniKingston/gipsDA")
```

The required development revision of `gips` is declared in `DESCRIPTION` and
will be installed automatically.

## Quick start

The example below creates a stratified train/test split of the iris data and
fits all three models.

```r
library(gipsDA)

set.seed(42)
train_id <- unlist(
  lapply(split(seq_len(nrow(iris)), iris$Species), sample, size = 35),
  use.names = FALSE
)

train <- iris[train_id, ]
test <- iris[-train_id, ]

lda_fit <- gipslda(Species ~ ., train, optimizer = "BF")
qda_fit <- gipsqda(Species ~ ., train, optimizer = "BF")
joint_qda_fit <- gipsmultqda(Species ~ ., train, optimizer = "BF")

lda_prediction <- predict(lda_fit, test)
qda_prediction <- predict(qda_fit, test)
joint_prediction <- predict(joint_qda_fit, test)

mean(lda_prediction$class == test$Species)
head(lda_prediction$posterior)
```

The equivalent matrix interface separates predictors from class labels:

```r
x <- as.matrix(iris[, 1:4])
grouping <- iris$Species

fit <- gipslda(x, grouping, optimizer = "BF")
predict(fit, x[1:5, ])$class
```

## Projection options

The principal projection arguments are shared across the models:

- `MAP = TRUE` uses the maximum a posteriori permutation.
- `MAP = FALSE` averages projections over retained permutations using their
  posterior probabilities.
- `optimizer = "BF"` performs deterministic brute-force optimization and is
  suitable for smaller problems.
- `optimizer = "MH"` uses Metropolis-Hastings optimization; use `max_iter` to
  control its runtime.

`gipslda()` also accepts `weighted_avg = TRUE`, which constructs the pooled
scatter estimate from a class-proportion-weighted average of class covariance
matrices.

```r
fit <- gipslda(
  Species ~ .,
  iris,
  MAP = FALSE,
  optimizer = "BF",
  weighted_avg = TRUE
)
```

## Prediction rules

LDA supports plug-in, predictive, and debiased prediction:

```r
predict(lda_fit, test, method = "plug-in")
predict(lda_fit, test, method = "predictive")
predict(lda_fit, test, method = "debiased")
```

Both QDA variants additionally support leave-one-out cross-validation when
`newdata` is omitted:

```r
predict(qda_fit, method = "looCV")
predict(joint_qda_fit, method = "looCV")
```

## LDA diagnostics and visualization

```r
coef(lda_fit)
plot(lda_fit)
pairs(lda_fit, type = "std")
pairs(lda_fit, type = "trellis")
```

## Development

The test suite is divided into focused unit tests and broader integration
validations:

```text
tests/testthat/unit/
tests/testthat/manual/
```

Run the complete suite and package checks with:

```r
devtools::test()
devtools::check()
```

Tests run automatically on pushes to the development branch through GitHub
Actions.

## References

Chojecki, A., et al. (2025). *Learning Permutation Symmetry of a Gaussian
Vector with gips in R*. Journal of Statistical Software, 112(7), 1–38.
[doi:10.18637/jss.v112.i07](https://doi.org/10.18637/jss.v112.i07)

The discriminant-analysis implementations are based on the interfaces and
algorithms in:

Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics with S*.
Fourth edition. Springer.

## License

`gipsDA` is licensed under GPL (>= 3).
