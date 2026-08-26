# Quadratic Discriminant Analysis with multiple gips-projected covariances

Quadratic Discriminant Analysis (QDA) in which each class covariance
matrix participates in one joint *gips* projection, allowing the
covariance matrices to share an estimated permutation symmetry.

## Usage

``` r
gipsmultqda(x, ...)

# S3 method for class 'formula'
gipsmultqda(formula, data, ..., subset, na.action)

# Default S3 method
gipsmultqda(x, grouping, prior = proportions,
  nu = 5, MAP = TRUE, optimizer = NULL, max_iter = NULL,
  show_progress_bar = FALSE, ...)

# S3 method for class 'data.frame'
gipsmultqda(x, ...)

# S3 method for class 'matrix'
gipsmultqda(x, grouping, ..., subset, na.action)
```

## Arguments

- x:

  (required if no formula is given as the principal argument) a matrix
  or data frame containing the explanatory variables.

- ...:

  Arguments passed to or from other methods.

- formula:

  A formula of the form `groups ~ x1 + x2 + ...`. The response is the
  grouping factor and the right-hand side specifies the (non-factor)
  discriminators.

- data:

  An optional data frame, list or environment from which variables
  specified in `formula` are preferentially taken.

- grouping:

  (required if no formula is given) a factor specifying the class for
  each observation.

- prior:

  Prior probabilities of class membership. If omitted, training-set
  class proportions are used. Supplied values must sum to one.

- nu:

  Reserved for compatibility with related QDA interfaces; it is not
  currently used by the covariance projection.

- MAP:

  Logical; if `TRUE`, a maximum a posteriori covariance projection is
  used. If `FALSE`, projections are averaged using retained posterior
  permutation probabilities.

- optimizer:

  Character string specifying the optimization method used for
  covariance projection. If `NULL`, a default choice is made based on
  the problem dimension.

- max_iter:

  Maximum number of iterations for stochastic optimizers.

- show_progress_bar:

  Logical; if `TRUE`, display the progress bar from the underlying gips
  optimizer. Defaults to `FALSE`.

- subset:

  An index vector specifying the cases to be used in the training
  sample. (NOTE: must be named.)

- na.action:

  A function specifying the action to be taken if `NA`s are found.

## Value

An object of class `"gipsmultqda"` containing:

- `prior`: prior probabilities of the groups

- `counts`: number of observations per group

- `means`: group means

- `scaling`: array of group-specific scaling matrices derived from the
  projected covariance matrices

- `ldet`: log-determinants of the projected covariance matrices

- `lev`: class labels

- `N`: total number of observations

- `optimization_info`: estimated probabilities of retained permutations
  returned by the joint gips optimization

- `call`: the matched call

- Formula fits additionally contain `terms`, `contrasts`, `xlevels`, and
  any recorded `na.action`.

## Details

This function is a modification of
[`qda`](https://rdrr.io/pkg/MASS/man/qda.html) in which the
class-specific covariance matrices are jointly projected to improve
numerical stability and exploit shared symmetry assumptions.

In contrast to classical QDA, which estimates each class covariance
matrix independently, `gipsmultqda` performs a joint projection of all
class covariance matrices using
[`gips`](https://przechoj.github.io/gips/reference/gips.html). This
allows the incorporation of shared permutation symmetries and can
improve classification performance in high-dimensional or small-sample
regimes.

Several classification rules are available via
[`predict.gipsmultqda`](https://antonikingston.github.io/gipsDA/reference/predict.gipsmultqda.md),
including plug-in, predictive, debiased, and leave-one-out
cross-validation.

## Note

This function is not a drop-in replacement for
[`qda`](https://rdrr.io/pkg/MASS/man/qda.html). The covariance
estimation, returned object, and classification rules differ
substantially.

The theoretical background and details of covariance projection are
documented by
[`gips`](https://przechoj.github.io/gips/reference/gips.html) and
Chojecki et al. (2025).

## See also

[`qda`](https://rdrr.io/pkg/MASS/man/qda.html),
[`predict.gipsmultqda`](https://antonikingston.github.io/gipsDA/reference/predict.gipsmultqda.md),
[`gipsqda`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md),
[`gipslda`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)

## Examples

``` r
tr <- sample(1:50, 25)
train <- rbind(iris3[tr, , 1], iris3[tr, , 2], iris3[tr, , 3])
test <- rbind(iris3[-tr, , 1], iris3[-tr, , 2], iris3[-tr, , 3])
cl <- factor(c(rep("s", 25), rep("c", 25), rep("v", 25)))
z <- gipsmultqda(train, cl)
predict(z, test)$class
#>  [1] s s s s s s s s s s s s s s s s s s s s s s s s s c c c c c c c c c c c c c
#> [39] c c c c c c c c c c c c v v v v v v v v v v v c v v v v v v v c v v v v v
#> Levels: c s v
```
