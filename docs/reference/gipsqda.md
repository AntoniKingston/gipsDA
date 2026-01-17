# Quadratic Discriminant Analysis with gips covariance projection

Quadratic discriminant analysis (QDA) using covariance matrices
projected via the *gips* framework to enforce permutation symmetry and
improve numerical stability.

## Usage

``` r
gipsqda(x, ...)

# S3 method for class 'formula'
gipsqda(formula, data, ..., subset, na.action)

# Default S3 method
gipsqda(x, grouping, prior = proportions,
  nu = 5, MAP = TRUE, optimizer = NULL, max_iter = NULL, ...)

# S3 method for class 'data.frame'
gipsqda(x, ...)

# S3 method for class 'matrix'
gipsqda(x, grouping, ..., subset, na.action)
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

  The prior probabilities of class membership. Must sum to one and have
  length equal to the number of groups.

- nu:

  Degrees of freedom parameter used internally by covariance projection.

- MAP:

  Logical; if `TRUE`, maximum a posteriori covariance projection is
  used.

- optimizer:

  Character string specifying the optimization method used for
  covariance projection. If `NULL`, a default choice depending on the
  problem dimension is used.

- max_iter:

  Maximum number of iterations for stochastic optimizers.

- subset:

  An index vector specifying the cases to be used in the training
  sample. (NOTE: must be named.)

- na.action:

  A function specifying the action to be taken if `NA`s are found.

## Value

An object of class `"gipsqda"` containing the following components:

- `prior`: prior probabilities of the groups

- `counts`: number of observations in each group

- `means`: group means

- `scaling`: group-specific scaling matrices derived from the projected
  covariance matrices

- `ldet`: log-determinants of the projected covariance matrices

- `lev`: class labels

- `N`: total number of observations

- `optimization_info`: information returned by the covariance projection
  optimizer

- `call`: the matched call

## Details

This function is a minor modification of
[`qda`](https://rdrr.io/pkg/MASS/man/qda.html), replacing the classical
sample covariance estimators by projected covariance matrices obtained
using `project_covs()`.

Quadratic discriminant analysis models each class with its own
covariance matrix. In `gipsqda`, these covariance matrices are projected
using the *gips* framework, which enforces permutation symmetry and
mitigates singularity and overfitting in high-dimensional or
small-sample settings.

Classification can be performed using plug-in, predictive, debiased, or
leave-one-out cross-validation rules via `predict.gipsqda`.

## Note

The function may be called with either a formula interface or with a
matrix and grouping factor. Arguments `subset` and `na.action`, if used,
must be named.

## References

Chojecki, A., et al. (2025). *Learning Permutation Symmetry of a
Gaussian Vector with gips in R*. Journal of Statistical Software,
**112**(7), 1–38.
[doi:10.18637/jss.v112.i07](https://doi.org/10.18637/jss.v112.i07)

Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics
with S*. Fourth edition. Springer.

## See also

[`qda`](https://rdrr.io/pkg/MASS/man/qda.html), `predict.gipsqda`,
[`gipslda`](https://AntoniKingston.github.io/gipsDA/reference/gipslda.md),
[`lda`](https://rdrr.io/pkg/MASS/man/lda.html)

## Examples

``` r
tr <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
z <- gipsqda(train, cl)
predict(z,test)$class
#>  [1] s s s s s s s s s s s s s s s s s s s s s s s s s c c c c c c c c v c c c c
#> [39] v c c c c c c c c c c c v v v v v v v v v v v v v v v v v v v v v v v v v
#> Levels: c s v
```
