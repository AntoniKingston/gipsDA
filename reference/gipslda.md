# Linear Discriminant Analysis with gips Covariance Projection

Linear discriminant analysis (LDA) using covariance matrices projected
via the *gips* framework to enforce permutation symmetry and improve
numerical stability.

## Usage

``` r
gipslda(x, ...)

# S3 method for class 'formula'
gipslda(formula, data, ..., subset, na.action)

# Default S3 method
gipslda(x, grouping, prior = proportions,
  tol = 1e-4, weighted_avg = FALSE,
  MAP = TRUE, optimizer = NULL, max_iter = NULL, ...)

# S3 method for class 'data.frame'
gipslda(x, ...)

# S3 method for class 'matrix'
gipslda(x, grouping, ..., subset, na.action)
```

## Arguments

- x:

  (required if no formula is given as the principal argument) a matrix
  or data frame or Matrix containing the explanatory variables.

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

  (required if no formula principal argument is given) a factor
  specifying the class for each observation.

- prior:

  The prior probabilities of class membership. If unspecified, the class
  proportions for the training set are used.

- tol:

  A tolerance to decide if a matrix is singular; variables whose
  variance is less than `tol^2` are rejected.

- subset:

  An index vector specifying the cases to be used in the training
  sample. (NOTE: must be named.)

- na.action:

  A function specifying the action for `NA`s. \#' @param weighted_avg
  Logical; if TRUE, uses a weighted average of class-specific covariance
  matrices instead of the pooled covariance.

- MAP:

  Logical; whether to compute a Maximum A Posteriori gips projection of
  the covariance matrix.

- optimizer:

  Character; optimization method used by gips (e.g. `"BF"` or `"MH"`).

- max_iter:

  Maximum number of iterations for the optimizer.

- weighted_avg:

  Logical; Whether to compute scatter from all classes at once or to
  compute them within classes and compute the main one as average
  weighted by class proportions.

## Value

An object of class `"gipslda"` containing:

- `prior`: prior class probabilities

- `counts`: number of observations per class

- `means`: group means

- `scaling`: linear discriminant coefficients

- `svd`: singular values of the between-class scatter

- `N`: number of observations

- `optimization_info`: information about the gips optimization

- `call`: matched call

## Details

This function is a minor modification of
[`lda`](https://rdrr.io/pkg/MASS/man/lda.html), replacing the classical
sample covariance estimators by projected covariance matrices obtained
using `project_covs()`.

Unlike classical LDA, the within-class covariance matrix is first
projected onto a permutation-invariant structure using the gips
framework. This can stabilize covariance estimation in high dimensions
or when symmetry assumptions are justified.

The choice of optimizer and MAP estimation affects both the covariance
estimate and the resulting discriminant directions.

See Chojecki et al. (2025) for theoretical background.

## Note

This function is inspired by
[`lda`](https://rdrr.io/pkg/MASS/man/lda.html) but is not a drop-in
replacement. The covariance estimator, optimization procedure, and
returned object differ substantially.

## References

Chojecki, A., et al. (2025). *Learning Permutation Symmetry of a
Gaussian Vector with gips in R*. Journal of Statistical Software,
**112**(7), 1–38.
[doi:10.18637/jss.v112.i07](https://doi.org/10.18637/jss.v112.i07)

## See also

[`lda`](https://rdrr.io/pkg/MASS/man/lda.html),
[`gips`](https://przechoj.github.io/gips/reference/gips.html)

## Examples

``` r
Iris <- data.frame(rbind(iris3[, , 1], iris3[, , 2], iris3[, , 3]),
  Sp = rep(c("s", "c", "v"), rep(50, 3))
)
train <- sample(1:150, 75)
z <- gipslda(Sp ~ ., Iris, prior = c(1, 1, 1) / 3, subset = train)
predict(z, Iris[-train, ])$class
#>  [1] s s s s s s s s s s s s s s s s s s s s s s s s s c c c c c c c c c c c c c
#> [39] c c c v c c c c c c v v v v v v v v v v v v v v v c v v v v v v v v v v v
#> Levels: c s v
(z1 <- update(z, . ~ . - Petal.W.))
#> Call:
#> gipslda(Sp ~ Sepal.L. + Sepal.W. + Petal.L., data = Iris, prior = c(1, 
#>     1, 1)/3, subset = train)
#> 
#> Prior probabilities of groups:
#>         c         s         v 
#> 0.3333333 0.3333333 0.3333333 
#> 
#> Group means:
#>   Sepal.L. Sepal.W. Petal.L.
#> c 5.885185 2.766667 4.166667
#> s 5.020000 3.476000 1.452000
#> v 6.630435 2.952174 5.508696
#> 
#> Coefficients of linear discriminants:
#>                 LD1       LD2
#> Sepal.L.  0.5096716 -1.145938
#> Sepal.W.  1.0822205  3.458572
#> Petal.L. -2.5577491  0.969699
#> 
#> Proportion of trace:
#>    LD1    LD2 
#> 0.9884 0.0116 
#> 
#> Permutations with their estimated probabilities:
#> [1] (23)
```
