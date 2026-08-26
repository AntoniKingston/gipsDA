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
  MAP = TRUE, optimizer = NULL, max_iter = NULL,
  show_progress_bar = FALSE, ...)

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

  A function specifying the action for `NA`s.

- MAP:

  Logical; whether to compute a Maximum A Posteriori gips projection. If
  `FALSE`, projected matrices are averaged using retained posterior
  permutation probabilities.

- optimizer:

  Character; optimization method used by gips (e.g. `"BF"` or `"MH"`).

- max_iter:

  Maximum number of iterations for the optimizer.

- show_progress_bar:

  Logical; if `TRUE`, display the progress bar from the underlying gips
  optimizer. Defaults to `FALSE`.

- weighted_avg:

  Logical; if `FALSE`, use the pooled within-class scatter matrix. If
  `TRUE`, use a class-proportion-weighted average of class-specific
  covariance matrices.

## Value

An object of class `"gipslda"` containing:

- `prior`: prior class probabilities

- `counts`: number of observations per class

- `means`: group means

- `scaling`: transformation matrix of linear discriminants

- `lev`: class labels

- `svd`: singular values of the between-class scatter

- `N`: number of observations

- `optimization_info`: estimated probabilities of retained permutations
  returned by the gips optimization

- `call`: matched call

- Formula fits additionally contain `terms`, `contrasts`, `xlevels`, and
  any recorded `na.action`.

## Details

This function is a minor modification of
[`lda`](https://rdrr.io/pkg/MASS/man/lda.html), replacing the classical
sample covariance estimators by projected covariance matrices obtained
using the [`gips`](https://przechoj.github.io/gips/reference/gips.html)
framework.

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
#>  [1] s s s s s s s s s s s s s s s s s s c c c c c c c c c c v c c c v c c c c c
#> [39] c c c c c v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v v
#> Levels: c s v
(z1 <- update(z, . ~ . - Petal.W.))
#> Call:
#> gipslda(Sp ~ Sepal.L. + Sepal.W. + Petal.L., data = Iris, prior = c(1, 
#>     1, 1)/3, subset = train)
#> 
#> Model: gipslda 
#> Number of observations: 75 
#> Number of groups: 3 
#> Number of predictors: 3 
#> 
#> Fitting options:
#> $MAP
#> [1] TRUE
#> 
#> $optimizer
#> [1] "BF"
#> 
#> $max_iter
#> NULL
#> 
#> $weighted_avg
#> [1] FALSE
#> 
#> 
#> Prior probabilities of groups:
#>         c         s         v 
#> 0.3333333 0.3333333 0.3333333 
#> 
#> Class counts:
#>  c  s  v 
#> 25 32 18 
#> 
#> Group means:
#>   Sepal.L. Sepal.W. Petal.L.
#> c 6.032000 2.836000 4.304000
#> s 4.978125 3.390625 1.468750
#> v 6.666667 3.005556 5.566667
#> 
#> Permutations with their estimated probabilities:
#>       (2,3)          ()       (1,3)     (1,2,3) 
#> 0.891542235 0.067187895 0.038375014 0.002848842 
#> 
#> Coefficients of linear discriminants:
#>                 LD1       LD2
#> Sepal.L.  0.7785584 -1.435165
#> Sepal.W.  0.8639359  3.865342
#> Petal.L. -3.7685179  1.017625
#> 
#> Proportion of trace:
#>    LD1    LD2 
#> 0.9949 0.0051 
```
