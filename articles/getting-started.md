# Getting started with gipsDA

## Overview

`gipsDA` extends classical Linear Discriminant Analysis (LDA) and
Quadratic Discriminant Analysis (QDA) by replacing standard empirical
covariance estimates with covariance estimates regularized by
permutation-invariant structures learned with the `gips` framework.

The package provides three model-fitting functions:

- [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
  — LDA with a projected pooled covariance matrix.
- [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
  — QDA with independently projected class-specific covariance matrices.
- [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
  — QDA with class-specific covariance matrices and a jointly estimated
  permutation structure.

The functions follow the familiar conventions of
[`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) and
[`MASS::qda()`](https://rdrr.io/pkg/MASS/man/qda.html): they support
formula, matrix, and data-frame interfaces, and their
[`predict()`](https://rdrr.io/r/stats/predict.html) methods return
predicted classes and posterior probabilities.

## Installation

Install the released version from CRAN:

``` r

install.packages("gipsDA")
```

Or install the development version from GitHub:

``` r

pak::pkg_install("AntoniKingston/gipsDA")
```

Load the package:

``` r

library(gipsDA)
```

## Quick start

We use the built-in `iris` data set and create a simple stratified
train/test split.

``` r

set.seed(42)

train_id <- unlist(
  lapply(split(seq_len(nrow(iris)), iris$Species), sample, size = 35),
  use.names = FALSE
)

train <- iris[train_id, ]
test <- iris[-train_id, ]
```

Fit a
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
model using the formula interface.

``` r

fit <- gipsqda(Species ~ ., data = train)

fit
#> Call:
#> gipsqda(Species ~ ., data = train)
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa         5.020000    3.420000     1.482857   0.2428571
#> versicolor     5.911429    2.771429     4.302857   1.3371429
#> virginica      6.725714    3.020000     5.654286   2.0685714
#> 
#> Permutations with their estimated probabilities:
#> [1] (13)(24)
```

Predict classes for the test set.

``` r

pred <- predict(fit, test)

head(pred$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(pred$posterior)
#>    setosa   versicolor    virginica
#> 6       1 2.909134e-18 1.033310e-27
#> 9       1 1.075425e-13 1.166576e-19
#> 14      1 2.143385e-16 2.068832e-22
#> 16      1 2.508036e-24 3.982454e-36
#> 17      1 4.556129e-22 4.934051e-32
#> 19      1 4.873816e-18 2.042588e-29

mean(pred$class == test$Species)
#> [1] 0.9777778
```

## Preparing data

`gipsDA` works with numeric predictors and a categorical grouping
variable.

Before fitting a model:

- remove or impute missing values,
- remove non-informative columns, such as identifiers,
- encode categorical predictors before passing them to the model,
- avoid mixing variables with incomparable scales when possible.

Because the method searches for symmetries between variables, predictors
should be meaningfully comparable. For example, variables measured in
the same units or representing analogous sensor readings are more
natural candidates for permutation symmetry.

Standardization can be useful when predictors are on very different
scales, but it should be done carefully. In particular, scaling
parameters should be estimated on the training data only and then
applied to new data.

## Model relationships

The three `gipsDA` models are related to classical LDA and QDA as shown
below.

![The diagram illustrates the hierarchical relationships between the
models.](figures/models_hierarchy.png)

The diagram illustrates the hierarchical relationships between the
models.

The figure illustrates the inclusion relations between the model
classes.
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
is contained in QDA,
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
is contained in LDA, and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
lies between
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
and
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md).

## Choosing a model

| Function | Covariance-matrix assumption | When this may be reasonable |
|----|----|----|
| [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) | all classes share one projected covariance matrix | classes differ mainly in their means, but have similar covariance structure |
| [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md) | each class has its own projected covariance matrix and its own permutation structure | classes may have different covariance patterns |
| [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) | each class has its own covariance matrix, but all classes share one permutation structure | classes differ in scale or variance, but have a common dependency pattern between features |

## Key arguments

| Argument | Meaning |
|----|----|
| `prior` | prior probabilities of classes |
| `MAP` | whether to use the Maximum A Posteriori permutation or posterior-weighted averaging |
| `optimizer` | permutation-search strategy: `"BF"` or `"MH"` |
| `max_iter` | number of Metropolis-Hastings iterations when `optimizer = "MH"` |
| `weighted_avg` | covariance pooling strategy for [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) |
| `tol` | numerical tolerance used during fitting |

## MAP: Maximum A Posteriori

`MAP` stands for **Maximum A Posteriori**.

When `MAP = TRUE`, the model selects the single most probable
permutation structure and uses the covariance matrix projected onto that
structure.

``` r

lda_map <- gipslda(
  Species ~ .,
  data = train,
  MAP = TRUE
)

lda_map
#> Call:
#> gipslda(Species ~ ., data = train, MAP = TRUE)
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa         5.020000    3.420000     1.482857   0.2428571
#> versicolor     5.911429    2.771429     4.302857   1.3371429
#> virginica      6.725714    3.020000     5.654286   2.0685714
#> 
#> Coefficients of linear discriminants:
#>                     LD1        LD2
#> Sepal.Length -0.2522284 -0.3624600
#> Sepal.Width   2.3832685 -2.2995034
#> Petal.Length -1.3623096  0.8891933
#> Petal.Width  -3.5000941 -2.3792848
#> 
#> Proportion of trace:
#>    LD1    LD2 
#> 0.9888 0.0112 
#> 
#> Permutations with their estimated probabilities:
#> [1] (1243)
```

When `MAP = FALSE`, the model uses posterior probabilities over retained
permutation structures and computes a posterior-weighted covariance
estimate.

``` r

lda_avg <- gipslda(
  Species ~ .,
  data = train,
  MAP = FALSE
)

lda_avg
#> Call:
#> gipslda(Species ~ ., data = train, MAP = FALSE)
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa         5.020000    3.420000     1.482857   0.2428571
#> versicolor     5.911429    2.771429     4.302857   1.3371429
#> virginica      6.725714    3.020000     5.654286   2.0685714
#> 
#> Coefficients of linear discriminants:
#>                      LD1        LD2
#> Sepal.Length  0.06678184 -0.5324123
#> Sepal.Width   2.14426388 -2.2137470
#> Petal.Length -1.58624789  0.9723419
#> Petal.Width  -3.35358888 -2.3945239
#> 
#> Proportion of trace:
#>    LD1    LD2 
#> 0.9886 0.0114 
#> 
#> Permutations with their estimated probabilities:
#>   (1,2,4,3)  (1,3)(2,4)   (1,2,3,4)  (1,2)(3,4)  (1,4)(2,3) 
#> 0.549732245 0.423533797 0.018304175 0.004073983 0.003658725
```

In printed model output, a result such as `(1243)` is written in cycle
notation. It means that the permutation maps feature `1` to `2`, `2` to
`4`, `4` to `3`, and `3` back to `1`. Features appearing in the same
cycle are treated as exchangeable under the selected symmetry structure.

## Optimizer

The `optimizer` argument controls how permutation structures are
searched.

| Value  | Meaning                    | Recommended use                     |
|--------|----------------------------|-------------------------------------|
| `"BF"` | brute-force search         | small problems, typically `p <= 10` |
| `"MH"` | Metropolis-Hastings search | larger problems, typically `p > 10` |

For `optimizer = "MH"`, use `max_iter` to control the number of
iterations.

``` r

fit_mh <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

For `optimizer = "BF"`, `max_iter` is ignored.

## `weighted_avg` in `gipslda()`

`weighted_avg` controls how the pooled covariance matrix is constructed
before the `gips` projection is applied.

Let `K` be the number of classes, `n_k` the number of observations in
class `k`, `n` the total number of observations, and `S_k` the sample
covariance matrix in class `k`.

With `weighted_avg = FALSE`,
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
uses the classic pooled covariance estimator:

``` math
S_{\mathrm{classic}}
=
\frac{1}{n - K}
\sum_{k = 1}^{K}
(n_k - 1) S_k.
```

With `weighted_avg = TRUE`, it uses the class-size-weighted average:

``` math
S_{\mathrm{weighted}}
=
\frac{1}{n}
\sum_{k = 1}^{K}
n_k S_k.
```

``` r

lda_classic <- gipslda(
  Species ~ .,
  data = train,
  weighted_avg = FALSE
)

lda_weighted <- gipslda(
  Species ~ .,
  data = train,
  weighted_avg = TRUE
)

c(
  classic = mean(predict(lda_classic, test)$class == test$Species),
  weighted = mean(predict(lda_weighted, test)$class == test$Species)
)
#>   classic  weighted 
#> 0.9777778 0.9333333
```

## Fitting all three models

``` r

lda_fit <- gipslda(Species ~ ., data = train)
qda_fit <- gipsqda(Species ~ ., data = train)
joint_qda_fit <- gipsmultqda(Species ~ ., data = train)
```

``` r

lda_pred <- predict(lda_fit, test)
qda_pred <- predict(qda_fit, test)
joint_qda_pred <- predict(joint_qda_fit, test)

c(
  gipslda = mean(lda_pred$class == test$Species),
  gipsqda = mean(qda_pred$class == test$Species),
  gipsmultqda = mean(joint_qda_pred$class == test$Species)
)
#>     gipslda     gipsqda gipsmultqda 
#>   0.9777778   0.9777778   1.0000000
```

## Inspecting fitted models

The most readable way to inspect a fitted model is to print it.

``` r

print(lda_fit)
#> Call:
#> gipslda(Species ~ ., data = train)
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa         5.020000    3.420000     1.482857   0.2428571
#> versicolor     5.911429    2.771429     4.302857   1.3371429
#> virginica      6.725714    3.020000     5.654286   2.0685714
#> 
#> Coefficients of linear discriminants:
#>                     LD1        LD2
#> Sepal.Length -0.2522284 -0.3624600
#> Sepal.Width   2.3832685 -2.2995034
#> Petal.Length -1.3623096  0.8891933
#> Petal.Width  -3.5000941 -2.3792848
#> 
#> Proportion of trace:
#>    LD1    LD2 
#> 0.9888 0.0112 
#> 
#> Permutations with their estimated probabilities:
#> [1] (1243)
```

``` r

print(qda_fit)
#> Call:
#> gipsqda(Species ~ ., data = train)
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa         5.020000    3.420000     1.482857   0.2428571
#> versicolor     5.911429    2.771429     4.302857   1.3371429
#> virginica      6.725714    3.020000     5.654286   2.0685714
#> 
#> Permutations with their estimated probabilities:
#> [1] (13)(24)
```

``` r

print(joint_qda_fit)
#> Call:
#> gipsmultqda(Species ~ ., data = train)
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa         5.020000    3.420000     1.482857   0.2428571
#> versicolor     5.911429    2.771429     4.302857   1.3371429
#> virginica      6.725714    3.020000     5.654286   2.0685714
#> 
#> Permutations with their estimated probabilities:
#> [1] (13)
```

The printed output shows the model call, class priors, group means, and
information about selected or averaged permutation structures.

For a structural summary of the object, use
[`summary()`](https://rdrr.io/r/base/summary.html).

``` r

summary(lda_fit)
#>                   Length Class     Mode     
#> prior              3     -none-    numeric  
#> counts             3     -none-    numeric  
#> means             12     -none-    numeric  
#> scaling            8     -none-    numeric  
#> lev                3     -none-    character
#> svd                2     -none-    numeric  
#> N                  1     -none-    numeric  
#> optimization_info  1     gips_perm list     
#> call               3     -none-    call     
#> terms              3     terms     call     
#> xlevels            0     -none-    list
summary(qda_fit)
#>                   Length Class     Mode     
#> prior              3     -none-    numeric  
#> counts             3     -none-    numeric  
#> means             12     -none-    numeric  
#> scaling           48     -none-    numeric  
#> ldet               3     -none-    numeric  
#> lev                3     -none-    character
#> N                  1     -none-    numeric  
#> call               3     -none-    call     
#> optimization_info  2     gips_perm list     
#> terms              3     terms     call     
#> xlevels            0     -none-    list
summary(joint_qda_fit)
#>                   Length Class     Mode     
#> prior              3     -none-    numeric  
#> counts             3     -none-    numeric  
#> means             12     -none-    numeric  
#> scaling           48     -none-    numeric  
#> ldet               3     -none-    numeric  
#> lev                3     -none-    character
#> N                  1     -none-    numeric  
#> call               3     -none-    call     
#> optimization_info  3     gips_perm list     
#> terms              3     terms     call     
#> xlevels            0     -none-    list
```

Fitted models are list-like S3 objects, so their components can also be
inspected directly.

``` r

names(lda_fit)
#>  [1] "prior"             "counts"            "means"            
#>  [4] "scaling"           "lev"               "svd"              
#>  [7] "N"                 "optimization_info" "call"             
#> [10] "terms"             "xlevels"
names(qda_fit)
#>  [1] "prior"             "counts"            "means"            
#>  [4] "scaling"           "ldet"              "lev"              
#>  [7] "N"                 "call"              "optimization_info"
#> [10] "terms"             "xlevels"
names(joint_qda_fit)
#>  [1] "prior"             "counts"            "means"            
#>  [4] "scaling"           "ldet"              "lev"              
#>  [7] "N"                 "call"              "optimization_info"
#> [10] "terms"             "xlevels"
```

The `optimization_info` component stores information returned by the
permutation optimization step.

``` r

lda_fit$optimization_info
#> [1] (1243)
qda_fit$optimization_info
#> [1] (13)(24)
joint_qda_fit$optimization_info
#> [1] (13)
```

## Prediction output

Prediction returns a list. The most important components are:

| Component   | Meaning                                  |
|-------------|------------------------------------------|
| `class`     | predicted class labels                   |
| `posterior` | posterior class probabilities            |
| `x`         | discriminant coordinates, when available |

``` r

names(lda_pred)
#> [1] "class"     "posterior" "x"

head(lda_pred$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(lda_pred$posterior)
#>    setosa   versicolor    virginica
#> 6       1 7.906717e-22 4.776377e-41
#> 9       1 2.674459e-16 1.730922e-35
#> 14      1 1.353910e-20 1.978581e-41
#> 16      1 3.431995e-28 5.926915e-49
#> 17      1 3.146059e-24 3.515450e-44
#> 19      1 8.270729e-22 2.610171e-41
```

## Matrix interface

The formula interface is usually the most convenient, but the matrix
interface is also available.

``` r

x <- as.matrix(iris[, 1:4])
grouping <- iris$Species

fit_matrix <- gipslda(x, grouping)

predict(fit_matrix, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

The same style can be used with
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md).

``` r

qda_matrix <- gipsqda(x, grouping)
joint_qda_matrix <- gipsmultqda(x, grouping)

predict(qda_matrix, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
predict(joint_qda_matrix, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

## Summary

For a first analysis, start with:

``` r

fit <- gipslda(Species ~ ., data = train)
pred <- predict(fit, test)

mean(pred$class == test$Species)
#> [1] 0.9777778
```

Use:

| Function | Use when |
|----|----|
| [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) | classes can reasonably share one covariance structure |
| [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md) | each class may have its own covariance structure |
| [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) | classes may have different covariance matrices but a shared dependency pattern |

The main tuning choices are:

| Argument | Practical meaning |
|----|----|
| `MAP = TRUE` | use one selected Maximum A Posteriori permutation |
| `MAP = FALSE` | average over retained permutations using posterior probabilities |
| `optimizer = "BF"` | exhaustive search, typically for `p <= 10` |
| `optimizer = "MH"` | stochastic search, typically for `p > 10` |
| `max_iter` | used only with `optimizer = "MH"` |
| `weighted_avg` | changes the pooled covariance estimator in [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) |
