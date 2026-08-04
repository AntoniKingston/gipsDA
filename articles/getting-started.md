# Getting started with gipsDA

## Overview

`gipsDA` extends classical Linear Discriminant Analysis (LDA) and
Quadratic Discriminant Analysis (QDA) by replacing standard empirical
covariance estimates with covariance estimates regularized by
permutation-invariant structures learned with the `gips` framework.

The package is especially useful when the covariance matrix is difficult
to estimate reliably, for example when the number of features is large
relative to the number of observations. In such settings, classical
covariance estimators may be unstable or singular.

The package provides three main model families:

- [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) -
  LDA with a projected pooled covariance matrix.
- [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md) -
  QDA with separately projected class-specific covariance matrices.
- [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) -
  QDA with class-specific covariance matrices and a jointly estimated
  permutation structure.

All three functions follow the familiar conventions of
[`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) and
[`MASS::qda()`](https://rdrr.io/pkg/MASS/man/qda.html). They support
formula, matrix, and data-frame interfaces, and their
[`predict()`](https://rdrr.io/r/stats/predict.html) methods return class
predictions and posterior probabilities.

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

## Example data

We use the built-in `iris` data set. To make the example closer to a
standard classification workflow, we create a stratified train/test
split.

``` r

set.seed(42)

train_id <- unlist(
  lapply(split(seq_len(nrow(iris)), iris$Species), sample, size = 35),
  use.names = FALSE
)

train <- iris[train_id, ]
test <- iris[-train_id, ]

table(train$Species)
#> 
#>     setosa versicolor  virginica 
#>         35         35         35
table(test$Species)
#> 
#>     setosa versicolor  virginica 
#>         15         15         15
```

## Main functions

The main model-fitting functions are:

| Function | Classical analogue | Covariance structure |
|----|----|----|
| [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) | [`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) | one projected pooled covariance matrix |
| [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md) | [`MASS::qda()`](https://rdrr.io/pkg/MASS/man/qda.html) | one projected covariance matrix per class |
| [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) | [`MASS::qda()`](https://rdrr.io/pkg/MASS/man/qda.html) | class-specific covariance matrices with one jointly estimated permutation structure |

The package also provides several S3 methods.

| Method | Objects | Purpose |
|----|----|----|
| [`predict()`](https://rdrr.io/r/stats/predict.html) | `gipslda`, `gipsqda`, `gipsmultqda` | class prediction and posterior probabilities |
| [`print()`](https://rdrr.io/r/base/print.html) | `gipslda`, `gipsqda`, `gipsmultqda` | compact model summary |
| [`coef()`](https://rdrr.io/r/stats/coef.html) | `gipslda` | LDA coefficients |
| [`plot()`](https://rdrr.io/r/graphics/plot.default.html) | `gipslda` | LDA diagnostic plot |
| [`pairs()`](https://rdrr.io/r/graphics/pairs.html) | `gipslda` | pairwise LDA visualization |
| [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) | `gipslda`, `gipsqda`, `gipsmultqda` | extract the model frame when available |

## Model-fitting arguments

The fitting functions follow the standard R modeling style. Depending on
the interface, the model can be fitted using a formula, a matrix, or a
data frame.

### Formula interface arguments

| Argument | Meaning |
|----|----|
| `formula` | model formula, for example `Species ~ .` |
| `data` | data frame containing predictors and the grouping variable |
| `subset` | optional subset of rows used for fitting |
| `na.action` | handling of missing values |
| `...` | additional arguments passed to the corresponding lower-level method |

### Matrix and data-frame interface arguments

| Argument   | Meaning                                           |
|------------|---------------------------------------------------|
| `x`        | matrix or data frame of predictors                |
| `grouping` | class labels                                      |
| `...`      | additional arguments passed to the default method |

### Main statistical and optimization arguments

| Argument | Used in | Meaning |
|----|----|----|
| `prior` | [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md), [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md), [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) | prior probabilities of classes |
| `tol` | mainly [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) | numerical tolerance used during fitting |
| `MAP` | all models | whether to use the single MAP permutation or posterior-weighted averaging |
| `optimizer` | all models | permutation-search strategy, usually `"BF"` or `"MH"` |
| `max_iter` | all models | number of iterations for the Metropolis-Hastings optimizer |
| `weighted_avg` | [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) only | covariance pooling strategy before projection |

The most important `gipsDA`-specific arguments are `MAP`, `optimizer`,
`max_iter`, and, for
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
`weighted_avg`.

## `MAP`: single best permutation vs posterior averaging

The `MAP` argument controls how the permutation-invariant covariance
estimate is constructed.

When `MAP = TRUE`, the model finds the most probable permutation
structure and projects the covariance matrix onto the space invariant
under that structure.

``` r

lda_map <- gipslda(
  Species ~ .,
  data = train,
  MAP = TRUE,
  optimizer = "BF"
)

lda_map
#> Call:
#> gipslda(Species ~ ., data = train, MAP = TRUE, optimizer = "BF")
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
permutations and computes a posterior-weighted covariance projection.

``` r

lda_avg <- gipslda(
  Species ~ .,
  data = train,
  MAP = FALSE,
  optimizer = "BF"
)

lda_avg
#> Call:
#> gipslda(Species ~ ., data = train, MAP = FALSE, optimizer = "BF")
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

Use `MAP = TRUE` when you want a single selected structure. Use
`MAP = FALSE` when you want to account for uncertainty in the selected
permutation structure.

## `optimizer`: brute force and Metropolis-Hastings

The `optimizer` argument controls how the permutation structure is
searched.

| Value  | Meaning                    | Recommended use           |
|--------|----------------------------|---------------------------|
| `"BF"` | brute-force search         | small number of features  |
| `"MH"` | Metropolis-Hastings search | larger number of features |

Brute force is deterministic and searches the relevant permutation space
exhaustively, but its cost grows quickly with the number of features.

``` r

lda_bf <- gipslda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)

lda_bf
#> Call:
#> gipslda(Species ~ ., data = train, optimizer = "BF")
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

For larger feature spaces, use `"MH"` and control runtime with
`max_iter`.

``` r

lda_mh <- gipslda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

## `max_iter`

`max_iter` controls the number of Metropolis-Hastings iterations when
`optimizer = "MH"`.

``` r

fit <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

Increasing `max_iter` may improve exploration of the permutation space,
but also increases runtime.

For `optimizer = "BF"`, `max_iter` is not the main computational
control.

## `prior`

The `prior` argument specifies class prior probabilities.

By default, class priors are estimated from the training data. You can
override them manually.

``` r

equal_prior <- rep(1 / length(levels(train$Species)), length(levels(train$Species)))
names(equal_prior) <- levels(train$Species)

lda_prior <- gipslda(
  Species ~ .,
  data = train,
  prior = equal_prior,
  optimizer = "BF"
)

lda_prior$prior
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333
```

The prior vector should have one value per class and should sum to one.

## `weighted_avg` in `gipslda()`

`weighted_avg` is specific to
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
and controls how the pooled covariance matrix is constructed before
applying the `gips` projection.

``` r

lda_classic <- gipslda(
  Species ~ .,
  data = train,
  weighted_avg = FALSE,
  optimizer = "BF"
)

lda_weighted <- gipslda(
  Species ~ .,
  data = train,
  weighted_avg = TRUE,
  optimizer = "BF"
)
```

Use:

- `weighted_avg = FALSE` for the classic LDA-style pooled covariance
  estimator.
- `weighted_avg = TRUE` for the class-size-weighted averaging variant.

## Linear discriminant analysis with `gipslda()`

The formula interface is the most convenient way to fit a model.

``` r

lda_fit <- gipslda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)

lda_fit
#> Call:
#> gipslda(Species ~ ., data = train, optimizer = "BF")
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

Prediction returns class labels and posterior probabilities.

``` r

lda_pred <- predict(lda_fit, test)

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

mean(lda_pred$class == test$Species)
#> [1] 0.9777778
```

## Prediction methods for `gipslda()`

For
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
objects, prediction can be performed with different methods.

``` r

pred_plugin <- predict(lda_fit, test, method = "plug-in")
pred_predictive <- predict(lda_fit, test, method = "predictive")
pred_debiased <- predict(lda_fit, test, method = "debiased")

head(pred_plugin$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(pred_predictive$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(pred_debiased$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

The available prediction methods are:

| Method         | Meaning                            |
|----------------|------------------------------------|
| `"plug-in"`    | standard plug-in discriminant rule |
| `"predictive"` | predictive rule                    |
| `"debiased"`   | debiased rule                      |

## Quadratic discriminant analysis with `gipsqda()`

[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
estimates a separate projected covariance matrix for each class.

``` r

qda_fit <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)

qda_fit
#> Call:
#> gipsqda(Species ~ ., data = train, optimizer = "BF")
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

Prediction works analogously to classical QDA.

``` r

qda_pred <- predict(qda_fit, test)

head(qda_pred$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(qda_pred$posterior)
#>    setosa   versicolor    virginica
#> 6       1 2.909134e-18 1.033310e-27
#> 9       1 1.075425e-13 1.166576e-19
#> 14      1 2.143385e-16 2.068832e-22
#> 16      1 2.508036e-24 3.982454e-36
#> 17      1 4.556129e-22 4.934051e-32
#> 19      1 4.873816e-18 2.042588e-29

mean(qda_pred$class == test$Species)
#> [1] 0.9777778
```

## Joint quadratic discriminant analysis with `gipsmultqda()`

[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
is an intermediate model between
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
and
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md).

Like QDA, it allows each class to have its own covariance matrix. Unlike
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md),
it estimates one shared permutation structure for all class covariance
matrices. This is useful when classes may differ in scale or variance,
but are expected to share a common dependency pattern between features.

``` r

joint_qda_fit <- gipsmultqda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)

joint_qda_fit
#> Call:
#> gipsmultqda(Species ~ ., data = train, optimizer = "BF")
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

Prediction returns class labels and posterior probabilities.

``` r

joint_qda_pred <- predict(joint_qda_fit, test)

head(joint_qda_pred$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(joint_qda_pred$posterior)
#>    setosa   versicolor    virginica
#> 6       1 1.809368e-18 3.630155e-32
#> 9       1 4.185357e-13 9.716563e-23
#> 14      1 8.534255e-16 4.000877e-26
#> 16      1 2.045692e-24 3.136984e-41
#> 17      1 1.110504e-21 8.419286e-37
#> 19      1 2.676609e-18 3.735764e-34

mean(joint_qda_pred$class == test$Species)
#> [1] 1
```

## Comparing the three models

``` r

acc <- c(
  gipslda = mean(lda_pred$class == test$Species),
  gipsqda = mean(qda_pred$class == test$Species),
  gipsmultqda = mean(joint_qda_pred$class == test$Species)
)

acc
#>     gipslda     gipsqda gipsmultqda 
#>   0.9777778   0.9777778   1.0000000
```

## Leave-one-out prediction for QDA models

For
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
objects, leave-one-out cross-validation can be requested by omitting
`newdata` and using `method = "looCV"`.

In this mode, each training observation is classified as if it had been
left out when fitting the model. This gives an internal estimate of
classification accuracy without creating a separate test set.

``` r

qda_loo <- predict(qda_fit, method = "looCV")
joint_qda_loo <- predict(joint_qda_fit, method = "looCV")

c(
  gipsqda_loo_accuracy = mean(qda_loo$class == train$Species),
  gipsmultqda_loo_accuracy = mean(joint_qda_loo$class == train$Species)
)
#>     gipsqda_loo_accuracy gipsmultqda_loo_accuracy 
#>                0.9714286                0.9714286
```

## Formula interface

The formula interface is usually the easiest one to use.

``` r

fit_formula <- gipslda(
  Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
  data = train,
  optimizer = "BF"
)

predict(fit_formula, test)$class[1:10]
#>  [1] setosa setosa setosa setosa setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

Formula methods also support standard arguments such as `subset` and
`na.action` through the usual R model-frame workflow.

``` r

fit_subset <- gipslda(
  Species ~ .,
  data = iris,
  subset = Species != "setosa",
  optimizer = "BF"
)
#> Warning in gipslda.default(x, grouping, ...): group setosa is empty

fit_subset
#> Call:
#> gipslda(Species ~ ., data = iris, optimizer = "BF", subset = Species != 
#>     "setosa")
#> 
#> Prior probabilities of groups:
#> versicolor  virginica 
#>        0.5        0.5 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> versicolor        5.936       2.770        4.260       1.326
#> virginica         6.588       2.974        5.552       2.026
#> 
#> Coefficients of linear discriminants:
#>                    LD1
#> Sepal.Length -1.094122
#> Sepal.Width  -1.229136
#> Petal.Length  1.864973
#> Petal.Width   3.324524
#> 
#> Permutations with their estimated probabilities:
#> [1] (13)(24)
```

## Matrix interface

The matrix interface separates predictors from class labels.

``` r

x <- as.matrix(iris[, 1:4])
grouping <- iris$Species

fit_matrix <- gipslda(
  x,
  grouping,
  optimizer = "BF"
)

predict(fit_matrix, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

The same interface can be used with
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md).

``` r

qda_matrix <- gipsqda(
  x,
  grouping,
  optimizer = "BF"
)

predict(qda_matrix, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

And with
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md).

``` r

joint_qda_matrix <- gipsmultqda(
  x,
  grouping,
  optimizer = "BF"
)

predict(joint_qda_matrix, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

## Data-frame interface

You can also pass predictors as a data frame and labels separately.

``` r

x_df <- iris[, 1:4]
y <- iris$Species

fit_df <- gipslda(
  x_df,
  y,
  optimizer = "BF"
)

predict(fit_df, x_df[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

The data-frame interface is also available for QDA-type models.

``` r

qda_df <- gipsqda(
  x_df,
  y,
  optimizer = "BF"
)

joint_qda_df <- gipsmultqda(
  x_df,
  y,
  optimizer = "BF"
)

predict(qda_df, x_df[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
predict(joint_qda_df, x_df[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

## Fitted objects

Fitted objects contain the usual discriminant-analysis components, such
as class priors, class counts, group means, and scaling information.

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

`gipsDA` objects additionally store optimization information.

``` r

lda_fit$optimization_info
#> [1] (1243)
qda_fit$optimization_info
#> [1] (13)(24)
joint_qda_fit$optimization_info
#> [1] (13)
```

When `MAP = TRUE`, this information corresponds to the selected
permutation structure. When `MAP = FALSE`, it contains
posterior-probability information used for averaging.

## Inspecting fitted model objects

Fitted `gipsDA` models are standard R objects. They contain the usual
discriminant-analysis components, such as class priors, class counts,
group means, scaling information, and additional information about the
permutation optimization.

The most readable way to inspect a fitted model is to print it.

``` r

print(lda_fit)
#> Call:
#> gipslda(Species ~ ., data = train, optimizer = "BF")
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
#> gipsqda(Species ~ ., data = train, optimizer = "BF")
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
#> gipsmultqda(Species ~ ., data = train, optimizer = "BF")
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

The printed output shows the main fitted quantities, including prior
probabilities, group means, and information about the permutation
structures selected or averaged by the model.

For a more structural view, use
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
#> call               4     -none-    call     
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
#> call               4     -none-    call     
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
#> call               4     -none-    call     
#> optimization_info  3     gips_perm list     
#> terms              3     terms     call     
#> xlevels            0     -none-    list
```

Because fitted models are list-like S3 objects, it is also useful to
inspect their component names directly.

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

The most important components are usually:

| Component | Meaning |
|----|----|
| `prior` | prior probabilities of classes |
| `counts` | number of observations in each class |
| `means` | class-wise feature means |
| `scaling` | scaling or decomposition information used for prediction |
| `ldet` | log-determinant information for QDA-type models |
| `lev` | class labels |
| `call` | original function call |
| `optimization_info` | information returned by the permutation optimization step |

The exact set of components can differ between LDA-type and QDA-type
models.

A compact way to inspect the internal structure is to create a small
helper.

``` r

inspect_model <- function(object) {
  data.frame(
    component = names(object),
    class = vapply(
      object,
      function(x) paste(class(x), collapse = ", "),
      character(1)
    ),
    length = vapply(object, length, integer(1)),
    dim = vapply(
      object,
      function(x) {
        d <- dim(x)
        if (is.null(d)) "" else paste(d, collapse = " x ")
      },
      character(1)
    ),
    row.names = NULL
  )
}
```

``` r

inspect_model(lda_fit)
#>            component          class length   dim
#> 1              prior        numeric      3      
#> 2             counts        integer      3      
#> 3              means  matrix, array     12 3 x 4
#> 4            scaling  matrix, array      8 4 x 2
#> 5                lev      character      3      
#> 6                svd        numeric      2      
#> 7                  N        integer      1      
#> 8  optimization_info      gips_perm      1      
#> 9               call           call      4      
#> 10             terms terms, formula      3      
#> 11           xlevels           list      0
```

``` r

inspect_model(qda_fit)
#>            component          class length       dim
#> 1              prior        numeric      3          
#> 2             counts        integer      3          
#> 3              means  matrix, array     12     3 x 4
#> 4            scaling          array     48 4 x 4 x 3
#> 5               ldet        numeric      3          
#> 6                lev      character      3          
#> 7                  N        integer      1          
#> 8               call           call      4          
#> 9  optimization_info      gips_perm      2          
#> 10             terms terms, formula      3          
#> 11           xlevels           list      0
```

``` r

inspect_model(joint_qda_fit)
#>            component          class length       dim
#> 1              prior        numeric      3          
#> 2             counts        integer      3          
#> 3              means  matrix, array     12     3 x 4
#> 4            scaling          array     48 4 x 4 x 3
#> 5               ldet        numeric      3          
#> 6                lev      character      3          
#> 7                  N        integer      1          
#> 8               call           call      4          
#> 9  optimization_info      gips_perm      3          
#> 10             terms terms, formula      3          
#> 11           xlevels           list      0
```

The `optimization_info` component stores the result of the `gips`
optimization.

``` r

lda_fit$optimization_info
#> [1] (1243)
```

``` r

qda_fit$optimization_info
#> [1] (13)(24)
```

``` r

joint_qda_fit$optimization_info
#> [1] (13)
```

When `MAP = TRUE`, this component describes the selected permutation
structure. When `MAP = FALSE`, it contains posterior-probability
information used for weighted averaging.

## LDA diagnostics and visualization

For
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
objects, standard diagnostic and visualization methods are available.

``` r

coef(lda_fit)
#>                     LD1        LD2
#> Sepal.Length -0.2522284 -0.3624600
#> Sepal.Width   2.3832685 -2.2995034
#> Petal.Length -1.3623096  0.8891933
#> Petal.Width  -3.5000941 -2.3792848
```

``` r

plot(lda_fit)
```

![](getting-started_files/figure-html/unnamed-chunk-44-1.png)

``` r

pairs(lda_fit, type = "std")
```

![](getting-started_files/figure-html/unnamed-chunk-45-1.png)

Depending on the object and method, additional plot types may be
available.

``` r

pairs(lda_fit, type = "trellis")
```

## Prediction output

Prediction methods return objects containing predicted classes and
posterior probabilities.

``` r

names(lda_pred)
#> [1] "class"     "posterior" "x"
names(qda_pred)
#> [1] "class"     "posterior"
names(joint_qda_pred)
#> [1] "class"     "posterior"
```

The most commonly used components are:

| Component   | Meaning                                  |
|-------------|------------------------------------------|
| `class`     | predicted class labels                   |
| `posterior` | posterior class probabilities            |
| `x`         | discriminant coordinates, when available |

For example:

``` r

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

## Practical recommendations

Start with
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
when a common covariance structure across classes is a reasonable
assumption.

``` r

fit <- gipslda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)
```

Use
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
when each class may have a different covariance structure.

``` r

fit <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)
```

Use
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
when class covariance matrices may differ but are expected to share a
common permutation structure.

``` r

fit <- gipsmultqda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)
```

For small numbers of features, use `optimizer = "BF"`.

``` r

fit <- gipslda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)
```

For larger numbers of features, use `optimizer = "MH"` and set
`max_iter`.

``` r

fit <- gipslda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

Use `MAP = FALSE` when you want to account for uncertainty in the
selected permutation structure.

``` r

fit <- gipslda(
  Species ~ .,
  data = train,
  MAP = FALSE,
  optimizer = "BF"
)
```

Use `weighted_avg = TRUE` only for
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
when you want the alternative class-size-weighted covariance pooling
strategy.

``` r

fit <- gipslda(
  Species ~ .,
  data = train,
  weighted_avg = TRUE,
  optimizer = "BF"
)
```

## Summary

`gipsDA` provides three main drop-in alternatives to classical
discriminant analysis:

- [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
  regularizes LDA by projecting the pooled covariance matrix onto a
  permutation-invariant structure.
- [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
  regularizes QDA by projecting each class-specific covariance matrix
  independently.
- [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
  regularizes QDA by estimating one shared permutation structure for all
  class-specific covariance matrices.

The most important arguments are:

| Argument       | Main use                                              |
|----------------|-------------------------------------------------------|
| `MAP`          | choose MAP projection or posterior-weighted averaging |
| `optimizer`    | choose `"BF"` or `"MH"`                               |
| `max_iter`     | control MH runtime                                    |
| `prior`        | set class priors                                      |
| `weighted_avg` | choose the LDA covariance pooling strategy            |
| `tol`          | numerical tolerance                                   |

For new users, a good starting point is:

``` r

fit <- gipslda(
  Species ~ .,
  data = train,
  MAP = TRUE,
  optimizer = "BF"
)

pred <- predict(fit, test)

mean(pred$class == test$Species)
#> [1] 0.9777778
```
