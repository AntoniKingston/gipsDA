# Advanced usage of gipsDA

## Overview

This vignette gives a more detailed overview of the main modeling
options in `gipsDA`.

It covers:

- data preparation,
- the differences between
  [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
  [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md),
  and
  [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md),
- the `MAP`, `optimizer`, `max_iter`, `prior`, and `weighted_avg`
  arguments,
- interpretation of permutation output,
- prediction methods,
- leave-one-out prediction,
- model inspection and diagnostics.

For a shorter first example, see the “Getting started” vignette.

``` r

library(gipsDA)
```

## Example data

We use the built-in `iris` data set.

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

## Preparing data

`gipsDA` models assume numeric predictors and a categorical grouping
variable.

Before fitting a model, it is usually useful to:

- remove identifier columns,
- remove or impute missing values,
- remove constant or almost constant columns,
- encode categorical predictors,
- check whether predictors are on comparable scales.

The last point is especially important for `gipsDA`. The method searches
for permutation symmetries between variables. Such symmetries are most
meaningful when variables are comparable, for example when they are
measured in the same units or represent analogous sensor readings.

### Checking predictor types

``` r

str(train)
#> 'data.frame':    105 obs. of  5 variables:
#>  $ Sepal.Length: num  5.3 5.5 5.1 4.8 4.9 5 5.1 5.1 4.6 5.1 ...
#>  $ Sepal.Width : num  3.7 3.5 3.5 3.4 3.1 3.2 3.5 3.3 3.4 3.8 ...
#>  $ Petal.Length: num  1.5 1.3 1.4 1.9 1.5 1.2 1.4 1.7 1.4 1.9 ...
#>  $ Petal.Width : num  0.2 0.2 0.2 0.2 0.1 0.2 0.3 0.5 0.3 0.4 ...
#>  $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
```

For formula-based usage, the response variable should be a factor or a
categorical variable.

``` r

is.factor(train$Species)
#> [1] TRUE
```

The predictors in `iris` are already numeric.

``` r

vapply(train[, 1:4], is.numeric, logical(1))
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>         TRUE         TRUE         TRUE         TRUE
```

### Scaling predictors

Scaling may be useful when predictors are measured on very different
scales. However, scaling should be done carefully: the center and scale
should be estimated only on the training data and then applied to new
data.

``` r

x_train <- train[, 1:4]
x_test <- test[, 1:4]

train_center <- vapply(x_train, mean, numeric(1))
train_scale <- vapply(x_train, sd, numeric(1))

train_scaled <- train
test_scaled <- test

train_scaled[, 1:4] <- scale(
  x_train,
  center = train_center,
  scale = train_scale
)

test_scaled[, 1:4] <- scale(
  x_test,
  center = train_center,
  scale = train_scale
)
```

Then fit the model on the scaled data.

``` r

fit_scaled <- gipslda(Species ~ ., data = train_scaled)
pred_scaled <- predict(fit_scaled, test_scaled)

mean(pred_scaled$class == test_scaled$Species)
#> [1] 0.9777778
```

Scaling is not always necessary. It depends on whether the original
variables are already comparable and whether scaling is meaningful for
the application.

## Model choice

The package provides three main classifiers.

| Function | Covariance-matrix assumption | Typical use case |
|----|----|----|
| [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) | all classes share one projected covariance matrix | classes differ mainly in their means |
| [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md) | each class has its own projected covariance matrix and its own permutation structure | classes may have different covariance patterns |
| [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) | each class has its own covariance matrix, but all classes share one permutation structure | classes may differ in scale, but share a dependency pattern |

Fit all three models on the same data.

``` r

lda_fit <- gipslda(Species ~ ., data = train)
qda_fit <- gipsqda(Species ~ ., data = train)
joint_qda_fit <- gipsmultqda(Species ~ ., data = train)
```

Compare test-set accuracy.

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

## Model hierarchy

The inclusion relations between the model classes are shown below.

![The diagram illustrates the hierarchical relationships between the
models.](figures/models_hierarchy.png)

The diagram illustrates the hierarchical relationships between the
models.

The figure shows that
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
is contained in QDA,
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
is contained in LDA, and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
lies between
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
and
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md).

In practical terms, moving inward in the diagram means imposing stronger
assumptions on the covariance structure. Stronger assumptions can reduce
estimation variance, especially when the number of observations is
small, but they may be too restrictive if the assumed structure is not
present in the data.

## MAP: Maximum A Posteriori

`MAP` stands for **Maximum A Posteriori**.

The `MAP` argument controls how the covariance matrix is projected after
the permutation search.

When `MAP = TRUE`, the model selects the single most probable
permutation structure and projects the covariance matrix onto the
invariant space determined by that structure.

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
permutation structures and computes a posterior-weighted projection.

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

Conceptually, if $`S`$ is an empirical covariance matrix and $`S_c`$ is
its projection under permutation structure $`c`$, then the two
approaches can be summarized as follows.

For `MAP = TRUE`:

``` math
\hat{S} = S_{c^*},
```

where

``` math
c^* = \arg\max_c P(c \mid X).
```

For `MAP = FALSE`:

``` math
\hat{S} =
\sum_c P(c \mid X) S_c.
```

The second option averages over several possible symmetry structures
instead of using only one selected structure.

## Interpreting permutation output

Printed model output may contain permutations such as:

``` text
(1,2)
(1,2)(3,4)
(1,2,4,3)
()
```

This notation is called cycle notation.

For example, the cycle

``` text
(1,2,4,3)
```

means that the permutation maps:

``` text
1 -> 2
2 -> 4
4 -> 3
3 -> 1
```

A product of cycles such as:

``` text
(1,2)(3,4)
```

means that feature 1 is swapped with feature 2, and feature 3 is swapped
with feature 4.

The empty permutation:

``` text
()
```

is the identity permutation. It means that no non-trivial permutation
symmetry was selected.

In the context of `gipsDA`, a selected permutation structure describes
which features are treated as exchangeable by the covariance model.

## Optimizer

The `optimizer` argument controls how permutation structures are
searched.

| Value  | Meaning                    | Typical use                         |
|--------|----------------------------|-------------------------------------|
| `"BF"` | brute-force search         | small problems, typically `p <= 10` |
| `"MH"` | Metropolis-Hastings search | larger problems, typically `p > 10` |

The brute-force optimizer searches the relevant permutation space
exhaustively. It is deterministic, but its cost grows quickly with the
number of features.

``` r

fit_bf <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "BF"
)

fit_bf
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

For larger problems, use the Metropolis-Hastings optimizer.

``` r

fit_mh <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

## `max_iter`

`max_iter` controls the number of Metropolis-Hastings iterations when
`optimizer = "MH"`.

Increasing `max_iter` gives the stochastic search more time to explore
the permutation space, but also increases runtime.

``` r

fit_mh_100 <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 100
)

fit_mh_1000 <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

For `optimizer = "BF"`, `max_iter` is ignored.

## Class priors

The `prior` argument specifies prior probabilities of classes.

By default, priors are estimated from the training data.

``` r

lda_fit$prior
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333
```

You can set them manually.

``` r

equal_prior <- rep(1 / length(levels(train$Species)), length(levels(train$Species)))
names(equal_prior) <- levels(train$Species)

lda_equal_prior <- gipslda(
  Species ~ .,
  data = train,
  prior = equal_prior
)

lda_equal_prior$prior
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333
```

The prior vector should contain one value per class and should sum to
one.

``` r

sum(equal_prior)
#> [1] 1
```

Class priors affect posterior probabilities and may affect predicted
classes, especially when classes overlap.

## `weighted_avg` in `gipslda()`

The `weighted_avg` argument is specific to
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md).

It controls how the pooled covariance matrix is constructed before the
`gips` projection is applied.

Let:

- $`K`$ be the number of classes,
- $`n_k`$ be the number of observations in class $`k`$,
- $`n = \sum_{k=1}^{K} n_k`$ be the total number of observations,
- $`S_k`$ be the sample covariance matrix in class $`k`$.

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

With `weighted_avg = TRUE`,
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
uses:

``` math
S_{\mathrm{weighted}}
=
\frac{1}{n}
\sum_{k = 1}^{K}
n_k S_k.
```

Fit both variants.

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
```

Compare predictions.

``` r

pred_classic <- predict(lda_classic, test)
pred_weighted <- predict(lda_weighted, test)

c(
  classic = mean(pred_classic$class == test$Species),
  weighted = mean(pred_weighted$class == test$Species)
)
#>   classic  weighted 
#> 0.9777778 0.9333333
```

The two variants can behave differently when class sizes are imbalanced
or when class-specific covariance estimates differ substantially.

## Formula interface

The formula interface is usually the most convenient interface.

``` r

fit_formula <- gipslda(
  Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
  data = train
)

fit_formula
#> Call:
#> gipslda(Species ~ Sepal.Length + Sepal.Width + Petal.Length + 
#>     Petal.Width, data = train)
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

The shorthand `Species ~ .` uses all remaining columns as predictors.

``` r

fit_formula_short <- gipslda(
  Species ~ .,
  data = train
)

fit_formula_short
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

Formula methods also support `subset`.

``` r

fit_subset <- gipslda(
  Species ~ .,
  data = iris,
  subset = Species != "setosa"
)
#> Warning in gipslda.default(x, grouping, ...): group setosa is empty

fit_subset
#> Call:
#> gipslda(Species ~ ., data = iris, subset = Species != "setosa")
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

fit_matrix_lda <- gipslda(x, grouping)
fit_matrix_qda <- gipsqda(x, grouping)
fit_matrix_joint <- gipsmultqda(x, grouping)
```

Prediction can then be performed on a matrix with the same columns.

``` r

predict(fit_matrix_lda, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
predict(fit_matrix_qda, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
predict(fit_matrix_joint, x[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

## Data-frame interface

Predictors can also be passed as a data frame, with labels supplied
separately.

``` r

x_df <- iris[, 1:4]
y <- iris$Species

fit_df_lda <- gipslda(x_df, y)
fit_df_qda <- gipsqda(x_df, y)
fit_df_joint <- gipsmultqda(x_df, y)
```

``` r

predict(fit_df_lda, x_df[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
predict(fit_df_qda, x_df[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
predict(fit_df_joint, x_df[1:5, ])$class
#> [1] setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
```

## Prediction output

Prediction returns a list.

``` r

pred <- predict(lda_fit, test)

names(pred)
#> [1] "class"     "posterior" "x"
```

The most important components are:

| Component   | Meaning                                  |
|-------------|------------------------------------------|
| `class`     | predicted class labels                   |
| `posterior` | posterior class probabilities            |
| `x`         | discriminant coordinates, when available |

``` r

head(pred$class)
#> [1] setosa setosa setosa setosa setosa setosa
#> Levels: setosa versicolor virginica
head(pred$posterior)
#>    setosa   versicolor    virginica
#> 6       1 7.906717e-22 4.776377e-41
#> 9       1 2.674459e-16 1.730922e-35
#> 14      1 1.353910e-20 1.978581e-41
#> 16      1 3.431995e-28 5.926915e-49
#> 17      1 3.146059e-24 3.515450e-44
#> 19      1 8.270729e-22 2.610171e-41
```

For
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md),
the output has the same main structure.

``` r

names(qda_pred)
#> [1] "class"     "posterior"
names(joint_qda_pred)
#> [1] "class"     "posterior"
```

## Prediction methods for `gipslda()`

For
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
objects, [`predict()`](https://rdrr.io/r/stats/predict.html) supports
the same main prediction-method names as
[`MASS::predict.lda()`](https://rdrr.io/pkg/MASS/man/predict.lda.html):

- `"plug-in"`,
- `"predictive"`,
- `"debiased"`.

The default method is `"plug-in"`.

``` r

pred_plugin <- predict(lda_fit, test, method = "plug-in")
pred_predictive <- predict(lda_fit, test, method = "predictive")
pred_debiased <- predict(lda_fit, test, method = "debiased")
```

The `"plug-in"` method uses estimated parameters directly in the
discriminant rule.

The `"predictive"` and `"debiased"` methods are alternative LDA
prediction rules following the
[`MASS::predict.lda()`](https://rdrr.io/pkg/MASS/man/predict.lda.html)
interface. See
[`?MASS::predict.lda`](https://rdrr.io/pkg/MASS/man/predict.lda.html)
for the original description of these prediction rules.

For easy data sets, the predicted classes may be identical across
methods.

``` r

c(
  plugin = mean(pred_plugin$class == test$Species),
  predictive = mean(pred_predictive$class == test$Species),
  debiased = mean(pred_debiased$class == test$Species)
)
#>     plugin predictive   debiased 
#>  0.9777778  0.9555556  0.9777778
```

Differences may be more visible in posterior probabilities.

``` r

head(pred_plugin$posterior)
#>    setosa   versicolor    virginica
#> 6       1 7.906717e-22 4.776377e-41
#> 9       1 2.674459e-16 1.730922e-35
#> 14      1 1.353910e-20 1.978581e-41
#> 16      1 3.431995e-28 5.926915e-49
#> 17      1 3.146059e-24 3.515450e-44
#> 19      1 8.270729e-22 2.610171e-41
head(pred_predictive$posterior)
#>    setosa   versicolor    virginica
#> 6       1 4.178535e-15 4.768416e-23
#> 9       1 3.998704e-12 5.873273e-21
#> 14      1 1.559157e-14 2.975005e-23
#> 16      1 4.033632e-17 1.834581e-24
#> 17      1 4.354432e-16 7.036180e-24
#> 19      1 3.498899e-15 2.957633e-23
head(pred_debiased$posterior)
#>    setosa   versicolor    virginica
#> 6       1 3.300984e-21 7.328371e-40
#> 9       1 7.678181e-16 1.822557e-34
#> 14      1 5.199420e-20 3.115444e-40
#> 16      1 2.204574e-27 1.553370e-47
#> 17      1 1.545286e-23 6.668687e-43
#> 19      1 3.448387e-21 4.076583e-40
```

## Leave-one-out prediction for QDA models

For
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
objects, leave-one-out cross-validation can be requested by omitting
`newdata` and using `method = "looCV"`.

In leave-one-out prediction, each training observation is classified as
if it had not been used to fit the model. This provides an internal
estimate of classification performance without creating a separate test
set.

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

Use leave-one-out results as an internal diagnostic, not as a
replacement for a proper independent test set when one is available.

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

The printed output shows the model call, prior probabilities, group
means, and information about selected or averaged permutation
structures.

For a structural view of the object, use
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

Fitted models are list-like S3 objects, so you can inspect their
component names.

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

The exact set of components depends on the model family, but the most
useful components are usually:

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

Inspect the optimization information directly.

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

A compact helper can be useful when debugging fitted objects.

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
#> 9               call           call      3      
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
#> 8               call           call      3          
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
#> 8               call           call      3          
#> 9  optimization_info      gips_perm      3          
#> 10             terms terms, formula      3          
#> 11           xlevels           list      0
```

## LDA diagnostics

For
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md)
objects, standard LDA-style diagnostics are available.

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

![](advanced-usage_files/figure-html/unnamed-chunk-48-1.png)

``` r

pairs(lda_fit, type = "std")
```

![](advanced-usage_files/figure-html/unnamed-chunk-49-1.png)

These methods are useful for inspecting the fitted discriminant
directions and class separation.

## Practical workflow

A typical workflow is:

1.  prepare numeric predictors and a categorical response,
2.  split data into training and test sets,
3.  start with
    [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
4.  compare with
    [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md)
    or
    [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
    if class-specific covariance structure may matter,
5.  inspect `optimization_info`,
6.  evaluate performance on held-out data.

``` r

fit_lda <- gipslda(Species ~ ., data = train)
fit_qda <- gipsqda(Species ~ ., data = train)
fit_joint <- gipsmultqda(Species ~ ., data = train)

pred_lda <- predict(fit_lda, test)
pred_qda <- predict(fit_qda, test)
pred_joint <- predict(fit_joint, test)

c(
  gipslda = mean(pred_lda$class == test$Species),
  gipsqda = mean(pred_qda$class == test$Species),
  gipsmultqda = mean(pred_joint$class == test$Species)
)
#>     gipslda     gipsqda gipsmultqda 
#>   0.9777778   0.9777778   1.0000000
```

## Troubleshooting

### The model is slow

Use `optimizer = "MH"` for larger numbers of features.

``` r

fit <- gipsqda(
  Species ~ .,
  data = train,
  optimizer = "MH",
  max_iter = 1000
)
```

Decrease `max_iter` for faster exploratory runs, then increase it for
final runs.

### The output selects `()`

The permutation `()` is the identity permutation. It means that the
selected structure did not impose a non-trivial permutation symmetry.

This can happen when:

- the data do not contain strong exchangeability patterns,
- the number of observations is too small,
- predictors are not meaningfully comparable,
- the optimizer did not find a more structured permutation.

### Predictions are identical across models

This may happen on simple data sets. Compare posterior probabilities and
inspect `optimization_info` for more detail.

``` r

head(pred_lda$posterior)
#>    setosa   versicolor    virginica
#> 6       1 7.906717e-22 4.776377e-41
#> 9       1 2.674459e-16 1.730922e-35
#> 14      1 1.353910e-20 1.978581e-41
#> 16      1 3.431995e-28 5.926915e-49
#> 17      1 3.146059e-24 3.515450e-44
#> 19      1 8.270729e-22 2.610171e-41
head(pred_qda$posterior)
#>    setosa   versicolor    virginica
#> 6       1 2.909134e-18 1.033310e-27
#> 9       1 1.075425e-13 1.166576e-19
#> 14      1 2.143385e-16 2.068832e-22
#> 16      1 2.508036e-24 3.982454e-36
#> 17      1 4.556129e-22 4.934051e-32
#> 19      1 4.873816e-18 2.042588e-29
head(pred_joint$posterior)
#>    setosa   versicolor    virginica
#> 6       1 1.809368e-18 3.630155e-32
#> 9       1 4.185357e-13 9.716563e-23
#> 14      1 8.534255e-16 4.000877e-26
#> 16      1 2.045692e-24 3.136984e-41
#> 17      1 1.110504e-21 8.419286e-37
#> 19      1 2.676609e-18 3.735764e-34
```

### Training and test data have different preprocessing

When preprocessing uses estimated quantities, such as means and standard
deviations for scaling, estimate them on the training data only and
apply the same transformation to test data.

## Summary

The main advanced controls are:

| Argument | Use |
|----|----|
| `MAP = TRUE` | use the single Maximum A Posteriori permutation |
| `MAP = FALSE` | average covariance projections using posterior probabilities |
| `optimizer = "BF"` | exhaustive search, typically for `p <= 10` |
| `optimizer = "MH"` | stochastic search, typically for `p > 10` |
| `max_iter` | number of Metropolis-Hastings iterations |
| `prior` | manually set class prior probabilities |
| `weighted_avg` | choose the pooled covariance estimator in [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) |

The three main models are:

| Function | Use when |
|----|----|
| [`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md) | classes can share one covariance structure |
| [`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md) | each class may have its own covariance structure |
| [`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md) | classes may have different covariance matrices but one shared permutation structure |
