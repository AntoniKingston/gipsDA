# The constructor of a `gipsmult` class.

Create a `gipsmult` object. This object will contain initial data and
all other information needed to find the most likely invariant
permutation. It will not perform optimization. One must call the
[`find_MAP()`](https://AntoniKingston.github.io/gipsDA/reference/find_MAP.md)
function to do it. See the examples below.

## Usage

``` r
gipsmult(
  Ss,
  numbers_of_observations,
  delta = 3,
  D_matrices = NULL,
  was_mean_estimated = TRUE,
  perm = ""
)

new_gipsmult(
  list_of_gips_perm,
  Ss,
  numbers_of_observations,
  delta,
  D_matrices,
  was_mean_estimated,
  optimization_info
)
```

## Arguments

- delta:

  A number, hyper-parameter of a Bayesian model. It has to be strictly
  bigger than 1. See the **Hyperparameters** section below.

- was_mean_estimated:

  A boolean.

  - Set `TRUE` (default) when your `S` parameter is a result of a
    [`stats::cov()`](https://rdrr.io/r/stats/cor.html) function.

  - Set FALSE when your `S` parameter is a result of a
    `(t(Z) %*% Z) / number_of_observations` calculation.

- perm:

  An optional permutation to be the base for the `gipsmult` object. It
  can be of a `gips_perm` or a `permutation` class, or anything the
  function
  [`permutations::permutation()`](https://robinhankin.github.io/permutations/reference/permutation.html)
  can handle. It can also be of a `gipsmult` class, but it will be
  interpreted as the underlying `gips_perm`.

- list_of_gips_perm:

  A list with a single element of a `gips_perm` class. The base object
  for the `gipsmult` object.

- optimization_info:

  For internal use only. `NULL` or the list with information about the
  optimization process.

- S:

  A matrix; empirical covariance matrix. When `Z` is the observed data:

  - if one does not know the theoretical mean and has to estimate it
    with the observed mean, use `S = cov(Z)`, and leave parameter
    `was_mean_estimated = TRUE` as default;

  - if one know the theoretical mean is 0, use
    `S = (t(Z) %*% Z) / number_of_observations`, and set parameter
    `was_mean_estimated = FALSE`.

- number_of_observations:

  A number of data points that `S` is based on.

- D_matrix:

  Symmetric, positive-definite matrix of the same size as `S`.
  Hyper-parameter of a Bayesian model. When `NULL`, the (hopefully)
  reasonable one is derived from the data. For more details, see the
  **Hyperparameters** section below.

## Value

`gipsmult()` returns an object of a `gipsmult` class after the safety
checks.

`new_gipsmult()` returns an object of a `gipsmult` class without the
safety checks.

## Functions

- `new_gipsmult()`: Constructor. It is only intended for low-level use.

## Methods for a `gipsmult` class

- [`plot.gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/plot.gipsmult.md)

- [`print.gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/print.gipsmult.md)

## Hyperparameters

We encourage the user to try `D_matrix = d * I`, where `I` is an
identity matrix of a size `p x p` and `d > 0` for some different `d`.
When `d` is small compared to the data (e.g.,
`d = 0.1 * mean(diag(S))`), bigger structures will be found. When `d` is
big compared to the data (e.g., `d = 100 * mean(diag(S))`), the
posterior distribution does not depend on the data.

Taking `D_matrix = d * I` is equivalent to setting `S <- S / d`.

The default for `D_matrix` is `D_matrix = d * I`, where
`d = mean(diag(S))`, which is equivalent to modifying `S` so that the
mean value on the diagonal is 1.

In the Bayesian model, the prior distribution for the covariance matrix
is a generalized case of [Wishart
distribution](https://en.wikipedia.org/wiki/Wishart_distribution).

## See also

- [`stats::cov()`](https://rdrr.io/r/stats/cor.html) – The `Ss`
  parameter, as a list of empirical covariance matrices, is most of the
  time a result of the [`cov()`](https://rdrr.io/r/stats/cor.html)
  function. For more information, see [Wikipedia - Estimation of
  covariance
  matrices](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices).

- [`find_MAP()`](https://AntoniKingston.github.io/gipsDA/reference/find_MAP.md)
  – The function that finds the Maximum A Posteriori (MAP) Estimator for
  a given `gipsmult` object.

- [`gips::gips_perm()`](https://przechoj.github.io/gips/reference/gips_perm.html)
  – The constructor of a `gips_perm` class. The `gips_perm` object is
  used as the base object for the `gipsmult` object.

## Examples

``` r
perm_size <- 5
numbers_of_observations <- c(15, 18, 19)
Sigma <- diag(rep(1, perm_size))
n_matrices <- 3
df <- 20
Ss <- rWishart(n = n_matrices, df = df, Sigma = Sigma)
Ss <- lapply(1:n_matrices, function(x) Ss[, , x])
g <- gipsmult(Ss, numbers_of_observations)

g_map <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
#> Error in find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force"): could not find function "find_MAP"
g_map
#> Error: object 'g_map' not found

print(g_map)
#> Error: object 'g_map' not found

if (require("graphics")) {
  plot(g_map, type = "MLE", logarithmic_x = TRUE)
}
#> Error: object 'g_map' not found
```
