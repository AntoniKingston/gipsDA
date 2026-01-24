# Find the Maximum A Posteriori Estimation

Use one of the optimization algorithms to find the permutation that
maximizes a posteriori probability based on observed data. Not all
optimization algorithms will always find the MAP, but they try to find a
significant value.

## Usage

``` r
find_MAP(
  g,
  max_iter = NA,
  optimizer = NA,
  show_progress_bar = TRUE,
  save_all_perms = FALSE,
  return_probabilities = FALSE
)
```

## Arguments

- g:

  Object of a `gipsmult` class.

- max_iter:

  The number of iterations for an algorithm to perform. At least 2. For
  `optimizer = "BF"`, it is not used; for `optimizer = "MH"`, it has to
  be finite; for `optimizer = "HC"`, it can be infinite.

- optimizer:

  The optimizer for the search of the maximum posteriori:

  - `"BF"` (the default for unoptimized `g` with `perm size <= 9`) -
    Brute Force;

  - `"MH"` (the default for unoptimized `g` with `perm size > 10`) -
    Metropolis-Hastings;

  - `"HC"` - Hill Climbing;

  - `"continue"` (the default for optimized `g`) - The same as the `g`
    was optimized by (see Examples).

  See the **Possible algorithms to use as optimizers** section below for
  more details.

- show_progress_bar:

  A boolean. Indicate whether or not to show the progress bar:

  - When `max_iter` is infinite, `show_progress_bar` has to be `FALSE`;

  - When `return_probabilities = TRUE`, then shows an additional
    progress bar for the time when the probabilities are calculated.

- save_all_perms:

  A boolean. `TRUE` indicates saving a list of all permutations visited
  during optimization. This can be useful sometimes but needs a lot more
  RAM.

- return_probabilities:

  A boolean. `TRUE` can only be provided only when
  `save_all_perms = TRUE`. For:

  - `optimizer = "MH"` - use Metropolis-Hastings results to estimate
    posterior probabilities;

  - `optimizer = "BF"` - use brute force results to calculate exact
    posterior probabilities.

  These additional calculations are costly, so a second and third
  progress bar is shown (when `show_progress_bar = TRUE`).

  To examine probabilities after optimization, call
  [`get_probabilities_from_gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/get_probabilities_from_gipsmult.md).

## Value

Returns an optimized object of a `gipsmult` class.

## Details

`find_MAP()` can produce a warning when:

- the optimizer "hill_climbing" gets to the end of its `max_iter`
  without converging.

- the optimizer will find the permutation with smaller `n0` than
  `number_of_observations`

## Possible algorithms to use as optimizers

For every algorithm, there are some aliases available.

- `"brute_force"`, `"BF"`, `"full"` - use the **Brute Force** algorithm
  that checks the whole permutation space of a given size. This
  algorithm will find the actual Maximum A Posteriori Estimation, but it
  is very computationally expensive for bigger spaces. We recommend
  Brute Force only for `p <= 9`.

- `"Metropolis_Hastings"`, `"MH"` - use the **Metropolis-Hastings**
  algorithm; [see
  Wikipedia](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm).
  The algorithm will draw a random transposition in every iteration and
  consider changing the current state (permutation). When the `max_iter`
  is reached, the algorithm will return the best permutation calculated
  as the MAP Estimator. This algorithm used in this context is a special
  case of the **Simulated Annealing** the user may be more familiar
  with; [see
  Wikipedia](https://en.wikipedia.org/wiki/Simulated_annealing).

- `"hill_climbing"`, `"HC"` - use the **hill climbing** algorithm; [see
  Wikipedia](https://en.wikipedia.org/wiki/Hill_climbing). The algorithm
  will check all transpositions in every iteration and go to the one
  with the biggest a posteriori value. The optimization ends when all
  *neighbors* will have a smaller a posteriori value. If the `max_iter`
  is reached before the end, then the warning is shown, and it is
  recommended to continue the optimization on the output of the
  `find_MAP()` with `optimizer = "continue"`; see examples. Remember
  that `p*(p-1)/2` transpositions will be checked in every iteration.
  For bigger `p`, this may be costly.

## See also

- [`gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/gipsmult.md) -
  The constructor of a `gipsmult` class. The `gipsmult` object is used
  as the `g` parameter of `find_MAP()`.

- [`plot.gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/plot.gipsmult.md) -
  Practical plotting function for visualizing the optimization process.

- [`get_probabilities_from_gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/get_probabilities_from_gipsmult.md) -
  When `find_MAP(return_probabilities = TRUE)` was called, probabilities
  can be extracted with this function.

- [`log_posteriori_of_gipsmult()`](https://AntoniKingston.github.io/gipsDA/reference/log_posteriori_of_gipsmult.md) -
  The function that the optimizers of `find_MAP()` tries to find the
  argmax of.

## Examples

``` r
require("MASS") # for mvrnorm()
#> Loading required package: MASS

perm_size <- 6
mu1 <- runif(6, -10, 10)
mu2 <- runif(6, -10, 10) # Assume we don't know the means
sigma1 <- matrix(
  data = c(
    1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
    0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
    0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
    0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
    0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
    0.8, 0.6, 0.4, 0.6, 0.8, 1.0
  ),
  nrow = perm_size, byrow = TRUE
)
sigma2 <- matrix(
  data = c(
    1.0, 0.5, 0.2, 0.0, 0.2, 0.5,
    0.5, 1.0, 0.5, 0.2, 0.0, 0.2,
    0.2, 0.5, 1.0, 0.5, 0.2, 0.0,
    0.0, 0.2, 0.5, 1.0, 0.5, 0.2,
    0.2, 0.0, 0.2, 0.5, 1.0, 0.5,
    0.5, 0.2, 0.0, 0.2, 0.5, 1.0
  ),
  nrow = perm_size, byrow = TRUE
)
# sigma1 and sigma2 are matrices invariant under permutation (1,2,3,4,5,6)
numbers_of_observations <- c(21, 37)
Z1 <- MASS::mvrnorm(numbers_of_observations[1], mu = mu1, Sigma = sigma1)
Z2 <- MASS::mvrnorm(numbers_of_observations[2], mu = mu2, Sigma = sigma2)
S1 <- cov(Z1)
S2 <- cov(Z2) # Assume we have to estimate the mean

g <- gipsmult(list(S1, S2), numbers_of_observations)

g_map <- find_MAP(g, max_iter = 5, show_progress_bar = FALSE, optimizer = "Metropolis_Hastings")
g_map
#> The permutation ():
#>  - was found after 5 posteriori calculations;
#>  - is 1 times more likely than the () permutation.

g_map2 <- find_MAP(g_map, max_iter = 5, show_progress_bar = FALSE, optimizer = "continue")

if (require("graphics")) {
  plot(g_map2, type = "both", logarithmic_x = TRUE)
}


g_map_BF <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
```
