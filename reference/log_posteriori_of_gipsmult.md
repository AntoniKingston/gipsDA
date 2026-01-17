# A log of a posteriori that the covariance matrix is invariant under permutation

More precisely, it is the logarithm of an unnormalized posterior
probability. It is the goal function for optimization algorithms in the
`find_MAP()` function. The `perm_proposal` that maximizes this function
is the Maximum A Posteriori (MAP) Estimator.

## Usage

``` r
log_posteriori_of_gipsmult(g)
```

## Arguments

- g:

  An object of a `gips` class.

## Value

Returns a value of the logarithm of an unnormalized A Posteriori.

## Details

It is calculated using [formulas (33) and (27) from
references](https://arxiv.org/abs/2004.03503).

If `Inf` or `NaN` is reached, it produces a warning.

## References

Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam. "Model
selection in the space of Gaussian models invariant by symmetry." The
Annals of Statistics, 50(3) 1747-1774 June 2022. [arXiv
link](https://arxiv.org/abs/2004.03503);
[doi:10.1214/22-AOS2174](https://doi.org/10.1214/22-AOS2174)

## See also

- `calculate_gamma_function()` - The function that calculates the value
  needed for `log_posteriori_of_gips()`.

- `get_structure_constants()` - The function that calculates the
  structure constants needed for `log_posteriori_of_gips()`.

- `find_MAP()` - The function that optimizes the
  `log_posteriori_of_gips` function.

- `compare_posteriories_of_perms()` - Uses `log_posteriori_of_gips()` to
  compare a posteriori of two permutations.

- [`vignette("Theory", package = "gips")`](https://przechoj.github.io/gips/articles/Theory.html)
  or its [pkgdown
  page](https://przechoj.github.io/gips/articles/Theory.html) - A place
  to learn more about the math behind the `gips` package.

## Examples

``` r
# In the space with p = 2, there is only 2 permutations:
perm1 <- permutations::as.cycle("(1)(2)")
perm2 <- permutations::as.cycle("(1,2)")
S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
g1 <- gips(S1, 100, perm = perm1)
#> Error in gips(S1, 100, perm = perm1): could not find function "gips"
g2 <- gips(S1, 100, perm = perm2)
#> Error in gips(S1, 100, perm = perm2): could not find function "gips"
log_posteriori_of_gips(g1) # -134.1615, this is the MAP Estimator
#> Error in log_posteriori_of_gips(g1): could not find function "log_posteriori_of_gips"
log_posteriori_of_gips(g2) # -138.1695
#> Error in log_posteriori_of_gips(g2): could not find function "log_posteriori_of_gips"

exp(log_posteriori_of_gips(g1) - log_posteriori_of_gips(g2)) # 55.0
#> Error in log_posteriori_of_gips(g1): could not find function "log_posteriori_of_gips"
# g1 is 55 times more likely than g2.
# This is the expected outcome because S[1,1] significantly differs from S[2,2].

compare_posteriories_of_perms(g1, g2)
#> Error in compare_posteriories_of_perms(g1, g2): could not find function "compare_posteriories_of_perms"
# The same result, but presented in a more pleasant way

# ========================================================================

S2 <- matrix(c(1, 0.5, 0.5, 1.1), nrow = 2, byrow = TRUE)
g1 <- gips(S2, 100, perm = perm1)
#> Error in gips(S2, 100, perm = perm1): could not find function "gips"
g2 <- gips(S2, 100, perm = perm2)
#> Error in gips(S2, 100, perm = perm2): could not find function "gips"
log_posteriori_of_gips(g1) # -98.40984
#> Error in log_posteriori_of_gips(g1): could not find function "log_posteriori_of_gips"
log_posteriori_of_gips(g2) # -95.92039, this is the MAP Estimator
#> Error in log_posteriori_of_gips(g2): could not find function "log_posteriori_of_gips"

exp(log_posteriori_of_gips(g2) - log_posteriori_of_gips(g1)) # 12.05
#> Error in log_posteriori_of_gips(g2): could not find function "log_posteriori_of_gips"
# g2 is 12 times more likely than g1.
# This is the expected outcome because S[1,1] is very close to S[2,2].

compare_posteriories_of_perms(g2, g1)
#> Error in compare_posteriories_of_perms(g2, g1): could not find function "compare_posteriories_of_perms"
# The same result, but presented in a more pleasant way
```
