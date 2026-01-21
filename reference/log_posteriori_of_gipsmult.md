# A log of a posteriori that the covariance matrix is invariant under permutation

More precisely, it is the logarithm of an unnormalized posterior
probability. It is the goal function for optimization algorithms in the
[`find_MAP()`](https://AntoniKingston.github.io/gipsDA/reference/find_MAP.md)
function. The `perm_proposal` that maximizes this function is the
Maximum A Posteriori (MAP) Estimator.

## Usage

``` r
log_posteriori_of_gipsmult(g)
```

## Arguments

- g:

  An object of a `gipsmult` class.

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

- [`find_MAP()`](https://AntoniKingston.github.io/gipsDA/reference/find_MAP.md) -
  The function that optimizes the `log_posteriori_of_gips` function.

- [`gips::compare_posteriories_of_perms()`](https://przechoj.github.io/gips/reference/compare_posteriories_of_perms.html) -
  Uses `log_posteriori_of_gips()` to compare a posteriori of two
  permutations.

## Examples

``` r
# In the space with p = 2, there is only 2 permutations:
perm1 <- permutations::as.cycle("(1)(2)")
perm2 <- permutations::as.cycle("(1,2)")
S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
S2 <- matrix(c(2, 1, 3, 7), nrow = 2, byrow = TRUE)
g1 <- gipsmult(list(S1, S2), c(100, 100), perm = perm1)
g2 <- gipsmult(list(S1, S2), c(100, 100), perm = perm2)
log_posteriori_of_gipsmult(g1) # -354.4394, this is the MAP Estimator
#> [1] -354.4394
log_posteriori_of_gipsmult(g2) # -380.0079
#> [1] -380.0079

exp(log_posteriori_of_gipsmult(g1) - log_posteriori_of_gipsmult(g2)) # 127131902082
#> [1] 127131902082
# g1 is 127131902082 times more likely than g2.
# This is the expected outcome because S1[1,1] and S2[1,1]
# differ significantly from S1[2,2] and S2[2,2] respectively.

# ========================================================================

S3 <- matrix(c(1, 0.5, 0.5, 1.1), nrow = 2, byrow = TRUE)
S4 <- matrix(c(2, 1, 3, 2.137), nrow = 2, byrow = TRUE)
g1 <- gipsmult(list(S3, S4), c(100, 100), perm = perm1)
g2 <- gipsmult(list(S3, S4), c(100, 100), perm = perm2)
log_posteriori_of_gipsmult(g1) # -148.6485
#> Warning: ℹ `project_matrix()` is designed for positive semi-definite matrices
#> ✖ You provided `S` that is not positive semi-definite matrix
#> • `gips` can still project this matrix on the provided permutation
#> ℹ Did You provide the wrong `S` matrix?
#> [1] -148.6485
log_posteriori_of_gipsmult(g2) # -145.3019, this is the MAP Estimator
#> Warning: ℹ `project_matrix()` is designed for positive semi-definite matrices
#> ✖ You provided `S` that is not positive semi-definite matrix
#> • `gips` can still project this matrix on the provided permutation
#> ℹ Did You provide the wrong `S` matrix?
#> [1] -145.3019

exp(log_posteriori_of_gipsmult(g2) - log_posteriori_of_gipsmult(g1)) # 28.406
#> Warning: ℹ `project_matrix()` is designed for positive semi-definite matrices
#> ✖ You provided `S` that is not positive semi-definite matrix
#> • `gips` can still project this matrix on the provided permutation
#> ℹ Did You provide the wrong `S` matrix?
#> Warning: ℹ `project_matrix()` is designed for positive semi-definite matrices
#> ✖ You provided `S` that is not positive semi-definite matrix
#> • `gips` can still project this matrix on the provided permutation
#> ℹ Did You provide the wrong `S` matrix?
#> [1] 28.406
# g2 is 28.406 times more likely than g1.
# This is the expected outcome because S1[1,1] and S2[1,1]
# are very close to S1[2,2] and S2[2,2] respectively.
```
