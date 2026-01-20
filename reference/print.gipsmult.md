# Printing `gipsmult` object

Printing function for a `gipsmult` class.

## Usage

``` r
# S3 method for class 'gipsmult'
print(
  x,
  digits = 3,
  compare_to_original = TRUE,
  log_value = FALSE,
  oneline = FALSE,
  ...
)
```

## Arguments

- x:

  An object of a `gipsmult` class.

- digits:

  The number of digits after the comma for a posteriori to be presented.
  It can be negative. By default, `Inf`. It is passed to
  [`base::round()`](https://rdrr.io/r/base/Round.html).

- compare_to_original:

  A logical. Whether to print how many times more likely is the current
  permutation compared to:

  - the identity permutation `()` (for unoptimized `gipsmult` object);

  - the starting permutation (for optimized `gipsmult` object).

- log_value:

  A logical. Whether to print the logarithmic value. Default to `FALSE`.

- oneline:

  A logical. Whether to print in one or multiple lines. Default to
  `FALSE`.

- ...:

  The additional arguments passed to
  [`base::cat()`](https://rdrr.io/r/base/cat.html).

## Value

Returns an invisible `NULL`.

## See also

- [`find_MAP()`](https://AntoniKingston.github.io/gipsDA/reference/find_MAP.md) -
  The function that makes an optimized `gipsmult` object out of the
  unoptimized one.

## Examples

``` r
Ss <- list(matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE),
matrix(c(2, 1, 3, 7), nrow = 2, byrow = TRUE))
noo <- c(10,13)
g <- gipsmult(Ss, noo, perm = "(12)")
print(g, digits = 4, oneline = TRUE)
#> The permutation (1,2): is 0.5413 times more likely than the () permutation.
```
