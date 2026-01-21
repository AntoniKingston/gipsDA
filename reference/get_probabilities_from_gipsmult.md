# Extract probabilities for `gipsmult` object optimized with `return_probabilities = TRUE`

After the `gipsmult` object was optimized with the
`find_MAP(return_probabilities = TRUE)` function, then those calculated
probabilities can be extracted with this function.

## Usage

``` r
get_probabilities_from_gipsmult(g)
```

## Arguments

- g:

  An object of class `gipsmult`. A result of a
  `find_MAP(return_probabilities = TRUE)`.

## Value

Returns a numeric vector, calculated values of probabilities. Names
contain permutations this probabilities represent. For `gipsmult` object
optimized with `find_MAP(return_probabilities = FALSE)`, it returns a
`NULL` object. It is sorted according to the probability.

## See also

- [`find_MAP()`](https://AntoniKingston.github.io/gipsDA/reference/find_MAP.md) -
  The `get_probabilities_from_gipsmult()` is called on the output of
  `find_MAP(return_probabilities = TRUE, save_all_perms = TRUE)`.

## Examples

``` r
Ss <- list(
  matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE),
  matrix(c(2, 1, 3, 7), nrow = 2, byrow = TRUE)
)
noo <- c(10, 13)
g <- gipsmult(Ss, noo)
g_map <- find_MAP(g,
  optimizer = "BF", show_progress_bar = FALSE,
  return_probabilities = TRUE, save_all_perms = TRUE
)

get_probabilities_from_gipsmult(g_map)
#>        ()     (1,2) 
#> 0.6487995 0.3512005 
```
