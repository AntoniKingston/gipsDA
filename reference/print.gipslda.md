# Print a fitted gipsDA model

Prints the main components of a fitted `gipsDA` model, including the
model call, fitting options, group means, class counts, selected MAP
permutation, and posterior probabilities of retained permutations when
stored.

## Usage

``` r
# S3 method for class 'gipslda'
print(x, ...)

# S3 method for class 'gipsqda'
print(x, ...)

# S3 method for class 'gipsmultqda'
print(x, ...)
```

## Arguments

- x:

  A fitted gipsDA model.

- ...:

  Further arguments passed to printing methods.

## Value

Invisibly returns `x`.

## Examples

``` r
fit <- gipslda(Species ~ ., data = iris, optimizer = "BF")
print(fit)
#> Call:
#> gipslda(Species ~ ., data = iris, optimizer = "BF")
#> 
#> Model: gipslda 
#> Number of observations: 150 
#> Number of groups: 3 
#> Number of predictors: 4 
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
#> $store_probabilities
#> [1] TRUE
#> 
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Class counts:
#>     setosa versicolor  virginica 
#>         50         50         50 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa            5.006       3.428        1.462       0.246
#> versicolor        5.936       2.770        4.260       1.326
#> virginica         6.588       2.974        5.552       2.026
#> 
#> Selected MAP permutation: (1,3)(2,4) 
#> 
#> Posterior probabilities of retained permutations:
#> (1,3)(2,4) 
#>  0.9995819 
#> 
#> Coefficients of linear discriminants:
#>                     LD1        LD2
#> Sepal.Length  0.8664772 -0.1303263
#> Sepal.Width   1.5080574 -2.0747392
#> Petal.Length -2.2199033  1.0452646
#> Petal.Width  -2.7397901 -2.9835335
#> 
#> Proportion of trace:
#>   LD1   LD2 
#> 0.991 0.009 
```
