# Summarize a fitted gipsDA model

Creates a compact summary object for a fitted `gipsDA` model. The
summary contains the most important fitted quantities and can be printed
using the corresponding [`print()`](https://rdrr.io/r/base/print.html)
method.

## Usage

``` r
# S3 method for class 'gipslda'
summary(object, ...)

# S3 method for class 'gipsqda'
summary(object, ...)

# S3 method for class 'gipsmultqda'
summary(object, ...)
```

## Arguments

- object:

  A fitted gipsDA model.

- ...:

  Further arguments passed to or from methods.

## Value

An object of class `"summary.gipsda"`.

## Examples

``` r
fit <- gipslda(Species ~ ., data = iris, optimizer = "BF")
summary(fit)
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
#> Class counts:
#>     setosa versicolor  virginica 
#>         50         50         50 
#> 
#> Prior probabilities of groups:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Group means:
#>            Sepal.Length Sepal.Width Petal.Length Petal.Width
#> setosa            5.006       3.428        1.462       0.246
#> versicolor        5.936       2.770        4.260       1.326
#> virginica         6.588       2.974        5.552       2.026
#> 
#> Proportion of trace:
#>   LD1   LD2 
#> 0.991 0.009 
#> 
#> Selected MAP permutation: (1,3)(2,4) 
#> 
#> Posterior probabilities of retained permutations:
#> (1,3)(2,4) 
#>  0.9995819 

summary_object <- summary(fit)
names(summary_object)
#>  [1] "model"                    "call"                    
#>  [3] "n"                        "p"                       
#>  [5] "groups"                   "counts"                  
#>  [7] "prior"                    "means"                   
#>  [9] "fit_info"                 "scaling"                 
#> [11] "svd"                      "proportion_trace"        
#> [13] "optimization_info"        "selected_map_permutation"
```
