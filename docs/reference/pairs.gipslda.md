# Pairs plot for a fitted gips LDA model

Displays discriminant coordinates using either base graphics or a
lattice scatterplot matrix.

## Usage

``` r
# S3 method for class 'gipslda'
pairs(
  x,
  labels = colnames(x),
  panel = panel.gipslda,
  dimen,
  abbrev = FALSE,
  ...,
  cex = 0.7,
  type = c("std", "trellis")
)
```

## Arguments

- x:

  A fitted `"gipslda"` object.

- labels:

  Labels for discriminant dimensions.

- panel:

  Panel function used by base graphics.

- dimen:

  Maximum number of discriminant dimensions to display.

- abbrev:

  Logical or integer controlling abbreviation of class labels.

- ...:

  Additional graphical arguments.

- cex:

  Character expansion used by the default panel.

- type:

  Either `"std"` for base graphics or `"trellis"` for a lattice display.

## Value

Invisibly returns `NULL`.

## See also

[`gipslda`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
[`plot.gipslda`](https://antonikingston.github.io/gipsDA/reference/plot.gipslda.md)
