# Plot a fitted gips LDA model

Displays training observations in discriminant space. One discriminant
dimension is shown as class-wise histograms or densities, two dimensions
as an equal-scale scatter plot, and higher dimensions as a pairs plot.

## Usage

``` r
# S3 method for class 'gipslda'
plot(
  x,
  panel = panel.gipslda,
  ...,
  cex = 0.7,
  dimen,
  abbrev = FALSE,
  xlab = "LD1",
  ylab = "LD2"
)
```

## Arguments

- x:

  A fitted `"gipslda"` object.

- panel:

  Panel function used for scatter plots.

- ...:

  Additional graphical arguments.

- cex:

  Character expansion used by the default panel.

- dimen:

  Maximum number of discriminant dimensions to display.

- abbrev:

  Logical or integer controlling abbreviation of class labels.

- xlab, ylab:

  Axis labels.

## Value

Invisibly returns `NULL`.

## See also

[`gipslda`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
[`pairs.gipslda`](https://antonikingston.github.io/gipsDA/reference/pairs.gipslda.md)
