# Predict from a gips LDA model

Computes class assignments, posterior probabilities, and linear
discriminant coordinates for new observations.

## Usage

``` r
# S3 method for class 'gipslda'
predict(
  object,
  newdata,
  prior = object$prior,
  dimen,
  method = c("plug-in", "predictive", "debiased"),
  ...
)
```

## Arguments

- object:

  A fitted `"gipslda"` object.

- newdata:

  An optional matrix or data frame of observations. If omitted, the
  training data are reconstructed from the fitted call.

- prior:

  Prior class probabilities used for prediction.

- dimen:

  Number of discriminant dimensions to use. The default uses all
  dimensions available in `object`.

- method:

  Prediction rule: `"plug-in"`, `"predictive"`, or `"debiased"`.

- ...:

  Further arguments passed to or from methods.

## Value

A list with `class` (predicted classes), `posterior` (posterior class
probabilities), and `x` (linear discriminant coordinates).

## See also

[`gipslda`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
[`predict.lda`](https://rdrr.io/pkg/MASS/man/predict.lda.html)
