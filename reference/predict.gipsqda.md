# Predict from a gips QDA model

Computes class assignments and posterior probabilities from a fitted
`"gipsqda"` model.

## Usage

``` r
# S3 method for class 'gipsqda'
predict(
  object,
  newdata,
  prior = object$prior,
  method = c("plug-in", "predictive", "debiased", "looCV"),
  ...
)
```

## Arguments

- object:

  A fitted `"gipsqda"` object.

- newdata:

  An optional matrix or data frame of observations. Omit this argument
  when using `method = "looCV"`.

- prior:

  Prior class probabilities used for prediction.

- method:

  Prediction rule: `"plug-in"`, `"predictive"`, `"debiased"`, or
  `"looCV"`. Leave-one-out prediction classifies the reconstructed
  training observations.

- ...:

  Further arguments passed to or from methods.

## Value

A list with `class` (predicted classes) and `posterior` (posterior class
probabilities).

## See also

[`gipsqda`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md),
[`predict.qda`](https://rdrr.io/pkg/MASS/man/predict.qda.html)
