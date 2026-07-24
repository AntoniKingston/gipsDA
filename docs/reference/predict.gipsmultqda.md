# Predict from a joint-projection gips QDA model

Computes class assignments and posterior probabilities from a fitted
`"gipsmultqda"` model. Prediction rules match those of
[`predict.gipsqda`](https://antonikingston.github.io/gipsDA/reference/predict.gipsqda.md);
the fitted covariance estimates differ because they were projected
jointly.

## Usage

``` r
# S3 method for class 'gipsmultqda'
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

  A fitted `"gipsmultqda"` object.

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

[`gipsmultqda`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md),
[`predict.gipsqda`](https://antonikingston.github.io/gipsDA/reference/predict.gipsqda.md)
