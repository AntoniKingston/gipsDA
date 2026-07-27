test_that("the public API exports the three discriminant-analysis generics", {
  exports <- getNamespaceExports("gipsDA")

  expect_setequal(exports, c("gipslda", "gipsqda", "gipsmultqda"))
  expect_true(all(vapply(
    mget(exports, envir = asNamespace("gipsDA")),
    is.function,
    logical(1)
  )))
})

test_that("documented S3 methods are registered", {
  fit_classes <- c("gipslda", "gipsqda", "gipsmultqda")
  fit_inputs <- c("default", "matrix", "data.frame", "formula")

  for (generic in fit_classes) {
    for (input in fit_inputs) {
      expect_true(is.function(getS3method(generic, input)))
    }
    expect_true(is.function(getS3method("predict", generic)))
    expect_true(is.function(getS3method("print", generic)))
    expect_true(is.function(getS3method("model.frame", generic)))
  }

  expect_true(is.function(getS3method("coef", "gipslda")))
  expect_true(is.function(getS3method("plot", "gipslda")))
  expect_true(is.function(getS3method("pairs", "gipslda")))
})
