make_four_group_lda_fixture <- function() {
  within <- cbind(
    x1 = c(-1.0, -0.4, 0.1, 0.5, 0.8, 1.2),
    x2 = c(0.3, -1.1, 0.7, -0.2, 1.4, -0.6),
    x3 = c(-0.8, 0.9, -0.3, 1.1, -1.2, 0.4),
    x4 = c(1.0, -0.7, -1.1, 0.2, 0.6, -0.1)
  )
  offsets <- rbind(
    c(0, 0, 0, 0),
    c(3, 0, 0, 0),
    c(0, 3, 0, 0),
    c(0, 0, 3, 0)
  )
  x <- do.call(rbind, lapply(seq_len(4), function(i) {
    sweep(within, 2, offsets[i, ], FUN = "+")
  }))
  grouping <- factor(rep(letters[1:4], each = nrow(within)))
  data <- data.frame(class = grouping, x, check.names = FALSE)
  list(x = x, grouping = grouping, data = data)
}

test_that("predict supports every discrimination method", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")

  for (method in c("plug-in", "predictive", "debiased")) {
    prediction <- predict(fit, fixture$x[1:7, , drop = FALSE], method = method)
    expect_named(prediction, c("class", "posterior", "x"))
    expect_s3_class(prediction$class, "factor")
    expect_equal(levels(prediction$class), fit$lev)
    expect_length(prediction$class, 7)
    expect_valid_posterior(prediction$posterior, n = 7, groups = 3)
    expect_equal(dim(prediction$x), c(7, length(fit$svd)))
  }
})

test_that("predict accepts matrix, data frame, row vector, and omitted newdata", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)
  fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")

  matrix_prediction <- predict(fit, fixture$x[1:4, , drop = FALSE])
  frame_prediction <- predict(
    fit,
    as.data.frame(fixture$x[1:4, , drop = FALSE])
  )
  row_prediction <- predict(fit, fixture$x[1, ])
  training_prediction <- predict(fit)

  expect_equal(frame_prediction$posterior, matrix_prediction$posterior)
  expect_length(row_prediction$class, 1)
  expect_valid_posterior(row_prediction$posterior, n = 1, groups = 3)
  expect_length(training_prediction$class, nrow(fixture$x))
  expect_valid_posterior(
    training_prediction$posterior,
    n = nrow(fixture$x),
    groups = 3
  )
})

test_that("formula prediction reconstructs training and new model frames", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)
  fit <- gipslda(class ~ x1 + x2, fixture$data)

  from_frame <- predict(fit, fixture$data[1:5, ])
  from_training <- predict(fit)

  expect_length(from_frame$class, 5)
  expect_equal(dim(from_frame$posterior), c(5, 3))
  expect_true(all(is.finite(from_frame$posterior)))
  expect_equal(unname(rowSums(from_frame$posterior)), rep(1, 5), tolerance = 1e-8)
  expect_length(from_training$class, nrow(fixture$data))
  expect_equal(rownames(from_frame$posterior), rownames(fixture$data)[1:5])
})

test_that("predict honors dimensions and alternative priors", {
  fixture <- make_classification_fixture(p = 3, n_per_class = 8)
  fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")
  prior <- c(setosa = 0.7, versicolor = 0.2, virginica = 0.1)

  one_dimension <- predict(
    fit,
    fixture$x[1:5, , drop = FALSE],
    prior = prior,
    dimen = 1
  )
  capped_dimension <- predict(
    fit,
    fixture$x[1:5, , drop = FALSE],
    dimen = 100
  )

  expect_equal(ncol(one_dimension$x), 1)
  expect_equal(ncol(capped_dimension$x), length(fit$svd))
  expect_valid_posterior(one_dimension$posterior, n = 5, groups = 3)
  expect_gt(mean(one_dimension$posterior[, "setosa"]), 0)
})

test_that("predict validates object, variables, priors, and method", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)
  fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")
  predict_method <- getS3method("predict", "gipslda")

  expect_error(
    predict_method(list(), fixture$x),
    "object not of class \"gipslda\"",
    fixed = TRUE
  )
  expect_error(
    predict(fit, fixture$x[, 1, drop = FALSE]),
    "wrong number of variables",
    fixed = TRUE
  )
  renamed <- fixture$x
  colnames(renamed) <- c("other1", "other2")
  expect_warning(
    predict(fit, renamed),
    "variable names in 'newdata' do not match"
  )
  expect_error(
    predict(fit, fixture$x, prior = c(0.5, 0.4, 0.1, 0)),
    "incorrect length"
  )
  expect_error(
    predict(fit, fixture$x, prior = c(0.5, 0.4, -0.1)),
    "invalid 'prior'",
    fixed = TRUE
  )
  expect_error(
    predict(fit, fixture$x, method = "unknown"),
    "'arg' should be one of"
  )
})

test_that("print, coef, and model.frame methods expose fitted content", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)
  matrix_fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")
  formula_fit <- gipslda(class ~ ., fixture$data)

  output <- capture.output(returned <- print(matrix_fit))
  expect_identical(returned, matrix_fit)
  expect_true(any(grepl("^Call:", output)))
  expect_true(any(grepl("Prior probabilities", output)))
  expect_true(any(grepl("Group means", output)))
  expect_true(any(grepl("Coefficients of linear discriminants", output)))
  expect_true(any(grepl(
    "Posterior probabilities of retained permutations:",
    output,
    fixed = TRUE
  )))

  expect_identical(coef(matrix_fit), matrix_fit$scaling)
  reconstructed <- model.frame(formula_fit)
  expect_equal(
    unname(as.data.frame(reconstructed)),
    unname(fixture$data),
    ignore_attr = TRUE
  )
  expect_equal(
    unname(model.response(reconstructed)),
    unname(fixture$grouping)
  )
})

test_that("plot.gipslda handles one, two, and more dimensions", {
  local_plot_device()
  one <- make_binary_fixture(p = 1, n_per_class = 6)
  two <- make_classification_fixture(p = 2, n_per_class = 6)
  many <- make_four_group_lda_fixture()

  fit_one <- gipslda(one$x, one$grouping, optimizer = "BF")
  fit_two <- gipslda(two$x, two$grouping, optimizer = "BF")
  fit_many <- gipslda(many$x, many$grouping, optimizer = "BF")

  expect_equal(ncol(fit_one$scaling), 1)
  expect_equal(ncol(fit_two$scaling), 2)
  expect_gt(ncol(fit_many$scaling), 2)
  expect_invisible(plot(fit_one))
  expect_invisible(plot(fit_two, abbrev = 1))
  expect_invisible(plot(fit_many))
  expect_invisible(plot(fit_many, dimen = 2))
})

test_that("pairs.gipslda supports standard and trellis displays", {
  local_plot_device()
  fixture <- make_four_group_lda_fixture()
  fit <- gipslda(class ~ ., fixture$data)

  expect_invisible(pairs(fit, type = "std", abbrev = 1))
  expect_invisible(pairs(fit, type = "std", dimen = 2))
  expect_invisible(pairs(fit, type = "trellis"))
  expect_error(pairs(fit, type = "invalid"), "'arg' should be one of")
})

test_that("ldahist supports histogram, density, and combined displays", {
  local_plot_device()
  data <- c(-2.0, -1.3, -0.8, -0.2, 0.1, 0.5, 1.1, 1.6, 2.0, 2.4)
  grouping <- factor(rep(c("low", "high"), each = 5))
  histogram <- gipsDA:::ldahist

  expect_invisible(histogram(data, grouping, type = "histogram"))
  expect_invisible(histogram(data, grouping, type = "density", width = 0.4))
  expect_invisible(
    histogram(data, grouping, type = "both", width = 0.4, sep = FALSE)
  )
})

test_that("ldahist validates type and histogram breaks", {
  local_plot_device()
  data <- c(-2, -1, 0, 1, 2, 3)
  grouping <- factor(rep(c("a", "b"), each = 3))
  histogram <- gipsDA:::ldahist

  expect_error(
    histogram(data, grouping, type = "invalid"),
    "'arg' should be one of"
  )
  expect_error(
    histogram(data, grouping, breaks = c(-3, 0, 0, 4)),
    "'breaks' must be strictly increasing",
    fixed = TRUE
  )
  expect_error(
    histogram(data, grouping, breaks = c(-1, 0, 1)),
    "'breaks' do not cover the data",
    fixed = TRUE
  )
})
