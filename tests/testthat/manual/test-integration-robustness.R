robustness_fitters <- list(
  gipslda = gipsDA::gipslda,
  gipsqda = gipsDA::gipsqda,
  gipsmultqda = gipsDA::gipsmultqda
)

fit_robust_model <- function(fitter, x, grouping, ...) {
  suppressWarnings(fitter(x, grouping, optimizer = "BF", ...))
}

expect_usable_predictions <- function(fit, newdata, levels) {
  prediction <- predict(fit, newdata)

  expect_length(prediction$class, nrow(newdata))
  expect_setequal(levels(prediction$class), levels)
  expect_false(anyNA(prediction$class))
  expect_valid_posterior(prediction$posterior, nrow(newdata), length(levels))
}

make_robust_binary_data <- function(n_per_class = 10L) {
  i <- seq_len(n_per_class)
  x1 <- c(i / n_per_class, 2 + i / n_per_class)
  x2 <- x1 + rep(c(-1, 1), length.out = 2L * n_per_class) * 1e-5

  list(
    x = cbind(x1 = x1, x2 = x2),
    grouping = factor(rep(c("left", "right"), each = n_per_class))
  )
}

test_that("near-singular and highly correlated variables remain usable", {
  data <- make_robust_binary_data()
  expect_gt(abs(stats::cor(data$x)[1, 2]), 0.999)

  for (fitter in robustness_fitters) {
    fit <- fit_robust_model(fitter, data$x, data$grouping)
    expect_true(all(is.finite(fit$scaling)))
    expect_usable_predictions(fit, data$x, levels(data$grouping))
  }
})

test_that("imbalanced classes and small groups can be classified", {
  grouping <- factor(rep(c("major", "minor", "tiny"), c(16, 6, 4)))
  within_group <- sequence(c(16, 6, 4))
  centers <- c(major = 0, minor = 3, tiny = 6)
  x <- cbind(
    x1 = unname(centers[as.character(grouping)]) + within_group / 20,
    x2 = unname(centers[as.character(grouping)]) -
      within_group / 25 + rep(c(-0.03, 0.03), 13)
  )

  for (fitter in robustness_fitters) {
    fit <- fit_robust_model(fitter, x, grouping)
    expect_equal(unname(fit$counts), c(16, 6, 4))
    expect_equal(fit$prior, prop.table(table(grouping)), ignore_attr = TRUE)
    expect_usable_predictions(fit, x, levels(grouping))
  }
})

test_that("large variable rescaling produces finite usable fits", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 10)
  multipliers <- c(1000, 0.01)
  scaled_x <- sweep(fixture$x, 2, multipliers, `*`)

  for (fitter in robustness_fitters) {
    original <- fit_robust_model(fitter, fixture$x, fixture$grouping)
    rescaled <- fit_robust_model(fitter, scaled_x, fixture$grouping)

    expect_equal(
      unname(rescaled$means),
      unname(sweep(original$means, 2, multipliers, `*`)),
      tolerance = 1e-10
    )
    expect_true(all(is.finite(rescaled$scaling)))
    expect_usable_predictions(rescaled, scaled_x, levels(fixture$grouping))
  }
})

test_that("isolated outliers do not produce invalid predictions", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 10)
  contaminated <- fixture$x
  contaminated[1, ] <- contaminated[1, ] + c(25, -20)

  for (fitter in robustness_fitters) {
    fit <- fit_robust_model(fitter, contaminated, fixture$grouping)
    prediction <- predict(fit, contaminated)

    expect_usable_predictions(fit, contaminated, levels(fixture$grouping))
    expect_gt(mean(prediction$class == fixture$grouping), 1 / 3)
  }
})

test_that("missing training data is rejected or explicitly omitted", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  incomplete <- fixture$x
  incomplete[c(2, 10, 18), 1] <- NA_real_

  for (fitter in robustness_fitters) {
    expect_error(
      fitter(incomplete, fixture$grouping, optimizer = "BF"),
      "infinite, NA or NaN"
    )

    fit <- suppressWarnings(fitter(
      incomplete,
      fixture$grouping,
      na.action = stats::na.omit,
      optimizer = "BF"
    ))
    expect_equal(fit$N, nrow(incomplete) - 3L)
    expect_equal(unname(fit$counts), rep(7, 3))
    expect_usable_predictions(
      fit,
      fixture$x[c(1, 9, 17), , drop = FALSE],
      levels(fixture$grouping)
    )
  }
})
