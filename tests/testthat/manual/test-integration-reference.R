make_integration_data <- function(n_per_class, seed) {
  set.seed(seed)
  means <- rbind(c(-2, 0), c(0, 2), c(2, -1))
  scales <- rbind(c(0.6, 2.5), c(0.8, 2.8), c(0.7, 2.2))
  groups <- c("a", "b", "c")

  x <- do.call(rbind, lapply(seq_along(groups), function(i) {
    noise <- matrix(stats::rnorm(2L * n_per_class), ncol = 2L)
    sweep(noise, 2L, scales[i, ], `*`) +
      matrix(means[i, ], n_per_class, 2L, byrow = TRUE)
  }))
  colnames(x) <- c("x1", "x2")
  grouping <- factor(rep(groups, each = n_per_class), levels = groups)

  list(
    x = x,
    grouping = grouping,
    data = data.frame(class = grouping, x, check.names = FALSE)
  )
}

normalize_discriminants <- function(scaling) {
  sweep(scaling, 2L, sqrt(colSums(scaling^2)), `/`)
}

qda_precision <- function(fit, group) {
  tcrossprod(fit$scaling[, , group])
}

test_that("gipslda agrees with MASS when covariance projection is unchanged", {
  train <- make_integration_data(15L, 20260724)
  test <- make_integration_data(5L, 20260725)

  gips_fit <- gipslda(
    train$x,
    train$grouping,
    optimizer = "BF"
  )
  mass_fit <- MASS::lda(train$x, train$grouping)
  gips_prediction <- predict(gips_fit, test$x)
  mass_prediction <- predict(mass_fit, test$x)

  expect_equal(gips_fit$means, mass_fit$means)
  expect_equal(gips_fit$prior, mass_fit$prior)
  expect_equal(
    abs(normalize_discriminants(gips_fit$scaling)),
    abs(normalize_discriminants(mass_fit$scaling)),
    tolerance = 1e-10
  )
  expect_identical(gips_prediction$class, mass_prediction$class)
  expect_equal(
    unname(gips_prediction$posterior),
    unname(mass_prediction$posterior),
    tolerance = 0.01
  )
})

test_that("gipsqda agrees with MASS when covariance projection is unchanged", {
  train <- make_integration_data(15L, 20260724)
  test <- make_integration_data(5L, 20260725)

  gips_fit <- gipsqda(
    train$x,
    train$grouping,
    optimizer = "BF"
  )
  mass_fit <- MASS::qda(train$x, train$grouping)
  gips_prediction <- predict(gips_fit, test$x)
  mass_prediction <- predict(mass_fit, test$x)

  expect_equal(gips_fit$means, mass_fit$means)
  expect_equal(gips_fit$prior, mass_fit$prior)
  expect_equal(gips_fit$ldet, mass_fit$ldet, tolerance = 1e-10)
  for (group in seq_along(gips_fit$prior)) {
    # The scaling factors may differ by column signs; their precision matrices
    # are the sign-invariant quantities used by the discriminant rule.
    expect_equal(
      qda_precision(gips_fit, group),
      qda_precision(mass_fit, group),
      tolerance = 1e-10
    )
  }
  expect_identical(gips_prediction$class, mass_prediction$class)
  expect_equal(
    unname(gips_prediction$posterior),
    unname(mass_prediction$posterior),
    tolerance = 1e-10
  )
})

test_that("all interfaces give consistent gipsDA predictions", {
  train <- make_integration_data(10L, 20260726)
  test <- make_integration_data(3L, 20260727)
  fitters <- list(
    gipslda = gipsDA::gipslda,
    gipsqda = gipsDA::gipsqda,
    gipsmultqda = gipsDA::gipsmultqda
  )

  for (name in names(fitters)) {
    fit <- fitters[[name]]
    matrix_fit <- fit(train$x, train$grouping, optimizer = "BF")
    frame_fit <- fit(
      as.data.frame(train$x),
      train$grouping,
      optimizer = "BF"
    )
    formula_fit <- fit(
      class ~ x1 + x2,
      data = train$data,
      optimizer = "BF"
    )

    matrix_prediction <- predict(matrix_fit, test$x)
    frame_prediction <- predict(
      frame_fit,
      as.data.frame(test$x)
    )
    formula_prediction <- predict(formula_fit, test$data)

    expect_equal(unname(frame_fit$means), unname(matrix_fit$means), info = name)
    expect_equal(unname(formula_fit$means), unname(matrix_fit$means), info = name)
    expect_identical(frame_prediction$class, matrix_prediction$class, info = name)
    expect_identical(formula_prediction$class, matrix_prediction$class, info = name)
    expect_equal(
      unname(frame_prediction$posterior),
      unname(matrix_prediction$posterior),
      tolerance = 1e-10,
      info = name
    )
    expect_equal(
      unname(formula_prediction$posterior),
      unname(matrix_prediction$posterior),
      tolerance = 1e-10,
      info = name
    )
  }
})
