gaussian_sample <- function(n, mean, covariance) {
  stopifnot(
    length(mean) == nrow(covariance),
    nrow(covariance) == ncol(covariance)
  )
  z <- matrix(stats::rnorm(n * length(mean)), nrow = n)
  sweep(z %*% chol(covariance), 2L, mean, FUN = "+")
}

synthetic_data <- function(n_per_class, means, covariances) {
  groups <- names(means)
  if (length(covariances) == 1L) {
    covariances <- rep(covariances, length(groups))
  }

  x <- do.call(rbind, lapply(seq_along(groups), function(i) {
    gaussian_sample(n_per_class, means[[i]], covariances[[i]])
  }))
  colnames(x) <- paste0("x", seq_len(ncol(x)))

  list(
    x = x,
    grouping = factor(rep(groups, each = n_per_class), levels = groups)
  )
}

expect_synthetic_prediction <- function(prediction, n, groups) {
  expect_length(prediction$class, n)
  expect_equal(levels(prediction$class), groups)
  expect_equal(dim(prediction$posterior), c(n, length(groups)))
  expect_true(all(is.finite(prediction$posterior)))
  expect_true(all(prediction$posterior >= 0 & prediction$posterior <= 1))
  expect_equal(rowSums(prediction$posterior), rep(1, n), tolerance = 1e-10)
}

test_that("shared-covariance synthetic Gaussian data are classified by LDA", {
  set.seed(4101)
  means <- list(
    west = c(-3, 0),
    north = c(0, 3),
    east = c(3, 0)
  )
  covariance <- matrix(c(1, 0.35, 0.35, 1), nrow = 2)
  train <- synthetic_data(45, means, list(covariance))
  test <- synthetic_data(120, means, list(covariance))

  fit <- gipslda(
    train$x,
    train$grouping,
    prior = rep(1 / 3, 3),
    optimizer = "BF"
  )
  prediction <- predict(fit, test$x, method = "plug-in")
  accuracy <- mean(prediction$class == test$grouping)

  expect_synthetic_prediction(prediction, nrow(test$x), names(means))
  expect_gt(accuracy, 0.88)

  near_means <- do.call(rbind, means)
  near_means <- near_means + matrix(c(0.08, -0.05), nrow = 3, ncol = 2, byrow = TRUE)
  expect_equal(
    as.character(predict(fit, near_means)$class),
    names(means)
  )
})

test_that("distinct-covariance Gaussian data are classified by both QDA fits", {
  set.seed(4102)
  means <- list(left = c(-2.6, 0), right = c(2.6, 0))
  covariances <- list(
    left = matrix(c(1.8, 0.55, 0.55, 0.45), nrow = 2),
    right = matrix(c(0.45, -0.25, -0.25, 1.7), nrow = 2)
  )
  train <- synthetic_data(60, means, covariances)
  test <- synthetic_data(180, means, covariances)

  for (fit_function in list(gipsqda, gipsmultqda)) {
    fit <- suppressWarnings(fit_function(
      train$x,
      train$grouping,
      prior = c(0.5, 0.5),
      optimizer = "BF"
    ))
    prediction <- predict(fit, test$x)

    expect_synthetic_prediction(prediction, nrow(test$x), names(means))
    expect_gt(mean(prediction$class == test$grouping), 0.88)
    expect_equal(
      as.character(predict(fit, do.call(rbind, means))$class),
      names(means)
    )
  }
})

test_that("exchangeable covariance gives stable symmetric classifications", {
  set.seed(4103)
  means <- list(low = c(-2.2, -2.2), high = c(2.2, 2.2))
  exchangeable_covariance <- matrix(c(1.1, 0.65, 0.65, 1.1), nrow = 2)
  train <- synthetic_data(70, means, list(exchangeable_covariance))
  test <- synthetic_data(160, means, list(exchangeable_covariance))
  probes <- rbind(
    low = c(-2, -2.15),
    high = c(2, 2.15),
    low_swapped = c(-2.15, -2),
    high_swapped = c(2.15, 2)
  )

  fits <- list(
    gipslda(train$x, train$grouping, optimizer = "BF"),
    suppressWarnings(gipsqda(train$x, train$grouping, optimizer = "BF")),
    suppressWarnings(gipsmultqda(train$x, train$grouping, optimizer = "BF"))
  )

  for (fit in fits) {
    prediction <- predict(fit, test$x)
    probe_prediction <- predict(fit, probes)

    expect_gt(mean(prediction$class == test$grouping), 0.95)
    expect_equal(
      as.character(probe_prediction$class),
      c("low", "high", "low", "high")
    )
    expect_lt(
      max(abs(probe_prediction$posterior[1, ] -
        probe_prediction$posterior[3, ])),
      0.08
    )
    expect_lt(
      max(abs(probe_prediction$posterior[2, ] -
        probe_prediction$posterior[4, ])),
      0.08
    )
  }
})

test_that("all prediction modes return normalized synthetic posteriors", {
  set.seed(4104)
  means <- list(a = c(-2, 0), b = c(2, 0))
  covariances <- list(
    a = matrix(c(1, 0.2, 0.2, 0.7), nrow = 2),
    b = matrix(c(0.7, -0.15, -0.15, 1.2), nrow = 2)
  )
  train <- synthetic_data(36, means, covariances)
  newdata <- rbind(c(-1.8, 0.1), c(1.9, -0.1), c(0, 0))

  lda_fit <- gipslda(train$x, train$grouping, optimizer = "BF")
  for (method in c("plug-in", "predictive", "debiased")) {
    prediction <- predict(lda_fit, newdata, method = method)
    expect_synthetic_prediction(prediction, nrow(newdata), names(means))
  }

  for (fit_function in list(gipsqda, gipsmultqda)) {
    fit <- suppressWarnings(fit_function(
      train$x,
      train$grouping,
      optimizer = "BF"
    ))
    for (method in c("plug-in", "predictive", "debiased")) {
      prediction <- predict(fit, newdata, method = method)
      expect_synthetic_prediction(prediction, nrow(newdata), names(means))
    }

    loo <- predict(fit, method = "looCV")
    expect_synthetic_prediction(loo, nrow(train$x), names(means))
    expect_gt(mean(loo$class == train$grouping), 0.80)
  }
})

test_that("priors move decisions for points near the class boundary", {
  set.seed(4105)
  means <- list(a = c(-1.7, 0), b = c(1.7, 0))
  covariance <- matrix(c(0.9, 0.25, 0.25, 0.9), nrow = 2)
  train <- synthetic_data(80, means, list(covariance))
  boundary <- matrix(c(0, 0), nrow = 1)

  fits <- list(
    gipslda(train$x, train$grouping, optimizer = "BF"),
    suppressWarnings(gipsqda(train$x, train$grouping, optimizer = "BF")),
    suppressWarnings(gipsmultqda(train$x, train$grouping, optimizer = "BF"))
  )

  for (fit in fits) {
    favor_a <- predict(fit, boundary, prior = c(0.9, 0.1))
    favor_b <- predict(fit, boundary, prior = c(0.1, 0.9))
    balanced <- predict(fit, boundary, prior = c(0.5, 0.5))

    expect_equal(as.character(favor_a$class), "a")
    expect_equal(as.character(favor_b$class), "b")
    expect_gt(favor_a$posterior[1, "a"], balanced$posterior[1, "a"])
    expect_gt(favor_b$posterior[1, "b"], balanced$posterior[1, "b"])
    expect_lt(max(balanced$posterior), 0.75)
  }
})
