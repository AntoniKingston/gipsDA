test_that("formula fits can be predicted after fitting-only arguments", {
  for (fitter in list(gipslda, gipsqda, gipsmultqda)) {
    fit <- fitter(
      Species ~ Sepal.Length + Sepal.Width,
      data = iris,
      optimizer = "BF",
      show_progress_bar = FALSE,
      store_probabilities = FALSE
    )

    expect_silent(pred <- predict(fit))
    expect_equal(length(pred$class), nrow(iris))
    expect_equal(nrow(pred$posterior), nrow(iris))
  }
})

test_that("matrix fits with na.action can be predicted without newdata", {
  x <- as.matrix(iris[, c("Sepal.Length", "Sepal.Width")])
  grouping <- iris$Species
  x[1, 1] <- NA

  for (fitter in list(gipslda, gipsqda, gipsmultqda)) {
    fit <- fitter(
      x,
      grouping,
      optimizer = "BF",
      na.action = na.omit
    )

    expect_silent(pred <- predict(fit))
    expect_equal(length(pred$class), nrow(iris) - 1L)
    expect_equal(nrow(pred$posterior), nrow(iris) - 1L)
  }
})

test_that("LDA diagnostics work after matrix fit with na.action", {
  x <- as.matrix(iris[, c("Sepal.Length", "Sepal.Width")])
  grouping <- iris$Species
  x[1, 1] <- NA

  fit <- gipslda(
    x,
    grouping,
    optimizer = "BF",
    na.action = na.omit
  )

  expect_silent(plot(fit))
  expect_silent(pairs(fit))
})
