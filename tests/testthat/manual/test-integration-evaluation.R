stratified_iris_split <- function(seed, train_per_class = 30L) {
  set.seed(seed)
  by_class <- split(seq_len(nrow(iris)), iris$Species)
  train <- unlist(lapply(
    by_class,
    sample,
    size = train_per_class,
    replace = FALSE
  ), use.names = FALSE)

  list(train = sort(train), test = setdiff(seq_len(nrow(iris)), train))
}

iris_accuracy <- function(predicted, truth) {
  mean(predicted == truth)
}

iris_confusion <- function(predicted, truth) {
  table(
    predicted = factor(predicted, levels = levels(iris$Species)),
    truth = factor(truth, levels = levels(iris$Species))
  )
}

fit_iris_models <- function(train_data) {
  list(
    gipslda = gipsDA::gipslda(
      Species ~ .,
      train_data,
      optimizer = "BF"
    ),
    gipsqda = suppressWarnings(gipsDA::gipsqda(
      Species ~ .,
      train_data,
      optimizer = "BF"
    )),
    gipsmultqda = suppressWarnings(gipsDA::gipsmultqda(
      Species ~ .,
      train_data,
      optimizer = "BF"
    )),
    MASSlda = MASS::lda(Species ~ ., train_data),
    MASSqda = MASS::qda(Species ~ ., train_data)
  )
}

test_that("stratified iris evaluation beats the majority baseline", {
  split <- stratified_iris_split(seed = 20260724)
  train_data <- iris[split$train, , drop = FALSE]
  test_data <- iris[split$test, , drop = FALSE]
  expect_equal(as.vector(table(train_data$Species)), rep(30L, 3))
  expect_equal(as.vector(table(test_data$Species)), rep(20L, 3))

  models <- fit_iris_models(train_data)
  predictions <- lapply(models, function(model) {
    predict(model, test_data)$class
  })
  confusion_matrices <- lapply(
    predictions,
    iris_confusion,
    truth = test_data$Species
  )
  accuracies <- vapply(
    predictions,
    iris_accuracy,
    truth = test_data$Species,
    numeric(1)
  )
  majority_baseline <- max(prop.table(table(test_data$Species)))

  for (confusion in confusion_matrices) {
    expect_equal(dim(confusion), c(3L, 3L))
    expect_equal(sum(confusion), nrow(test_data))
    expect_identical(rownames(confusion), levels(iris$Species))
    expect_identical(colnames(confusion), levels(iris$Species))
  }
  expect_true(all(accuracies > majority_baseline))
  expect_gte(unname(accuracies["gipslda"]), unname(accuracies["MASSlda"]) - 0.10)
  expect_gte(unname(accuracies["gipsqda"]), unname(accuracies["MASSqda"]) - 0.10)
  expect_gte(
    unname(accuracies["gipsmultqda"]),
    unname(accuracies["MASSqda"]) - 0.10
  )
})

test_that("stratified iris accuracy is stable across repeated seeds", {
  testthat::skip_on_cran()
  seeds <- c(104L, 305L, 709L)
  model_names <- c("gipslda", "gipsqda", "gipsmultqda")
  accuracies <- matrix(
    NA_real_,
    nrow = length(seeds),
    ncol = length(model_names),
    dimnames = list(seeds, model_names)
  )

  for (i in seq_along(seeds)) {
    split <- stratified_iris_split(seeds[i])
    train_data <- iris[split$train, , drop = FALSE]
    test_data <- iris[split$test, , drop = FALSE]
    models <- fit_iris_models(train_data)[model_names]

    accuracies[i, ] <- vapply(models, function(model) {
      iris_accuracy(predict(model, test_data)$class, test_data$Species)
    }, numeric(1))
  }

  expect_true(all(accuracies > 1 / 3))
  expect_true(all(colMeans(accuracies) >= 0.80))
  expect_true(all(apply(accuracies, 2, stats::sd) <= 0.10))
})
