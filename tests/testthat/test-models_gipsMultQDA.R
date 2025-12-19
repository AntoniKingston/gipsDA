# ==============================================================================
# Global Setup for this file
# ==============================================================================
data(iris)

# Split into Matrix (X) and Grouping Factor (Y) for testing matrix interfaces
X <- as.matrix(iris[, 1:4])
Y <- iris$Species

# ==============================================================================
# 1. Testing Model Fitting Interfaces (S3 Dispatch)
# ==============================================================================

test_that("gipsmultqda formula interface works correctly", {
  fit_formula <- gipsmultqda(Species ~ ., data = iris)

  # Check class and type
  expect_s3_class(fit_formula, "gipsmultqda")
  expect_type(fit_formula, "list")

  # Check if the result contains key components
  required_components <- c("prior", "means", "scaling", "counts")
  expect_true(all(required_components %in% names(fit_formula)))

  # Check dimensions of the means (3 classes x 4 variables)
  expect_equal(dim(fit_formula$means), c(3, 4))
})

test_that("gipsmultqda matrix interface works and matches formula", {
  # Re-fit formula for comparison
  fit_formula <- gipsmultqda(Species ~ ., data = iris)

  # Action
  fit_matrix <- gipsmultqda(x = X, grouping = Y)

  expect_s3_class(fit_matrix, "gipsmultqda")

  # The results between formula and matrix interface should be consistent
  expect_equal(fit_matrix$prior, fit_formula$prior)
  expect_equal(fit_matrix$counts, fit_formula$counts)
})

test_that("gipsmultqda data.frame interface works", {
  fit_df <- gipsmultqda(x = iris[, 1:4], grouping = Y)
  expect_s3_class(fit_df, "gipsmultqda")
})

test_that("gipsmultqda handles custom arguments", {
  fit_opts <- gipsmultqda(X, Y, MAP = FALSE, weighted_avg = TRUE)
  expect_s3_class(fit_opts, "gipsmultqda")
})

# ==============================================================================
# 2. Testing the 'print' method
# ==============================================================================

test_that("print.gipsmultqda works", {
  fit_formula <- gipsmultqda(Species ~ ., data = iris)

  # expect_output asserts that something is printed to the console
  expect_output(print(fit_formula))
})

# ==============================================================================
# 3. Testing the 'predict' method
# ==============================================================================

test_that("predict.gipsmultqda works correctly", {
  # Setup
  fit_formula <- gipsmultqda(Species ~ ., data = iris)

  # Action
  pred <- predict(fit_formula, newdata = iris)

  expect_type(pred, "list")

  # Check components
  expect_equal(names(pred), c("class", "posterior"))

  # Check dimensions matching the input data
  expect_length(pred$class, nrow(iris))
  expect_equal(nrow(pred$posterior), nrow(iris))
  expect_equal(ncol(pred$posterior), 3) # 3 classes

  # Check if posterior probabilities sum to 1 (row-wise)
  row_sums <- rowSums(pred$posterior)

  # expect_equal with a vector handles tolerance automatically
  expect_equal(unname(row_sums), rep(1, nrow(iris)), tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("predict.gipsmultqda works on new subset", {
  # Setup
  fit_formula <- gipsmultqda(Species ~ ., data = iris)
  new_data <- iris[1:5, ]

  # Action
  pred_small <- predict(fit_formula, newdata = new_data)

  # Assert
  expect_length(pred_small$class, 5)
})