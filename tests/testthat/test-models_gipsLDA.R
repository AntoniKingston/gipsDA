# ==============================================================================
# Global Setup for this file
# ==============================================================================

# We use the standard Iris dataset for testing
data(iris)

# Split into Matrix (X) and Grouping Factor (Y) for testing matrix interfaces
X <- as.matrix(iris[, 1:4])
Y <- iris$Species

# ==============================================================================
# 1. Testing Model Fitting Interfaces (S3 Dispatch)
# ==============================================================================

test_that("gipslda formula interface works correctly", {
  # Action
  fit_formula <- gipslda(Species ~ ., data = iris)

  # Check class and type
  expect_s3_class(fit_formula, "gipslda")
  expect_type(fit_formula, "list")

  # Check key components
  required_components <- c("prior", "means", "scaling", "counts")
  expect_true(all(required_components %in% names(fit_formula)))

  # Check dimensions (3 classes x 4 variables)
  expect_equal(dim(fit_formula$means), c(3, 4))
})

test_that("gipslda matrix interface works and matches formula", {
  # Setup formula fit again for comparison
  fit_formula <- gipslda(Species ~ ., data = iris)

  # Action
  fit_matrix <- gipslda(x = X, grouping = Y)

  expect_s3_class(fit_matrix, "gipslda")

  # Consistency check: priors and counts should be identical
  expect_equal(fit_matrix$prior, fit_formula$prior)
  expect_equal(fit_matrix$counts, fit_formula$counts)
})

test_that("gipslda data.frame interface works", {
  fit_df <- gipslda(x = iris[, 1:4], grouping = Y)
  expect_s3_class(fit_df, "gipslda")
})

test_that("gipslda handles custom arguments", {
  # Verify no error on custom args
  fit_opts <- gipslda(X, Y, MAP = FALSE, weighted_avg = TRUE)
  expect_s3_class(fit_opts, "gipslda")
})

# ==============================================================================
# 2. Testing the 'print' method
# ==============================================================================

test_that("print.gipslda works", {
  fit_formula <- gipslda(Species ~ ., data = iris)

  # expect_output checks if printing produces any text to console
  # You can also add a regexp argument to check for specific words
  expect_output(print(fit_formula))
})

# ==============================================================================
# 3. Testing the 'predict' method
# ==============================================================================

test_that("predict.gipslda works correctly", {
  # Setup
  fit_formula <- gipslda(Species ~ ., data = iris)

  # Action: Basic prediction
  pred <- predict(fit_formula, newdata = iris)

  expect_type(pred, "list")

  # Check components
  expect_equal(names(pred), c("class", "posterior", "x"))

  # Check dimensions
  expect_length(pred$class, nrow(iris))
  expect_equal(nrow(pred$posterior), nrow(iris))
  expect_equal(ncol(pred$posterior), 3) # 3 classes

  # Check row sums (probabilities sum to 1)
  row_sums <- rowSums(pred$posterior)
  # expect_equal handles tolerance automatically for doubles
  expect_equal(row_sums, rep(1, nrow(iris)), tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("predict.gipslda works on new subset", {
  # Setup
  fit_formula <- gipslda(Species ~ ., data = iris)
  new_data <- iris[1:5, ]

  # Action
  pred_small <- predict(fit_formula, newdata = new_data)

  # Assert
  expect_length(pred_small$class, 5)
})