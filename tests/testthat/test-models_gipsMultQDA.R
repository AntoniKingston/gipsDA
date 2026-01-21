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
  expect_equal(unname(row_sums), rep(1, nrow(iris)),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )
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

# ==============================================================================
# 4. Testing Input Validation and Error Handling (Covering gipsmultqda.default)
# ==============================================================================

test_that("gipsmultqda throws errors for invalid inputs", {
  # 1. Test infinite/NA values check
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(gipsmultqda(X_na, Y), "infinite, NA or NaN values in 'x'")

  X_inf <- X
  X_inf[1, 1] <- Inf
  expect_error(gipsmultqda(X_inf, Y), "infinite, NA or NaN values in 'x'")

  # 2. Test dimension mismatch
  expect_error(
    gipsmultqda(X[1:10, ], Y),
    "nrow\\(x\\) and length\\(grouping\\) are different"
  )

  # 3. Test 'x' is not a matrix (passed as vector without dim)
  expect_error(gipsmultqda.default(1:10, rep(1, 10)), "'x' is not a matrix")

  # 4. Test small group size
  # Create a dataset where one group has fewer samples than variables (p=4)
  # Iris has 50 per group. Let's make a tiny subset.
  # We need counts < p+1 (so < 5).
  small_X <- X[c(1:4, 51:54, 101:104), ]
  small_Y <- factor(c(rep("s", 4), rep("ve", 4), rep("vi", 4)))
  expect_no_error(gipsmultqda(small_X, small_Y))
})

test_that("gipsmultqda validates priors", {
  # 1. Priors don't sum to 1
  expect_error(gipsmultqda(X, Y, prior = c(0.5, 0.5, 0.5)), "invalid 'prior'")

  # 2. Negative priors
  expect_error(gipsmultqda(X, Y, prior = c(1.2, -0.1, -0.1)), "invalid 'prior'")

  # 3. Incorrect length
  expect_error(
    gipsmultqda(X, Y, prior = c(0.5, 0.5)),
    "'prior' is of incorrect length"
  )
})

# ==============================================================================
# 5. Testing Matrix Interface Specifics (Subset & NA Action)
# ==============================================================================

test_that("gipsmultqda.matrix handles subset and na.action", {
  subset_idx <- c(1:10, 51:60, 101:110)
  fit_sub <- gipsmultqda(X, Y, subset = subset_idx)
  expect_equal(fit_sub$N, 30)

  X_dirty <- X
  X_dirty[1, 1] <- NA
  safe_na_omit <- function(object) {
    attr(object, "row.names") <- seq_along(object$g)
    stats::na.omit(object)
  }
  fit_na <- gipsmultqda(X_dirty, Y, na.action = safe_na_omit)

  expect_equal(fit_na$N, 149)
})

# ==============================================================================
# 6. Testing Optimizer Logic
# ==============================================================================

test_that("optimizer selection logic works", {
  # 1. Default BF for p < 10 (Iris has p=4)
  fit_bf <- gipsmultqda(X, Y, optimizer = NULL)
  # Check internals if available, or just ensure it runs without warning
  expect_s3_class(fit_bf, "gipsmultqda")

  # 2. Manual MH optimizer selection
  fit_mh <- gipsmultqda(X, Y, optimizer = "MH", max_iter = 10)
  expect_s3_class(fit_mh, "gipsmultqda")

  # 3. Warning when MH is chosen but max_iter is missing
  expect_warning(
    gipsmultqda(X, Y, optimizer = "MH"),
    "MH optimizer set but 'max_iter' argument is unspecified"
  )
})

# ==============================================================================
# 7. Testing Prediction Methods & Edge Cases
# ==============================================================================

test_that("predict handles different methods", {
  fit <- gipsmultqda(X, Y)

  # 1. Predictive method (corresponds to 'else' block in predict function)
  pred_pred <- predict(fit, newdata = X, method = "predictive")
  expect_equal(nrow(pred_pred$posterior), 150)

  # 2. Debiased method
  pred_deb <- predict(fit, newdata = X, method = "debiased")
  expect_equal(nrow(pred_deb$posterior), 150)

  # 3. looCV (Leave-One-Out Cross Validation)
  # Note: looCV works on the training data, usually without 'newdata'
  pred_loo <- predict(fit, method = "looCV")
  expect_equal(nrow(pred_loo$posterior), 150)
})

test_that("predict throws errors for invalid scenarios", {
  fit <- gipsmultqda(X, Y)

  # 1. looCV with newdata provided (should fail)
  expect_error(
    predict(fit, newdata = X, method = "looCV"),
    "cannot have leave-one-out CV with 'newdata'"
  )

  # 2. Wrong dimensions in newdata (e.g., only 2 columns instead of 4)
  expect_error(predict(fit, newdata = X[, 1:2]), "wrong number of variables")

  # 3. Invalid prior in predict
  expect_error(
    predict(fit, newdata = X, prior = c(0.5, 0.5)),
    "'prior' is of incorrect length"
  )

  # 4. Variable name mismatch warning
  X_bad_names <- X
  colnames(X_bad_names) <- c("A", "B", "C", "D")
  expect_warning(
    predict(fit, newdata = X_bad_names),
    "variable names in 'newdata' do not match those in 'object'"
  )
})

test_that("predict works without newdata (re-substitution)", {
  # This tests the logic inside predict where it tries to reconstruct data
  # from the call object when newdata is missing.

  # Case A: Formula interface
  fit_f <- gipsmultqda(Species ~ ., data = iris)
  pred_f <- predict(fit_f) # No newdata
  expect_equal(nrow(pred_f$posterior), 150)

  # Case B: Matrix interface
  fit_m <- gipsmultqda(X, Y)
  pred_m <- predict(fit_m) # No newdata
  expect_equal(nrow(pred_m$posterior), 150)
})

test_that("predict reconstructs data when using subset (covers eval.parent)", {
  # 1. Defining subset indices
  subset_idx <- c(1:10, 51:60, 101:110)

  # 2. Creating matrix fit with subset
  fit_subset <- gipsmultqda(X, Y, subset = subset_idx)

  # 3. Without newdata, predict should reconstruct data using subset
  pred <- predict(fit_subset)

  expect_equal(nrow(pred$posterior), 30)
  expect_length(pred$class, 30)
})
