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

test_that("gipsqda formula interface works correctly", {
  fit_formula <- gipsqda(Species ~ ., data = iris)

  # Check class and type
  expect_s3_class(fit_formula, "gipsqda")
  expect_type(fit_formula, "list")

  # Check if the result contains key components
  required_components <- c("prior", "means", "scaling", "counts")
  expect_true(all(required_components %in% names(fit_formula)))

  # Check dimensions of the means (3 classes x 4 variables)
  expect_equal(dim(fit_formula$means), c(3, 4))
})

test_that("gipsqda matrix interface works and matches formula", {
  # Re-fit formula for comparison
  fit_formula <- gipsqda(Species ~ ., data = iris)

  # Action
  fit_matrix <- gipsqda(x = X, grouping = Y)

  expect_s3_class(fit_matrix, "gipsqda")

  # The results between formula and matrix interface should be consistent
  expect_equal(fit_matrix$prior, fit_formula$prior)
  expect_equal(fit_matrix$counts, fit_formula$counts)
})

test_that("gipsqda data.frame interface works", {
  fit_df <- gipsqda(x = iris[, 1:4], grouping = Y)
  expect_s3_class(fit_df, "gipsqda")
})

test_that("gipsqda handles custom arguments", {
  fit_opts <- gipsqda(X, Y, MAP = FALSE, weighted_avg = TRUE)
  expect_s3_class(fit_opts, "gipsqda")
})

# ==============================================================================
# 2. Testing the 'print' method
# ==============================================================================

test_that("print.gipsqda works", {
  fit_formula <- gipsqda(Species ~ ., data = iris)

  # expect_output asserts that something is printed to the console
  expect_output(print(fit_formula))
})

# ==============================================================================
# 3. Testing the 'predict' method
# ==============================================================================

test_that("predict.gipsqda works correctly", {
  # Setup
  fit_formula <- gipsqda(Species ~ ., data = iris)

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

test_that("predict.gipsqda works on new subset", {
  # Setup
  fit_formula <- gipsqda(Species ~ ., data = iris)
  new_data <- iris[1:5, ]

  # Action
  pred_small <- predict(fit_formula, newdata = new_data)

  # Assert
  expect_length(pred_small$class, 5)
})

# ==============================================================================
# 4. Input Validation & Error Handling
# ==============================================================================

test_that("gipsqda correctly validates bad inputs", {
  # Infinite or NA values
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(gipsqda(X_na, Y), "infinite, NA or NaN values in 'x'")

  X_inf <- X
  X_inf[1, 1] <- Inf
  expect_error(gipsqda(X_inf, Y), "infinite, NA or NaN values in 'x'")

  # Dimension mismatch
  expect_error(gipsqda(X[1:10, ], Y), "nrow\\(x\\) and length\\(grouping\\) are different")

  # 'x' is not a matrix (vector input)
  expect_error(gipsqda.default(1:10, rep(1, 10)), "'x' is not a matrix")

  # Small group size (group count < p+1)
  # Create subset where groups have < 5 observations (since p=4)
  small_X <- X[c(1:4, 51:54, 101:104), ]
  small_Y <- factor(c(rep("s", 4), rep("ve", 4), rep("vi", 4)))
  expect_no_error(gipsqda(small_X, small_Y))
})

test_that("gipsqda validates prior probabilities", {
  # Sum != 1
  expect_error(gipsqda(X, Y, prior = c(0.5, 0.5, 0.5)), "invalid 'prior'")

  # Negative values
  expect_error(gipsqda(X, Y, prior = c(1.2, -0.1, -0.1)), "invalid 'prior'")

  # Wrong length
  expect_error(gipsqda(X, Y, prior = c(0.5, 0.5)), "'prior' is of incorrect length")
})

test_that("gipsqda handles missing dimnames", {
  X_no_names <- X
  dimnames(X_no_names) <- NULL

  fit <- gipsqda(X_no_names, Y)

  expect_s3_class(fit, "gipsqda")
  # Should generate default colnames "1", "2", "3", "4"
  expect_equal(colnames(fit$means), as.character(1:4))
})

# ==============================================================================
# 5. Matrix Interface: Subset & NA Action
# ==============================================================================

test_that("gipsqda.matrix handles subset and na.action properly", {
  # 1. Test subset
  # Must pick samples from ALL groups to avoid "group too small" error
  subset_idx <- c(1:10, 51:60, 101:110)
  fit_sub <- gipsqda(X, Y, subset = subset_idx)

  expect_equal(fit_sub$N, 30) # 3 groups * 10 obs

  # 2. Test na.action
  X_dirty <- X
  X_dirty[1, 1] <- NA

  # Workaround for missing row.names in the source code
  safe_na_omit <- function(obj) {
    attr(obj, "row.names") <- seq_along(obj$g)
    stats::na.omit(obj)
  }

  fit_na <- gipsqda(X_dirty, Y, na.action = safe_na_omit)

  # Expect 1 row removed (150 -> 149)
  expect_equal(fit_na$N, 149)
})

# ==============================================================================
# 6. Optimizer Selection Logic
# ==============================================================================

test_that("optimizer logic and warnings work", {
  # Default behavior
  fit_def <- gipsqda(X, Y, optimizer = NULL)
  expect_s3_class(fit_def, "gipsqda")

  # Warning when max_iter is missing for MH
  expect_warning(
    gipsqda(X, Y, optimizer = "MH"),
    "MH optimizer set but 'max_iter' argument is unspecified"
  )

  # No warning when explicitly set
  expect_no_warning(
    gipsqda(X, Y, optimizer = "MH", max_iter = 10)
  )
})

# ==============================================================================
# 7. Prediction Methods & Edge Cases
# ==============================================================================

test_that("predict supports different methods", {
  fit <- gipsqda(X, Y)

  # Predictive (Bayesian) - covers the 'else' block
  pred_bay <- predict(fit, newdata = X, method = "predictive")
  expect_equal(nrow(pred_bay$posterior), 150)
  expect_false(any(is.na(pred_bay$posterior)))

  # Debiased
  pred_deb <- predict(fit, newdata = X, method = "debiased")
  expect_equal(nrow(pred_deb$posterior), 150)

  # looCV (Leave-One-Out)
  pred_loo <- predict(fit, method = "looCV")
  expect_equal(nrow(pred_loo$posterior), 150)
})

test_that("predict catches invalid arguments", {
  fit <- gipsqda(X, Y)

  # Forbidden combination: looCV + newdata
  expect_error(predict(fit, newdata = X, method = "looCV"),
               "cannot have leave-one-out CV with 'newdata'")

  # Wrong dimensions
  expect_error(predict(fit, newdata = X[, 1:2]), "wrong number of variables")

  # Invalid prior length
  expect_error(predict(fit, newdata = X, prior = c(0.5, 0.5)),
               "'prior' is of incorrect length")

  # Variable name mismatch warning
  X_bad <- X
  colnames(X_bad) <- c("A", "B", "C", "D")
  expect_warning(predict(fit, newdata = X_bad),
                 "variable names in 'newdata' do not match")
})

test_that("predict reconstructs data from subset when newdata is missing", {
  # Covers the complex eval.parent block

  # Use indices from all groups
  subset_idx <- c(1:10, 51:60, 101:110)
  fit_sub <- gipsqda(X, Y, subset = subset_idx)

  # Call predict without newdata -> triggers reconstruction
  pred <- predict(fit_sub)

  expect_equal(nrow(pred$posterior), 30)
  expect_length(pred$class, 30)
})
