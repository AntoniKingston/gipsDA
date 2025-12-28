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
  expect_equal(unname(row_sums), rep(1, nrow(iris)), tolerance = 1e-6,
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

# ==============================================================================
# 4. Input Validation & Error Handling
# ==============================================================================

test_that("gipslda correctly validates bad inputs", {
  # Infinite or NA values
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(gipslda(X_na, Y), "infinite, NA or NaN values in 'x'")

  # Dimension mismatch
  expect_error(gipslda(X[1:10, ], Y), "nrow\\(x\\) and length\\(grouping\\) are different")

  # 'x' is not a matrix
  expect_error(gipslda.default(1:10, rep(1, 10)), "'x' is not a matrix")

  # Invalid priors
  expect_error(gipslda(X, Y, prior = c(0.5, 0.5, 0.5)), "invalid 'prior'")
  expect_error(gipslda(X, Y, prior = c(1.2, -0.1, -0.1)), "invalid 'prior'")
  expect_error(gipslda(X, Y, prior = c(0.5, 0.5)), "'prior' is of incorrect length")
})

test_that("gipslda handles empty groups and constant variables", {
  # Empty group warning
  Y_empty <- factor(Y, levels = c("setosa", "versicolor", "virginica", "ghost_group"))
  expect_warning(gipslda(X, Y_empty), "group ghost_group is empty")

  # Constant variable within groups
  X_const <- X
  X_const[, 1] <- 1
  expect_error(gipslda(X_const, Y), "variable .* appears to be constant within groups")
})

test_that("gipslda handles rank deficiency", {
  # Case where variables are numerically constant (0 variance)
  X_zero <- matrix(0, nrow = 150, ncol = 4)

  # Updated regex to handle both singular ("variable ... appears")
  # and plural ("variables ... appear") error messages.
  expect_error(gipslda(X_zero, Y), "variables? .* appear(s)? to be constant")
})

# ==============================================================================
# 5. Matrix Interface: Subset & NA Action
# ==============================================================================

test_that("gipslda.matrix handles subset and na.action properly", {
  # 1. Test subset
  # Selecting 10 from each species.
  subset_idx <- c(1:10, 51:60, 101:110)

  fit_sub <- gipslda(X, Y, subset = subset_idx)
  expect_equal(fit_sub$N, 30)

  # 2. Test na.action
  X_dirty <- X
  X_dirty[1, 1] <- NA

  # Workaround for missing row.names in source code
  safe_na_omit <- function(obj) {
    attr(obj, "row.names") <- seq_along(obj$g)
    stats::na.omit(obj)
  }

  fit_na <- gipslda(X_dirty, Y, na.action = safe_na_omit)
  expect_equal(fit_na$N, 149)
})

# ==============================================================================
# 6. Optimizer Selection & Special Arguments
# ==============================================================================

test_that("optimizer logic and special args work", {
  # Default behavior (BF for p < 10)
  fit_def <- gipslda(X, Y, optimizer = NULL)
  expect_s3_class(fit_def, "gipslda")

  # Warning when max_iter is missing for MH
  expect_warning(
    gipslda(X, Y, optimizer = "MH"),
    "MH optimizer set but 'max_iter' argument is unspecified"
  )

  # Weighted average path (weighted_avg = TRUE)
  # This triggers the "if (weighted_avg)" block in gipslda.default
  fit_weighted <- gipslda(X, Y, weighted_avg = TRUE)
  expect_s3_class(fit_weighted, "gipslda")
})

# ==============================================================================
# 7. Prediction Methods & Edge Cases
# ==============================================================================

test_that("predict supports different methods", {
  fit <- gipslda(X, Y)

  # Plug-in method (default check done in standard test, here specific)
  pred_plug <- predict(fit, newdata = X, method = "plug-in")
  expect_equal(nrow(pred_plug$posterior), 150)

  # Debiased method
  pred_deb <- predict(fit, newdata = X, method = "debiased")
  expect_equal(nrow(pred_deb$posterior), 150)

  # Predictive method (Bayesian) - covers the 'else' block in predict logic
  pred_pred <- predict(fit, newdata = X, method = "predictive")
  expect_equal(nrow(pred_pred$posterior), 150)
  expect_false(any(is.na(pred_pred$posterior)))
})

test_that("predict catches invalid arguments", {
  fit <- gipslda(X, Y)

  # Dimension mismatch in newdata
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

test_that("predict reconstructs data when newdata is missing (Subset logic)", {
  subset_idx <- c(1:10, 51:60, 101:110)

  fit_sub <- gipslda(X, Y, subset = subset_idx)

  # Call predict without newdata -> triggers reconstruction
  pred <- predict(fit_sub)

  expect_equal(nrow(pred$posterior), 30)
  expect_length(pred$class, 30)
})

# ==============================================================================
# 8. Plotting & Auxiliary Methods (Smoke Tests)
# ==============================================================================

test_that("coef.gipslda returns scaling matrix", {
  fit <- gipslda(X, Y)
  expect_equal(coef(fit), fit$scaling)
})

test_that("plot.gipslda runs without error", {
  # The plot method uses MASS::eqscplot. We must ensure MASS is available.
  skip_if_not_installed("MASS")
  # We attach MASS just for this test so the function eqscplot is found
  library(MASS)

  fit <- gipslda(X, Y)

  # Use temp file to avoid opening graphics window
  pdf(NULL)
  on.exit(dev.off())

  # Standard plot
  expect_error(plot(fit), NA)

  # Plot with abbreviation
  expect_error(plot(fit, abbrev = TRUE), NA)

  # Plot with subset
  fit_sub <- gipslda(X, Y, subset = c(1:10, 51:60, 101:110))
  expect_error(plot(fit_sub), NA)
})

test_that("pairs.gipslda runs without error", {
  fit <- gipslda(X, Y)

  pdf(NULL)
  on.exit(dev.off())

  # Standard pairs
  expect_error(pairs(fit), NA)

  # Trellis type (requires lattice)
  if (requireNamespace("lattice", quietly = TRUE)) {
    expect_error(pairs(fit, type = "trellis"), NA)
  }
})

test_that("model.frame.gipslda works", {
  fit <- gipslda(Species ~ ., data = iris)
  mf <- model.frame(fit)

  expect_s3_class(mf, "data.frame")
  expect_equal(nrow(mf), 150)
})
