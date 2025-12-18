# ==============================================================================
# Global Setup for this file
# ==============================================================================
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
S2 <- matrix(c(2, -0.2, -0.2, 2), nrow = 2, byrow = TRUE)
Ss <- list(S1, S2)
ns <- c(10L, 15L) # Integers

# ==============================================================================
# 1. Constructor and Validation Tests
# ==============================================================================

test_that("gipsmult constructor creates a valid object", {
  # Action
  g <- gipsmult(Ss, ns)

  # Assertions
  expect_s3_class(g, "gipsmult")
  expect_type(g, "list")
  expect_length(g, 1)

  # Check attributes
  expect_false(is.null(attr(g, "Ss")))
  expect_false(is.null(attr(g, "numbers_of_observations")))
  expect_false(is.null(attr(g, "D_matrices")))
})

test_that("gipsmult validation catches input errors", {
  # Test: Validation catches mismatch between matrices and observations
  expect_error(
    gipsmult(Ss, c(10L)),
    regexp = "cannot be created"
  )

  # Test: Validation catches matrices with different shapes
  S_wrong <- matrix(c(1), nrow = 1)
  Ss_mismatch <- list(S1, S_wrong)

  expect_error(
    gipsmult(Ss_mismatch, ns),
    regexp = "cannot be created"
  )

  # Test: Validation catches non-square matrices
  S_rect <- matrix(1:6, nrow = 2)

  expect_error(
    gipsmult(list(S_rect), c(10L)),
    regexp = "cannot be created"
  )
})

# ==============================================================================
# 2. Print Method Tests
# ==============================================================================

test_that("print.gipsmult outputs expected text", {
  # Setup
  g <- gipsmult(Ss, ns)

  # In testthat, we use expect_output instead of capture.output
  expect_output(print(g), regexp = "permutation")
})

# ==============================================================================
# 3. Plot Method Tests
# ==============================================================================

test_that("plot.gipsmult works correctly", {
  # Setup
  g <- gipsmult(Ss, ns)

  # Test successful plot execution
  # We use a temporary pdf to prevent RPlots.pdf creation or window popping up
  pdf(file = NULL)
    expect_no_error(plot(g, type = "heatmap"))
  dev.off()

  # Test error for invalid type
  expect_error(
    plot(g, type = "non_existent_type"),
    regexp = "must be one of"
  )
})