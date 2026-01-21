# ==============================================================================
# Global Setup for this file
# ==============================================================================
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
Ss <- list(S1)
ns <- c(100L) # Integer

# ==============================================================================
# 1. Calculation Logic Tests
# ==============================================================================

test_that("log_posteriori_of_gipsmult returns a valid scalar", {
  # Setup
  g <- gipsmult(Ss, ns) # By default was_mean_estimated = TRUE

  val <- log_posteriori_of_gipsmult(g)

  # Assertions
  expect_type(val, "double") # check for numeric type
  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("internal log_posteriori_of_perm matches the wrapper function", {
  # Setup
  g <- gipsmult(Ss, ns)
  g2 <- gipsmult(Ss, ns)
  attr(g2, "was_mean_estimated") <- FALSE
  val <- log_posteriori_of_gipsmult(g)
  val2 <- log_posteriori_of_gipsmult(g2)

  # Internal extraction
  D_mats <- attr(g, "D_matrices")
  delta <- attr(g, "delta")
  perm <- g[[1]]

  # FIX: Because gipsmult() sets was_mean_estimated = TRUE by default,
  # the log_posteriori_of_gipsmult function internally subtracts 1 from ns.
  # To match the result manually, we must also subtract 1 here.

  # Note: Assuming log_posteriori_of_perm is available in the testing namespace
  val_internal <- log_posteriori_of_perm(perm, Ss, ns - 1L, delta, D_mats)
  val_internal2 <- log_posteriori_of_perm(perm, Ss, ns, delta, D_mats)

  # Assertions
  expect_type(val_internal, "double")
  expect_equal(val, val_internal)
  expect_equal(val2, val_internal2)
})
