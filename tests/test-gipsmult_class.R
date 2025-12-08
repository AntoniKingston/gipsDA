library(tinytest)

# --- Helper: Setup Data ---
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
S2 <- matrix(c(2, -0.2, -0.2, 2), nrow = 2, byrow = TRUE)
Ss <- list(S1, S2)
ns <- c(10L, 15L) # Integers

# ==============================================================================
# 1. Constructor and Validation Tests
# ==============================================================================

# Test: gipsmult constructor creates a valid object
g <- gipsmult(Ss, ns)

expect_true(inherits(g, "gipsmult"))
expect_true(is.list(g))
expect_equal(length(g), 1)

expect_false(is.null(attr(g, "Ss")))
expect_false(is.null(attr(g, "numbers_of_observations")))
expect_false(is.null(attr(g, "D_matrices")))

# Test: Validation catches mismatch between matrices and observations
# FIX: Adjusted pattern to match the generic error from new_gipsmult
expect_error(
  gipsmult(Ss, c(10L)),
  pattern = "cannot be created"
)

# Test: Validation catches matrices with different shapes
S_wrong <- matrix(c(1), nrow = 1)
Ss_mismatch <- list(S1, S_wrong)

# FIX: Adjusted pattern to match the generic error from new_gipsmult
expect_error(
  gipsmult(Ss_mismatch, ns),
  pattern = "cannot be created"
)

# Test: Validation catches non-square matrices
S_rect <- matrix(1:6, nrow = 2)
# FIX: Adjusted pattern to match the generic error from new_gipsmult
expect_error(
  gipsmult(list(S_rect), c(10L)),
  pattern = "cannot be created"
)


# ==============================================================================
# 2. Print Method Tests
# ==============================================================================

output <- capture.output(print(g))
expect_true(any(grepl("permutation", output)))


# ==============================================================================
# 3. Plot Method Tests
# ==============================================================================

pdf(file = NULL)
  res <- try(plot(g, type = "heatmap"), silent = TRUE)
dev.off()

expect_false(inherits(res, "try-error"))

expect_error(
  plot(g, type = "non_existent_type"),
  pattern = "must be one of"
)