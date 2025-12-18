# ==============================================================================
# Global Setup for this file
# ==============================================================================
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
S2 <- matrix(c(2, -0.2, -0.2, 2), nrow = 2, byrow = TRUE)
Ss <- list(S1, S2)
ns <- c(10L, 15L) # Integers

g <- gipsmult(Ss, ns)

# ==============================================================================
# 1. Optimizer Functionality Tests
# ==============================================================================

test_that("find_MAP works with different optimizers", {

  # --- Brute Force (BF) ---
  g_bf <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)

  expect_s3_class(g_bf, "gipsmult")

  opt_info <- attr(g_bf, "optimization_info")
  expect_type(opt_info, "list")
  expect_equal(opt_info$optimization_algorithm_used, "brute_force")
  expect_true(opt_info$did_converge)

  # --- Hill Climbing (HC) ---
  g_hc <- find_MAP(g, optimizer = "HC", max_iter = 5, show_progress_bar = FALSE)

  expect_s3_class(g_hc, "gipsmult")
  opt_info_hc <- attr(g_hc, "optimization_info")
  expect_equal(opt_info_hc$optimization_algorithm_used, "hill_climbing")

  # --- Metropolis-Hastings (MH) ---
  g_mh <- find_MAP(g, optimizer = "MH", max_iter = 5, show_progress_bar = FALSE)

  expect_s3_class(g_mh, "gipsmult")
  opt_info_mh <- attr(g_mh, "optimization_info")
  expect_equal(opt_info_mh$optimization_algorithm_used, "Metropolis_Hastings")
})

# ==============================================================================
# 2. Probability Calculation Tests
# ==============================================================================

test_that("probability calculations work correctly", {

  # Success case
  g_probs <- find_MAP(g, optimizer = "BF", return_probabilities = TRUE, show_progress_bar = FALSE)
  probs <- get_probabilities_from_gipsmult(g_probs)

  expect_type(probs, "double") # expect_true(is.numeric(...))
  expect_gt(length(probs), 0)
  expect_equal(sum(probs), 1, tolerance = 1e-5)

  # Failure case (not optimized)
  expect_error(
    get_probabilities_from_gipsmult(g),
    regexp = "not optimized"
  )
})

# ==============================================================================
# 3. Argument Validation (Error Handling)
# ==============================================================================

test_that("find_MAP catches invalid arguments", {

  expect_error(
    find_MAP(g, optimizer = "NON_EXISTENT"),
    regexp = "must be one of"
  )

  expect_error(
    find_MAP(g, optimizer = "MH", max_iter = NA),
    regexp = "max_iter = NA"
  )

  expect_error(
    find_MAP(g, optimizer = "HC", max_iter = 5, return_probabilities = TRUE),
    regexp = "Probabilities can only be returned"
  )

  expect_error(
    find_MAP(g, optimizer = "HC", max_iter = Inf, show_progress_bar = TRUE),
    regexp = "not yet supported"
  )
})