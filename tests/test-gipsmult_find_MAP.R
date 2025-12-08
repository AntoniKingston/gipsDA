library(tinytest)

# --- Helper: Setup Data ---
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
S2 <- matrix(c(2, -0.2, -0.2, 2), nrow = 2, byrow = TRUE)
Ss <- list(S1, S2)
ns <- c(10L, 15L) # Integers

g <- gipsmult(Ss, ns)

# ==============================================================================
# 1. Optimizer Functionality Tests
# ==============================================================================

g_bf <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)

expect_true(inherits(g_bf, "gipsmult"))
opt_info <- attr(g_bf, "optimization_info")
expect_true(is.list(opt_info))
expect_equal(opt_info$optimization_algorithm_used, "brute_force")
expect_true(opt_info$did_converge)


g_hc <- find_MAP(g, optimizer = "HC", max_iter = 5, show_progress_bar = FALSE)

expect_true(inherits(g_hc, "gipsmult"))
opt_info_hc <- attr(g_hc, "optimization_info")
expect_equal(opt_info_hc$optimization_algorithm_used, "hill_climbing")


g_mh <- find_MAP(g, optimizer = "MH", max_iter = 5, show_progress_bar = FALSE)

expect_true(inherits(g_mh, "gipsmult"))
opt_info_mh <- attr(g_mh, "optimization_info")
expect_equal(opt_info_mh$optimization_algorithm_used, "Metropolis_Hastings")


# ==============================================================================
# 2. Probability Calculation Tests
# ==============================================================================

g_probs <- find_MAP(g, optimizer = "BF", return_probabilities = TRUE, show_progress_bar = FALSE)
probs <- get_probabilities_from_gipsmult(g_probs)

expect_true(is.numeric(probs))
expect_true(length(probs) > 0)
expect_equal(sum(probs), 1, tolerance = 1e-5)


expect_error(
  get_probabilities_from_gipsmult(g),
  pattern = "not optimized"
)


# ==============================================================================
# 3. Argument Validation (Error Handling)
# ==============================================================================

expect_error(
  find_MAP(g, optimizer = "NON_EXISTENT"),
  pattern = "must be one of"
)

expect_error(
  find_MAP(g, optimizer = "MH", max_iter = NA),
  pattern = "max_iter = NA"
)

expect_error(
  find_MAP(g, optimizer = "HC", max_iter = 5, return_probabilities = TRUE),
  pattern = "Probabilities can only be returned"
)

expect_error(
  find_MAP(g, optimizer = "HC", max_iter = Inf, show_progress_bar = TRUE),
  pattern = "not yet supported"
)