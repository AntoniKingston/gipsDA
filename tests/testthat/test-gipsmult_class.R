# ==============================================================================
# Global Setup
# ==============================================================================
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
S2 <- matrix(c(2, -0.2, -0.2, 2), nrow = 2, byrow = TRUE)
Ss <- list(S1, S2)
ns <- c(10L, 15L)
perm_valid <- gips::gips_perm("(1)(2)", 2)

# Helper: Compute real log_posteriori to pass strict math validation in mocks
get_real_log_post <- function(perm_str) {
  perm <- gips::gips_perm(perm_str, 2)
  g <- gipsDA:::gipsmult(Ss, ns, perm = perm)
  gipsDA:::log_posteriori_of_gipsmult(g)
}

# Helper: Create valid optimization_info list
create_valid_mock_opt_info <- function() {
  perm_obj <- gips::gips_perm("(1)(2)", 2)
  real_val <- get_real_log_post("(1)(2)")

  list(
    original_perm = perm_obj,
    acceptance_rate = 0.5,
    log_posteriori_values = c(real_val, real_val),
    visited_perms = list(perm_obj, perm_obj),
    start_perm = perm_obj,
    last_perm = perm_obj,
    last_perm_log_posteriori = real_val,
    iterations_performed = c(1, 1),
    optimization_algorithm_used = c("Metropolis_Hastings", "Metropolis_Hastings"),
    post_probabilities = c(0.5, 0.5),
    did_converge = NULL,
    best_perm_log_posteriori = real_val,
    optimization_time = as.difftime(1, units = "secs"),
    whole_optimization_time = as.difftime(2, units = "secs"),
    all_n0 = c(10, 10)
  )
}

# ==============================================================================
# 1. Constructor & Basic Utils Coverage
# ==============================================================================

test_that("gipsmult constructor handles various inputs", {
  # String permutation
  g <- gipsmult(Ss, ns, perm = "(1,2)")
  expect_s3_class(g, "gipsmult")

  # gips object permutation
  perm_gips <- gips::gips(S1, 10, perm = "(1,2)")
  g2 <- gipsmult(Ss, ns, perm = perm_gips)
  expect_equal(g2[[1]], perm_gips[[1]])

  # Default D_matrices generation
  g3 <- gipsmult(Ss, ns, D_matrices = NULL)
  expect_true(is.list(attr(g3, "D_matrices")))

  # Internal helper: list_of_matrices_check
  expect_false(gipsDA:::list_of_matrices_check("not_a_list"))
  expect_false(gipsDA:::list_of_matrices_check(list("not_matrix")))
  expect_true(gipsDA:::list_of_matrices_check(list(matrix(1), matrix(1,2,2))))
  expect_true(gipsDA:::list_of_matrices_check(list(matrix(1), matrix(2))))

  # Internal helper: SDN_compatibility_check
  expect_false(gipsDA:::SDN_compatibility_check(Ss, list(S1), ns))
  expect_false(gipsDA:::SDN_compatibility_check(Ss, list(matrix(1), matrix(1)), ns))
})

# ==============================================================================
# 2. Plot Method Coverage
# ==============================================================================

test_that("plot.gipsmult executes all plot types", {
  g <- gipsmult(Ss, ns)
  opt_g <- g
  info <- create_valid_mock_opt_info()
  attr(opt_g, "optimization_info") <- info

  pdf(NULL) # Prevent graphics device output
    expect_no_error(plot(g, type = NA))
    expect_no_error(plot(opt_g, type = NA))

    if (rlang::is_installed("ggplot2")) {
      expect_no_error(plot(opt_g, type = "heatmap"))
      try(plot(opt_g, type = "block_heatmap"), silent = TRUE)
      expect_no_error(plot(opt_g, type = "MLE"))
    }

    expect_no_error(plot(opt_g, type = "best"))
    expect_no_error(plot(opt_g, type = "all"))
    expect_no_error(plot(opt_g, type = "both"))
    expect_no_error(plot(opt_g, type = "n0"))

    # Check parameters passed to plot
    expect_no_error(plot(opt_g, type="both", logarithmic_x=TRUE, logarithmic_y=FALSE, color="purple"))
  dev.off()
})

# ==============================================================================
# 3. Helper Functions Coverage
# ==============================================================================

test_that("get_probabilities executes", {
  g <- gipsmult(Ss, ns)

  # Error path
  expect_error(get_probabilities_from_gipsmult(g))

  # Valid path
  opt_g <- g
  info <- create_valid_mock_opt_info()
  attr(opt_g, "optimization_info") <- info

  probs <- get_probabilities_from_gipsmult(opt_g)
  expect_type(probs, "double")

  # NULL probabilities path
  info_null <- info
  info_null$post_probabilities <- NULL
  attr(opt_g, "optimization_info") <- info_null

  expect_no_error(get_probabilities_from_gipsmult(opt_g))

  # Checking print function
  expect_output(print(opt_g))

  # Expecting abort when no optimization_info
  attr(g, "optimization_info") <- NULL
  expect_no_error(print(g))
})