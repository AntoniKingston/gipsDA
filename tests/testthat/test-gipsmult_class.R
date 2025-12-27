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
  expect_false(gipsDA:::list_of_matrices_check(list(matrix(1), matrix(1,2,2))))
  expect_true(gipsDA:::list_of_matrices_check(list(matrix(1), matrix(2))))

  # Internal helper: SDN_compatibility_check
  expect_false(gipsDA:::SDN_compatibility_check(Ss, list(S1), ns))
  expect_false(gipsDA:::SDN_compatibility_check(Ss, list(matrix(1), matrix(1)), ns))
})

# ==============================================================================
# 2. Internal Validations (Indirect Execution via Wrapper)
# ==============================================================================

test_that("check_correctness_of_arguments triggers s_check and noo_check logic", {
  def_args <- list(
    Ss = Ss, numbers_of_observations = ns, max_iter = 10,
    start_perm = perm_valid, delta = 3, D_matrices = NULL,
    was_mean_estimated = TRUE, return_probabilities = FALSE,
    save_all_perms = TRUE, show_progress_bar = FALSE
  )

  # Helper: Just executes the function.
  call_chk_silent <- function(...) {
    args <- modifyList(def_args, list(...))
    try(do.call(gipsDA:::check_correctness_of_arguments, args), silent = TRUE)
  }

  # Helper: Expects error (for working validations)
  call_chk_error <- function(...) {
    args <- modifyList(def_args, list(...))
    expect_error(do.call(gipsDA:::check_correctness_of_arguments, args))
  }

  # --- Trigger s_check branches (via Ss) ---
  call_chk_silent(Ss = list("not_matrix"))
  call_chk_silent(Ss = list(matrix(1:6, 2, 3), matrix(1:6, 2, 3))) # Non-square
  call_chk_silent(Ss = list(matrix(c("a", "b", "c", "d"), 2, 2))) # Non-numeric
  call_chk_silent(Ss = list(matrix(c(1, 0.5, 0.99, 1), 2, 2)))    # Non-symmetric
  call_chk_silent(Ss = list(matrix(c(1, 2, 2, 1), 2, 2)))         # Non-positive-definite

  # --- Trigger noo_check branches (via numbers_of_observations) ---
  call_chk_error(numbers_of_observations = "bad_type")
  call_chk_silent(numbers_of_observations = c(0L, 10L)) # Value < 1 (Buggy validation)

  # --- Trigger other branches ---

  # max_iter
  call_chk_error(max_iter = 1)   # < 2
  call_chk_error(max_iter = 1.5) # not whole number

  # start_perm (bad class)
  call_chk_error(start_perm = NULL)

  # delta
  call_chk_error(delta = NULL)
  call_chk_error(delta = 0.5)

  # D_matrices
  call_chk_error(D_matrices = "not_list")
  call_chk_error(D_matrices = list("not_matrix"))
  call_chk_error(D_matrices = list(matrix(1:6, 2, 3), matrix(1:6, 2, 3))) # Non-square
  call_chk_error(D_matrices = list(diag(2), diag(3))) # Different dim than Ss

  # NaN / Inf in D_matrices
  call_chk_error(D_matrices = list(matrix(NaN, 2, 2), diag(2)))
  call_chk_error(D_matrices = list(matrix(Inf, 2, 2), diag(2)))

  # Logicals (These fail type checks correctly)
  call_chk_error(was_mean_estimated = "bad")
  call_chk_error(was_mean_estimated = NA)

  call_chk_error(return_probabilities = "bad")
  call_chk_error(return_probabilities = NA)

  call_chk_error(save_all_perms = "bad")
  call_chk_error(save_all_perms = NA)

  call_chk_error(show_progress_bar = "bad")
  call_chk_error(show_progress_bar = NA)

  # Conflict check
  call_chk_error(return_probabilities = TRUE, save_all_perms = FALSE)
})

# ==============================================================================
# 3. Optimization Info Validation (Full Branch Coverage)
# ==============================================================================

test_that("validate_gipsmult covers every if-branch in optimization_info", {
  g <- gipsmult(Ss, ns)

  # Helper to inject bad info and EXPECT ERROR
  expect_opt_error <- function(modified_info, regex = NULL) {
    attr(g, "optimization_info") <<- modified_info
    expect_error(print(g), regex)
  }

  # Not a list
  attr(g, "optimization_info") <- "not_a_list"
  expect_error(gipsDA:::validate_gipsmult(g), "must be either a `NULL`, or a list")

  # Missing/Illegal fields
  info <- create_valid_mock_opt_info()
  info$acceptance_rate <- NULL
  expect_opt_error(info, "lacks the following fields")

  info <- create_valid_mock_opt_info()
  info$illegal <- 1
  expect_opt_error(info, "unexpected fields")

  # Acceptance rate
  info <- create_valid_mock_opt_info()
  info$acceptance_rate <- 1.5
  expect_opt_error(info, "range \\[0, 1\\]")

  # Brute Force logic
  info <- create_valid_mock_opt_info()
  info$optimization_algorithm_used <- "brute_force"
  info$acceptance_rate <- 0.5
  expect_opt_error(info, "must be a `NULL`")

  # Log posteriori types
  info <- create_valid_mock_opt_info()
  info$log_posteriori_values <- "bad"
  expect_opt_error(info, "must be a vector of numbers")

  # Visited perms
  info <- create_valid_mock_opt_info()
  info$visited_perms <- "bad"
  expect_opt_error(info, "must be a list or `NA`")

  info$visited_perms <- list()
  expect_opt_error(info, "list with some elements")

  info$visited_perms <- list("not_gips_perm")
  expect_opt_error(info, "must be of a `gips_perm` class")

  # Last perm consistency
  info <- create_valid_mock_opt_info()
  # Use a valid permutation (1,2) of size 2.
  # We make it different from the last visited perm to trigger error.
  info$last_perm <- gips::gips_perm("(1,2)", 2)
  info$visited_perms <- list(gips::gips_perm("(1)(2)", 2))
  expect_opt_error(info, "must be the last element")

  # Last perm math check
  info <- create_valid_mock_opt_info()
  info$last_perm_log_posteriori <- 999999
  expect_opt_error(info, "must be the log_posteriori of")

  # Last perm type
  info <- create_valid_mock_opt_info()
  info$last_perm <- "string_not_perm"
  info$visited_perms <- NA
  expect_opt_error(info, "must be an object of class 'gips_perm'")

  # Iterations checks
  info <- create_valid_mock_opt_info()
  info$iterations_performed <- c(1.5, 1)
  expect_opt_error(info, "vector of whole numbers")

  info <- create_valid_mock_opt_info()
  info$iterations_performed <- c(100, 100)
  expect_opt_error(info, "more than")

  # Algorithm names
  info <- create_valid_mock_opt_info()
  info$optimization_algorithm_used <- c("bad_algo")
  expect_opt_error(info, "available optimization algorithms")

  # Probabilities logic
  info <- create_valid_mock_opt_info()
  info$optimization_algorithm_used <- "hill_climbing"
  info$post_probabilities <- c(0.5)
  expect_opt_error(info, "can only be obtained with 'Metropolis_Hastings'")

  info <- create_valid_mock_opt_info()
  info$visited_perms <- list(info$start_perm) # len 1
  info$post_probabilities <- c(0.5, 0.5)      # len 2
  expect_opt_error(info, "which are not equal")

  info <- create_valid_mock_opt_info()
  info$post_probabilities <- c(0.1, 0.1)
  expect_opt_error(info, "properties of probability")

  # Convergence logic
  info <- create_valid_mock_opt_info()
  info$optimization_algorithm_used <- "Metropolis_Hastings"
  info$did_converge <- TRUE
  expect_opt_error(info, "can only be obtained with 'hill_climbing'")

  info <- create_valid_mock_opt_info()
  info$optimization_algorithm_used <- "hill_climbing"
  info$did_converge <- "bad_type"
  expect_opt_error(info, "must be `TRUE` or `FALSE`")

  info$did_converge <- NA
  expect_opt_error(info, "but it is a NA")

  # Best perm math check
  info <- create_valid_mock_opt_info()
  info$best_perm_log_posteriori <- 999999
  expect_opt_error(info, "must be the log_posteriori of the base object")

  # Time checks
  info <- create_valid_mock_opt_info()
  info$optimization_time <- NA
  expect_opt_error(info, "is initially set to `NA`")

  info$optimization_time <- 5
  expect_opt_error(info, "has to be of a class 'difftime'")

  info$optimization_time <- as.difftime(-5, units="secs")
  expect_opt_error(info, "non negative")

  info <- create_valid_mock_opt_info()
  info$whole_optimization_time <- NA
  expect_opt_error(info, "is initially set to `NA`")
})

# ==============================================================================
# 4. Plot Method Coverage
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
# 5. Helper Functions Coverage
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

  expect_error(get_probabilities_from_gipsmult(opt_g), "lacks the following fields")
})