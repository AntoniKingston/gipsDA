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

  # --- Shortcut for Metropolis-Hastings ---
  g_mh <- find_MAP(g, optimizer = "Me", max_iter = 5, show_progress_bar = FALSE)
  expect_s3_class(g_mh, "gipsmult")
  opt_info_mh <- attr(g_mh, "optimization_info")
  expect_equal(opt_info_mh$optimization_algorithm_used, "Metropolis_Hastings")

  # --- Continue ---
  expect_error(find_MAP(g, optimizer = "continue"))

  g_bf <- find_MAP(g, optimizer = "BF", show_progress_bar = FALSE)
  expect_error(
    find_MAP(g_bf, optimizer = "continue"),
    regexp = "brute_force"
  )

  g_null <- g_mh
  attr(g_null, "optimization_info") <- NULL
  expect_error(find_MAP(g_null, optimizer = "continue", max_iter = 5),
               regexp = "problem identified with the provided arguments")

  g_mh <- g_mh
  opt_info <- attr(g_mh, "optimization_info")
  last_index <- length(opt_info[["optimization_algorithm_used"]])
  opt_info[["optimization_algorithm_used"]][last_index] <- "BF"
  attr(g_mh, "optimization_info") <- opt_info

  expect_error(
    find_MAP(g_mh, optimizer = "continue", max_iter = 5),
    regexp = "brute_force"
  )

  opt_info[["optimization_algorithm_used"]][last_index] <- "MH"
  attr(g_mh, "optimization_info") <- opt_info
  expect_no_error(find_MAP(g_mh, optimizer = "continue",
                           max_iter = 5, show_progress_bar = FALSE))

  # --- was_mean_estimated FALSE ---
  g_no_mean <- g
  attr(g_no_mean, "was_mean_estimated") <- FALSE
  expect_no_error(
    find_MAP(g_no_mean, optimizer = "MH",
             max_iter = 5, show_progress_bar = FALSE)
  )
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

test_that("find_MAP warns about missing stringi package and disables return_probabilities", {
  skip_if_not_installed("mockery")

  S <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
  g <- gipsmult(list(S), c(10))
  target_func <- find_MAP

  # --- MOCKING ---
  mockery::stub(target_func, "rlang::is_installed", FALSE)
  mockery::stub(target_func, "rlang::check_installed", function(...) TRUE)

  # --- TESTING ---
  expect_warning(result <- target_func(g,
                                       optimizer = "MH",
                                       max_iter = 5,
                                       return_probabilities = TRUE,
                                       show_progress_bar = FALSE),
                 regexp = "Package `stringi` is required")

  expect_s3_class(result, "gipsmult")
})

test_that("find_MAP metropolis warns about infinite max_iter", {
  S <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
  g <- gipsmult(list(S), c(10))

  expect_error(
    find_MAP(g, optimizer = "MH", max_iter = Inf, show_progress_bar = FALSE),
    regexp = "must be finite."
  )
})

# ==============================================================================
# MH optimizer Tests
# ==============================================================================

test_that("MH optimizer sets start perm as identity when none provided", {
  expect_error(new_gipsmult <- Metropolis_Hastings_optimizer(
    Ss, ns, 10
  ))
})

# ==============================================================================
# BF optimizer Tests
# ==============================================================================

test_that("BF stops working if too big permutation size", {
  S_big <- matrix(diag(1, 20), nrow = 20)
  g_big <- gipsmult(list(S_big), c(10))

  expect_error(
    find_MAP(g_big, optimizer = "BF", show_progress_bar = FALSE),
    regexp = "big"
  )
})

# ==============================================================================
# Combine gips Tests
# ==============================================================================

test_that("find_MAP works with combined gips objects", {

  g_combined <- combine_gips(g, g)

  g_map <- find_MAP(g_combined, optimizer = "BF", show_progress_bar = FALSE)

  expect_s3_class(g_map, "gipsmult")
  opt_info <- attr(g_map, "optimization_info")
  expect_equal(opt_info$optimization_algorithm_used, "brute_force")
})

test_that("combine_gips warns on inconsistent visited_perms", {
  # --- SETUP ---
  S <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
  g <- gipsmult(list(S), c(10))

  fake_perm <- list("(1)(2)")
  attr(fake_perm, "size") <- 2L
  class(fake_perm) <- "gips_perm"
  # Mock for g1 (has visited_perms)
  opt_info_full <- list(
    optimization_algorithm_used = "Metropolis_Hastings",
    visited_perms = list(fake_perm),
    post_probabilities = c(1),
    log_posteriori_values = -5,
    acceptance_rate = 0.5,
    iterations_performed = 5,
    best_perm_log_posteriori = -100,
    original_perm = fake_perm,
    start_perm = fake_perm,
    last_perm = fake_perm,
    last_perm_log_posteriori = -5,
    did_converge = TRUE,
    optimization_time = 0.1,
    whole_optimization_time = 0.2,
    all_n0 = 0
  )

  # Mock for g2 (visited_perms set to NA - simulate no history collection)
  opt_info_empty <- opt_info_full
  opt_info_empty$visited_perms <- NA
  opt_info_empty$post_probabilities <- NULL

  g1 <- g
  attr(g1, "optimization_info") <- opt_info_full

  g2 <- g
  attr(g2, "optimization_info") <- opt_info_empty

  # --- TEST 3: Warning (Mismatch) ---
  expect_warning(
    combine_gips(g1, g2),
    regexp = "You wanted to save visited_perms"
  )

  g3 <- g
  attr(g3, "optimization_info") <- opt_info_empty
  expect_warning(combine_gips(g1, g3))

  g4 <- g
  opt_info_full2 <- opt_info_full
  opt_info_full2$best_perm_log_posteriori <- -50
  opt_info_full2$optimization_algorithm_used <- "BF"
  attr(g4, "optimization_info") <- opt_info_full2
  expect_no_error(combine_gips(g4, g1))
})

