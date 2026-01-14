find_MAP <- function(g, max_iter = NA, optimizer = NA,
                     show_progress_bar = TRUE,
                     save_all_perms = FALSE,
                     return_probabilities = FALSE) {

  possible_optimizers <- c(
    "MH", "Metropolis_Hastings", "HC", "hill_climbing",
    "BF", "brute_force", "full", "continue"
  )


  # get a chosen optimizer, even with part of the name:
  chosen_optimizer_number <- pmatch(optimizer, possible_optimizers)

  if (is.na(chosen_optimizer_number)) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = paste0(
        "`optimizer` must be one of: c('",
        paste0(possible_optimizers, collapse = "', '"), "')."
      ),
      "x" = paste0("You provided `optimizer == '", optimizer, "'`."),
      "i" = "Did You misspell the optimizer name?"
    ))
  }

  if (optimizer != possible_optimizers[chosen_optimizer_number]) {
    rlang::inform(c(
      "You provided a shortcut for the optimization method's name:",
      "i" = paste0("You provided `optimizer == '", optimizer, "'`"),
      "i" = paste0(
        "This will be changed to `optimizer == '",
        possible_optimizers[chosen_optimizer_number], "'`"
      )
    ))

    optimizer <- possible_optimizers[chosen_optimizer_number]
  }

  if (!(optimizer %in% c("BF", "brute_force", "full")) &&
    is.na(max_iter)) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`max_iter = NA` can be provided only for `optimizer` one of: c('BF', 'brute_force', 'full'). For any other, `max_iter` must be a whole number, strictly bigger than 1.",
      "x" = paste0("You provided `optimizer == ", optimizer, "` and `max_iter = NA`."),
      "i" = "Did You forget to set the `max_iter`?",
      "i" = "Did You misspell the optimizer name?"
    ))
  }

  continue_optimization <- (optimizer == "continue")
  if (continue_optimization) {
    if (is.null(attr(g, "optimization_info"))) {
      rlang::abort(c("There was a problem identified with the provided arguments:",
        "i" = "`optimizer == 'continue'` can be provided only with an optimized gips object `g`.",
        "x" = "You provided `optimizer == 'continue'`, but the gips object `g` is not optimized.",
        "i" = "Did You provide wrong `gips` object?",
        "i" = "Did You want to call another optimizer like 'MH' or 'HC'?"
      ))
    }

    optimizer <- attr(g, "optimization_info")[["optimization_algorithm_used"]][length(attr(g, "optimization_info")[["optimization_algorithm_used"]])] # this is the last used optimizer
    if (optimizer %in% c("BF", "brute_force", "full")) {
      rlang::abort(c("There was a problem identified with the provided arguments:",
        "i" = "`optimizer == 'continue'` cannot be provided after optimizating with `optimizer == 'brute_force'`, because the whole space was already browsed.",
        "x" = "You provided `optimizer == 'continue'`, but the gips object `g` was optimized with brute_force optimizer. Better permutation will not be found."
      ))
    }
  }

  if (!(optimizer %in% c("MH", "Metropolis_Hastings", "BF", "brute_force", "full")) && return_probabilities) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "Probabilities can only be returned with the `optimizer == 'Metropolis_Hastings'` or `optimizer == 'brute_force'`",
      "x" = "You provided both `!(optimizer %in% c('Metropolis_Hastings', 'brute_force'))` and `return_probabilities == TRUE`!",
      "i" = "Did You want to use `optimizer == 'Metropolis_Hastings'` or `optimizer == 'brute_force'`, or `return_probabilities == FLASE`?"
    ))
  }

  if (return_probabilities && optimizer %in% c("MH", "Metropolis_Hastings")) {
    rlang::check_installed("stringi",
      reason = "to return probabilities in `find_MAP(optimizer = 'Metropolis_Hastings', return_probabilities = TRUE)`; without this package, probabilities cannot be returned"
    )
    if (!rlang::is_installed("stringi")) {
      rlang::warn(c("There was a problem with return_probabilities:",
        "i" = "Package `stringi` is required to successfully call `find_MAP(optimizer = 'Metropolis_Hastings', return_probabilities = TRUE)`.",
        "x" = "You do not have package `stringi` installed.",
        "i" = "Optimization will proceed as `find_MAP(optimizer = 'Metropolis_Hastings', return_probabilities = FALSE)`."
      ))

      return_probabilities <- FALSE
    }
  }

  # inform that user can consider "BF"
  if ((optimizer %in% c("MH", "Metropolis_Hastings")) &&
    (max_iter * 10 >= prod(1:ncol(attr(g, "Ss")[[1]]))) &&
    is.finite(max_iter)) { # infinite max_iter is illegal, but additional check will not hurt
    rlang::inform(c(
      paste0(
        "You called optimization with Metropolis_Hastings algorith with ",
        max_iter, " iterations."
      ),
      "i" = paste0(
        "Consider using `optimizer = 'brute_force'`, because it will use ",
        ncol(attr(g, "Ss")), "! (factorial) = ", prod(1:ncol(attr(g, "Ss")[[1]])),
        " iterations and will browse all permutations, therefore it will definitely find the maximum a posteriori estimator."
      )
    ))
  }

  # extract parameters
  Ss <- attr(g, "Ss")
  numbers_of_observations <- attr(g, "numbers_of_observations")
  if (continue_optimization) {
    start_perm <- attr(g, "optimization_info")[["last_perm"]]
  } else {
    start_perm <- g[[1]]
  }
  delta <- attr(g, "delta")
  D_matrices <- attr(g, "D_matrices")
  was_mean_estimated <- attr(g, "was_mean_estimated")

  if (was_mean_estimated) { # one degree of freedom is lost; we will return this 1 to numbers_of_observations after optimization in `combine_gips()`
    edited_numbers_of_observations <- numbers_of_observations - 1
  } else {
    edited_numbers_of_observations <- numbers_of_observations
  }

  start_time <- Sys.time()

  if (optimizer %in% c("MH", "Metropolis_Hastings")) {
    gips_optimized <- Metropolis_Hastings_optimizer(
      Ss = Ss, numbers_of_observations = edited_numbers_of_observations,
      max_iter = max_iter, start_perm = start_perm,
      delta = delta, D_matrices = D_matrices,
      return_probabilities = return_probabilities,
      save_all_perms = save_all_perms,
      show_progress_bar = show_progress_bar
    )
  } else if (optimizer %in% c("HC", "hill_climbing")) {
    gips_optimized <- hill_climbing_optimizer(
      Ss = Ss, numbers_of_observations = edited_numbers_of_observations,
      max_iter = max_iter, start_perm = start_perm,
      delta = delta, D_matrices = D_matrices,
      save_all_perms = save_all_perms,
      show_progress_bar = show_progress_bar
    )
  } else if (optimizer %in% c("BF", "brute_force", "full")) {
    gips_optimized <- brute_force_optimizer(
      Ss = Ss, numbers_of_observations = edited_numbers_of_observations,
      delta = delta, D_matrices = D_matrices,
      return_probabilities = return_probabilities,
      save_all_perms = save_all_perms,
      show_progress_bar = show_progress_bar
    )
  }

  end_time <- Sys.time()
  attr(gips_optimized, "optimization_info")[["optimization_time"]] <- end_time - start_time
  attr(gips_optimized, "optimization_info")[["whole_optimization_time"]] <- end_time - start_time

  structure_constants <- gips::get_structure_constants(gips_optimized[[1]])
  n0 <- max(structure_constants[["r"]] * structure_constants[["d"]] / structure_constants[["k"]])
  if (attr(g, "was_mean_estimated")) { # correction for estimating the mean
    n0 <- n0 + 1
    attr(gips_optimized, "optimization_info")[["all_n0"]] <- attr(gips_optimized, "optimization_info")[["all_n0"]] + 1 # when all_n0 is NA, all_n0 + 1 is also an NA
  }
  if (any(n0 > numbers_of_observations)) {
    rlang::warn(c(
      paste0(
        "The found permutation has n0 = ", n0,
        ", which is bigger than the numbers_of_observations = ",
        numbers_of_observations, "."
      ),
      "i" = "The covariance matrix invariant under the found permutation does not have the likelihood properly defined.",
      "i" = "For a more in-depth explanation, see the 'Project Matrix - Equation (6)' section in the `vignette('Theory', package = 'gips')` or its pkgdown page: https://przechoj.github.io/gips/articles/Theory.html."
    ))
  }

  return(combine_gips(g, gips_optimized))
}


Metropolis_Hastings_optimizer <- function(Ss,
    numbers_of_observations, max_iter, start_perm = NULL,
    delta = 3, D_matrices = NULL, return_probabilities = FALSE,
    save_all_perms = FALSE, show_progress_bar = TRUE) {
  if (is.null(start_perm)) {
    start_perm <- permutations::id
  }

  # check_correctness_of_arguments(
  #   Ss = Ss, numbers_of_observations = numbers_of_observations,
  #   max_iter = max_iter, start_perm = start_perm,
  #   delta = delta, D_matrices = D_matrices, was_mean_estimated = FALSE,
  #   return_probabilities = return_probabilities,
  #   save_all_perms = save_all_perms,
  #   show_progress_bar = show_progress_bar
  # )

  if (!inherits(start_perm, "gips_perm")) {
    start_perm <- gips::gips_perm(start_perm, nrow(Ss[[1]])) # now we know the `Ss` is a matrix
  }

  if (is.infinite(max_iter)) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`max_iter` in `Metropolis_Hastings_optimizer` must be finite.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    ))
  }

  perm_size <- dim(Ss[[1]])[1]
  if (permutations::is.cycle(start_perm)) {
    start_perm <- gips::gips_perm(start_perm, perm_size)
  }
  if (is.null(D_matrices)) {
    D_matrices <- list(rep(diag(nrow = perm_size), length(numbers_of_observations)))
  }

  my_goal_function <- function(perm, i) {
    out_val <- log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm()`, remember to edit `numbers_of_observations` if the mean was estimated!
      Ss = Ss, numbers_of_observations = numbers_of_observations,
      delta = delta, D_matrices = D_matrices
    )

    if (is.nan(out_val) || is.infinite(out_val)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::abort(c(
        "gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(out_val), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrices. If it is not the case for You, please get in touch with us on ISSUE#5.",
        "x" = paste0("The Metropolis Hastings algorithm was stopped after ", i, " iterations.")
      ))
    }

    out_val
  }

  acceptance <- rep(FALSE, max_iter)
  log_posteriori_values <- rep(0, max_iter)
  all_n0 <- rep(0, max_iter)
  if (save_all_perms) {
    visited_perms <- list()
    visited_perms[[1]] <- start_perm
  } else {
    visited_perms <- NA
  }
  current_perm <- start_perm
  all_n0[1] <- gips:::get_n0_from_perm(current_perm, was_mean_estimated = FALSE) # was_mean_estimated will be corrected in find_MAP()

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  }
  log_posteriori_values[1] <- my_goal_function(current_perm, 0)

  found_perm <- start_perm
  found_perm_log_posteriori <- log_posteriori_values[1]

  Uniformly_drawn_numbers <- stats::runif(max_iter, min = 0, max = 1)

  # main loop
  for (i in 1:(max_iter - 1)) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }

    e <- gips:::runif_transposition(perm_size)
    perm_proposal <- gips:::compose_with_transposition(current_perm, e)

    goal_function_perm_proposal <- my_goal_function(perm_proposal, i)

    # if goal_function_perm_proposal > log_posteriori_values[i], then it is true, because Uniformly_drawn_numbers[i] \in [0,1]
    if (Uniformly_drawn_numbers[i] < exp(goal_function_perm_proposal - log_posteriori_values[i])) { # the probability of drawing e such that g' = g*e is the same as the probability of drawing e' such that g = g'*e. This probability is 1/(p choose 2). That means this is Metropolis algorithm, not necessary Metropolis-Hastings.
      current_perm <- perm_proposal
      if (save_all_perms) {
        visited_perms[[i + 1]] <- current_perm
      }
      log_posteriori_values[i + 1] <- goal_function_perm_proposal
      acceptance[i] <- TRUE
      all_n0[i+1] <- gips:::get_n0_from_perm(current_perm, was_mean_estimated = FALSE) # was_mean_estimated will be corrected in find_MAP()

      if (found_perm_log_posteriori < log_posteriori_values[i + 1]) {
        found_perm_log_posteriori <- log_posteriori_values[i + 1]
        found_perm <- current_perm
      }
    } else {
      if (save_all_perms) {
        visited_perms[[i + 1]] <- current_perm
      }
      log_posteriori_values[i + 1] <- log_posteriori_values[i]
      all_n0[i+1] <- all_n0[i]
    }
  }

  if (show_progress_bar) {
    close(progressBar)
  }

  function_calls <- length(log_posteriori_values)

  # visited_perms are already either a list of things or a `NULL` object

  if (return_probabilities) {
    probabilities <- gips:::estimate_probabilities(visited_perms, show_progress_bar)
  } else {
    probabilities <- NULL
  }

  optimization_info <- list(
    "acceptance_rate" = mean(acceptance),
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = visited_perms,
    "start_perm" = start_perm,
    "last_perm" = current_perm,
    "last_perm_log_posteriori" = log_posteriori_values[function_calls],
    "iterations_performed" = i,
    "optimization_algorithm_used" = "Metropolis_Hastings",
    "post_probabilities" = probabilities,
    "did_converge" = NULL,
    "best_perm_log_posteriori" = found_perm_log_posteriori,
    "optimization_time" = NA,
    "whole_optimization_time" = NA,
    "all_n0" = all_n0
  )


  new_gipsmult(
    list(found_perm), Ss, numbers_of_observations,
    delta, D_matrices,
    was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}


hill_climbing_optimizer <- function(Ss,
    numbers_of_observations, max_iter = 5,
    start_perm = NULL, delta = 3, D_matrices = NULL,
    save_all_perms = FALSE, show_progress_bar = TRUE) {
  if (is.null(start_perm)) {
    start_perm <- permutations::id
  }

  # check_correctness_of_arguments(
  #   Ss = Ss, numbers_of_observations = numbers_of_observations,
  #   max_iter = max_iter, start_perm = start_perm,
  #   delta = delta, D_matrices = D_matrices, was_mean_estimated = FALSE,
  #   return_probabilities = FALSE, save_all_perms = save_all_perms,
  #   show_progress_bar = show_progress_bar
  # )

  if (!inherits(start_perm, "gips_perm")) {
    start_perm <- gips::gips_perm(start_perm, nrow(Ss[[1]])) # now we know the `S` is a matrix
  }

  if (show_progress_bar && is.infinite(max_iter)) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "x" = "You tried to run `find_MAP(show_progress_bar=TRUE, max_iter=Inf)`.",
      "i" = "Progress bar is not yet supported for infinite max_iter.",
      "i" = "Do You want to use `show_progress_bar=FALSE` or a finite `max_iter`?",
      "i" = "For more information on progress bar see ISSUE#8."
    ))
  }

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = max_iter, initial = 1)
  }

  perm_size <- dim(Ss[[1]])[1]

  if (is.null(D_matrices)) {
    D_matrices <- list(rep(diag(nrow = perm_size), length(numbers_of_observations)))
  }


  my_goal_function <- function(perm, i) {
    out_val <- log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `numbers_of_observations` if the mean was estimated!
      Ss = Ss, numbers_of_observations = numbers_of_observations,
      delta = delta, D_matrices = D_matrices
    )

    if (is.nan(out_val) || is.infinite(out_val)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::abort(c(
        "gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(out_val), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for D_matrices with huge values. If it is not the case for You, please get in touch with us on ISSUE#5.",
        "x" = paste0("The Hill Climbing algorithm was stopped after ", i, " iterations.")
      ))
    }

    out_val
  }

  goal_function_best_logvalues <- numeric(0)
  log_posteriori_values <- numeric(0)

  # init
  if (save_all_perms) {
    visited_perms <- list()
    visited_perms[[1]] <- start_perm
  } else {
    visited_perms <- NA
  }
  current_perm <- start_perm

  goal_function_best_logvalues[1] <- my_goal_function(current_perm, 0)
  log_posteriori_values[1] <- goal_function_best_logvalues[1]

  # mail loop
  iteration <- 0
  did_converge <- FALSE
  while (iteration <= max_iter - 1) {
    iteration <- iteration + 1
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, iteration)
    }

    best_neighbour <- NULL
    best_neighbour_value <- -Inf
    for (i in 1:(perm_size - 1)) {
      for (j in (i + 1):perm_size) {
        neighbour <- gips:::compose_with_transposition(current_perm, c(i, j))
        neighbour_value <- my_goal_function(neighbour, iteration)
        log_posteriori_values[length(log_posteriori_values) + 1] <- neighbour_value

        if (neighbour_value > best_neighbour_value) {
          best_neighbour_value <- neighbour_value
          best_neighbour <- neighbour
        }
      }
    }

    if (best_neighbour_value > goal_function_best_logvalues[iteration]) {
      goal_function_best_logvalues[iteration + 1] <- best_neighbour_value
      current_perm <- best_neighbour
      if (save_all_perms) {
        visited_perms[[iteration + 1]] <- best_neighbour
      }
    } else {
      did_converge <- TRUE
      break
    }
  }

  last_perm <- current_perm

  if (show_progress_bar) {
    close(progressBar)
  }

  if (!did_converge) {
    rlang::warn(c(paste0("Hill Climbing algorithm did not converge in ", iteration, " iterations!"), # now, iteration == max_iter
      "i" = "We recommend to run the `find_MAP(optimizer = 'continue')` on the acquired output."
    ))
    iteration <- iteration + 1 # the very first was the starting perm
  } else {
    goal_function_best_logvalues <- goal_function_best_logvalues[1:iteration]
    if (show_progress_bar) {
      print(paste0("Algorithm did converge in ", iteration, " iterations"))
    }
  }

  function_calls <- length(log_posteriori_values)

  optimization_info <- list(
    "acceptance_rate" = 1 / choose(perm_size, 2),
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = visited_perms,
    "start_perm" = start_perm,
    "last_perm" = last_perm,
    "last_perm_log_posteriori" = goal_function_best_logvalues[iteration],
    "iterations_performed" = iteration,
    "optimization_algorithm_used" = "hill_climbing",
    "post_probabilities" = NULL,
    "did_converge" = did_converge,
    "best_perm_log_posteriori" = goal_function_best_logvalues[iteration],
    "optimization_time" = NA,
    "whole_optimization_time" = NA,
    "all_n0" = NA
  )


  new_gipsmult(
    list(last_perm), Ss, numbers_of_observations,
    delta, D_matrices,
    was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}


brute_force_optimizer <- function(
    Ss,
    numbers_of_observations,
    delta = 3, D_matrices = NULL,
    return_probabilities = return_probabilities,
    save_all_perms = FALSE, show_progress_bar = TRUE) {
  # check_correctness_of_arguments(
  #   Ss = Ss, numbers_of_observations = numbers_of_observations,
  #   max_iter = 5, start_perm = permutations::id, # max_iter, was_mean_estimated and start_perm are not important for optimization with brute_force
  #   delta = delta, D_matrices = D_matrices, was_mean_estimated = FALSE,
  #   return_probabilities = return_probabilities, save_all_perms = save_all_perms,
  #   show_progress_bar = show_progress_bar
  # )

  perm_size <- dim(Ss[[1]])[1]

  if (perm_size > 18) {
    rlang::abort(c("Optimizer 'brute_force' cannot browse such a big permutational space.",
      "x" = paste0(
        "You provided a space with size ", perm_size,
        "! (factorial), which has ", prod(1:perm_size),
        " elements."
      ),
      "i" = "Do You want to use other optimizer for such a big space? For example 'Metropolis_Hastings' or 'hill_climbing'?"
    ))
  }

  if (perm_size > 9) { # I don't know how to test this without running the optimization...
    rlang::warn(c("Optimizer 'brute_force' will take very long time to browse such a big permutational space.",
      "x" = paste0(
        "You provided a space with size ", perm_size,
        "! (factorial), which has ", prod(1:perm_size),
        " elements."
      ),
      "i" = "Do You want to use other optimizer for such a big space? For example 'Metropolis_Hastings' or 'hill_climbing'?"
    ))
  }

  iterations_to_perform <-
    if ((3 <= perm_size) && (perm_size <= 9)) {
      # Only the generators are interesting for us:
      # We precalculated perm_group_generators only for up to perm_size = 9
      # See ISSUE#21 for more information
      gips:::OEIS_A051625[perm_size]
    } else {
      prod(1:perm_size)
    }

  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = iterations_to_perform, initial = 1)
  }

  if (is.null(D_matrices)) {
    D_matrices <- list(rep(diag(nrow = perm_size), length(numbers_of_observations)))
  }

  my_goal_function <- function(perm, i) {
    out_val <- log_posteriori_of_perm(perm, # We recommend to use the `log_posteriori_of_gips()` function. If You really want to use `log_posteriori_of_perm`, remember to edit `numbers_of_observations` if the mean was estimated!
      Ss = Ss, numbers_of_observations = numbers_of_observations,
      delta = delta, D_matrices = D_matrices
    )

    if (is.nan(out_val) || is.infinite(out_val)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::abort(c(
        "gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(out_val), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrices. If it is not the case for You, please get in touch with us on ISSUE#5.",
        "x" = paste0("The Brute Force algorithm was stopped after ", i, " iterations.")
      ))
    }

    out_val
  }

  # main loop
  all_perms_list <- permutations::allperms(perm_size)
  all_perms_list <- permutations::as.cycle(all_perms_list)
  if ((3 <= perm_size) && (perm_size <= 9)) {
    # Only the generators are interesting for us:
    # perm_group_generators are calculated only for up to perm_size = 9
    # See ISSUE#21 for more information
    all_perms_list <- all_perms_list[gips:::perm_group_generators_list[[perm_size - 2]]]
  }
  log_posteriori_values <- sapply(1:length(all_perms_list), function(i) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }
    this_perm <- permutations::cycle(list(all_perms_list[[i]]))
    my_goal_function(this_perm, i)
  })

  if (show_progress_bar) {
    close(progressBar)
  }

  if (return_probabilities) { # calculate exact probabilities
    probabilities <- gips:::calculate_probabilities(all_perms_list, log_posteriori_values, show_progress_bar)
  } else {
    probabilities <- NULL
  }

  best_perm <- gips::gips_perm(permutations::cycle(list(all_perms_list[[which.max(log_posteriori_values)]])), perm_size)

  if (save_all_perms) {
    visited_perms <- all_perms_list
  } else {
    visited_perms <- NA
  }

  optimization_info <- list(
    "acceptance_rate" = NULL,
    "log_posteriori_values" = log_posteriori_values,
    "visited_perms" = visited_perms,
    "start_perm" = permutations::id,
    "last_perm" = NULL,
    "last_perm_log_posteriori" = NULL,
    "iterations_performed" = iterations_to_perform,
    "optimization_algorithm_used" = "brute_force",
    "post_probabilities" = probabilities,
    "did_converge" = TRUE,
    "best_perm_log_posteriori" = log_posteriori_values[which.max(log_posteriori_values)],
    "optimization_time" = NA,
    "whole_optimization_time" = NA,
    "all_n0" = NA
  )


  new_gipsmult(
    list(best_perm), Ss, numbers_of_observations,
    delta, D_matrices,
    was_mean_estimated = FALSE, optimization_info
  ) # was_mean_estimated will be changed in the `find_MAP` function
}



#' Combining 2 gips objects
#'
#' g2 was optimized with a single optimization method. g1 was potentially non-optimized or optimized once, or optimized multiple times.
#' If g2 was optimized with "brute_force", forget the g1.
#'
#' @noRd
combine_gips <- function(g1, g2, show_progress_bar = FALSE) {
  # first, adjust the number of observations:
  attr(g2, "numbers_of_observations") <- attr(g1, "numbers_of_observations")
  attr(g2, "was_mean_estimated") <- attr(g1, "was_mean_estimated")

  if (is.null(attr(g1, "optimization_info")) ||
    attr(g2, "optimization_info")[["optimization_algorithm_used"]] == "brute_force") { # when brute_force was used, forget the initial optimization

    attr(g2, "optimization_info")[["original_perm"]] <- g1[[1]]

    return(g2)
  }

  # g1 is also an effect of optimization.
  optimization_info1 <- attr(g1, "optimization_info")
  optimization_info2 <- attr(g2, "optimization_info")

  n1 <- length(optimization_info1[["log_posteriori_values"]])
  n2 <- length(optimization_info2[["log_posteriori_values"]])

  if (all(is.na(optimization_info1[["visited_perms"]])) || all(is.na(optimization_info2[["visited_perms"]]))) {
    if (!all(is.na(optimization_info1[["visited_perms"]])) || !all(is.na(optimization_info2[["visited_perms"]]))) {
      rlang::warn("You wanted to save visited_perms on one of the optimized `gips` objects but forget it for the other. This is not possible, so both will be forgotten.")
      optimization_info2[["post_probabilities"]] <- NULL
      optimization_info1[["post_probabilities"]] <- NULL
    }
    visited_perms <- NA
  } else {
    visited_perms <- c(optimization_info1[["visited_perms"]], optimization_info2[["visited_perms"]]) # WoW, one can use `c()` to combine lists!
  }
  optimization_algorithm_used <- c(optimization_info1[["optimization_algorithm_used"]], optimization_info2[["optimization_algorithm_used"]])

  if (all(optimization_algorithm_used == "Metropolis_Hastings") &&
    !is.null(optimization_info2[["post_probabilities"]])) {
    post_probabilities <- estimate_probabilities(visited_perms, show_progress_bar) # TODO(This can be combined more optimally when !is.null(optimization_info1[["post_probabilities"]]). It is significant, because those calculations are like the same speed as the MH itself. However, I (Adam) think this will be rarely done nevertheless.)
  } else {
    post_probabilities <- NULL
  }

  optimization_info_new <- list(
    "original_perm" = optimization_info1[["original_perm"]],
    "acceptance_rate" = (n1 * optimization_info1[["acceptance_rate"]] + n2 * optimization_info2[["acceptance_rate"]]) / (n1 + n2),
    "log_posteriori_values" = c(optimization_info1[["log_posteriori_values"]], optimization_info2[["log_posteriori_values"]]),
    "visited_perms" = visited_perms,
    "start_perm" = optimization_info1[["start_perm"]],
    "last_perm" = optimization_info2[["last_perm"]],
    "last_perm_log_posteriori" = optimization_info2[["last_perm_log_posteriori"]],
    "iterations_performed" = c(optimization_info1[["iterations_performed"]], optimization_info2[["iterations_performed"]]),
    "optimization_algorithm_used" = optimization_algorithm_used,
    "post_probabilities" = post_probabilities,
    "did_converge" = optimization_info2[["did_converge"]],
    "best_perm_log_posteriori" = max(optimization_info1[["best_perm_log_posteriori"]], optimization_info2[["best_perm_log_posteriori"]]),
    "optimization_time" = c(optimization_info1[["optimization_time"]], optimization_info2[["optimization_time"]]),
    "whole_optimization_time" = optimization_info1[["whole_optimization_time"]] + optimization_info2[["whole_optimization_time"]],
    "all_n0" = c(optimization_info1[["all_n0"]], optimization_info2[["all_n0"]])
  )

  if (optimization_info1[["best_perm_log_posteriori"]] > optimization_info2[["best_perm_log_posteriori"]]) {
    g_out <- g1 # in the continuation the new best was not found
  } else {
    g_out <- g2
  }

  attr(g_out, "optimization_info") <- optimization_info_new

  g_out
}