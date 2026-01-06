#' The constructor of a `gipsmult` class.
#'
#' Create a `gipsmult` object.
#' This object will contain initial data and all other information
#' needed to find the most likely invariant permutation.
#' It will not perform optimization. One must call
#' the [find_MAP()] function to do it. See the examples below.
#'
#' @param S A matrix; empirical covariance matrix.
#'     When `Z` is the observed data:
#' * if one does not know the theoretical mean and has to
#'     estimate it with the observed mean, use `S = cov(Z)`,
#'     and leave parameter `was_mean_estimated = TRUE` as default;
#' * if one know the theoretical mean is 0, use
#'     `S = (t(Z) %*% Z) / number_of_observations`, and set
#'     parameter `was_mean_estimated = FALSE`.
#' @param number_of_observations A number of data points
#'     that `S` is based on.
#' @param delta A number, hyper-parameter of a Bayesian model.
#'     It has to be strictly bigger than 1.
#'     See the **Hyperparameters** section below.
#' @param D_matrix Symmetric, positive-definite matrix of the same size as `S`.
#'     Hyper-parameter of a Bayesian model.
#'     When `NULL`, the (hopefully) reasonable one is derived from the data.
#'     For more details, see the **Hyperparameters** section below.
#' @param was_mean_estimated A boolean.
#' * Set `TRUE` (default) when your `S` parameter is a result of
#'     a [stats::cov()] function.
#' * Set FALSE when your `S` parameter is a result of
#'     a `(t(Z) %*% Z) / number_of_observations` calculation.
#' @param perm An optional permutation to be the base for the `gipsmult` object.
#'     It can be of a `gips_perm` or a `permutation` class, or anything
#'     the function [permutations::permutation()] can handle.
#'     It can also be of a `gipsmult` class, but
#'     it will be interpreted as the underlying `gips_perm`.
#'
#' @section Methods for a `gipsmult` class:
#' * [plot.gipsmult()]
#' * [print.gipsmult()]
#'
#' @section Hyperparameters:
#' We encourage the user to try `D_matrix = d * I`, where `I` is an identity matrix of a size
#' `p x p` and `d > 0` for some different `d`.
#' When `d` is small compared to the data (e.g., `d = 0.1 * mean(diag(S))`),
#'     bigger structures will be found.
#' When `d` is big compared to the data (e.g., `d = 100 * mean(diag(S))`),
#'     the posterior distribution does not depend on the data.
#'
#' Taking `D_matrix = d * I` is equivalent to setting \code{S <- S / d}.
#'
#' The default for `D_matrix` is `D_matrix = d * I`, where
#' `d = mean(diag(S))`, which is equivalent to modifying `S`
#' so that the mean value on the diagonal is 1.
#'
#' In the Bayesian model, the prior distribution for
#' the covariance matrix is a generalized case of
#' [Wishart distribution](https://en.wikipedia.org/wiki/Wishart_distribution).
#'
#' @returns `gipsmult()` returns an object of
#'     a `gipsmult` class after the safety checks.
#'
#' @export
#' @seealso
#' * [stats::cov()] – The `Ss` parameter, as a list of empirical covariance matrices,
#'     is most of the time a result of the `cov()` function.
#'     For more information, see
#'     [Wikipedia - Estimation of covariance matrices](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices).
#' * [find_MAP()] – The function that finds
#'     the Maximum A Posteriori (MAP) Estimator
#'     for a given `gipsmult` object.
#' * [gips_perm()] – The constructor of a `gips_perm` class.
#'     The `gips_perm` object is used as the base object for
#'     the `gipsmult` object.
#'
#' @examples
#' perm_size <- 5
#' numbers_of_observations <- c(15, 18, 19)
#' Sigma <- diag(rep(1, perm_size))
#' n_matrices <- 3
#' df <- 20
#' Ss <- rWishart(n = n_matrices, df = df, Sigma = Sigma)
#' Ss <- lapply(1:n_matrices, function(x) Ss[, , x])
#' g <- gipsmult(Ss, numbers_of_observations)
#'
#' g_map <- find_MAP(g, show_progress_bar = FALSE, optimizer = "brute_force")
#' g_map
#'
#' print(g_map)
#'
#' if (require("graphics")) {
#'   plot(g_map, type = "MLE", logarithmic_x = TRUE)
#' }


gipsmult <- function(Ss, numbers_of_observations, delta = 3, D_matrices = NULL,
                     was_mean_estimated = TRUE, perm = "") {
  if (inherits(perm, "gips")) {
      gips::validate_gips(perm)
      perm <- perm[[1]]
    }
  if (!inherits(perm, c("gips_perm", "permutation"))) {
    perm <- permutations::permutation(perm)
  }
  if (inherits(perm, "gips_perm")) {
    gips_perm_object <- perm # it is already a `gips_perm`
  } else {
    gips_perm_object <- gips::gips_perm(perm, nrow(Ss[[1]])) # it is of a `cycle` class from permutations package (it was checked in `check_correctness_of_arguments()`. Make it 'gips_perm' class
  }


  if (is.null(D_matrices)) {
    D_matrices <- lapply(Ss, function(y) diag(x = mean(diag(y)), nrow = ncol(y)))

  }

  new_gipsmult(
    list(gips_perm_object),
    Ss, numbers_of_observations, delta = delta,
    D_matrices = D_matrices,
    was_mean_estimated = was_mean_estimated,
    optimization_info = NULL
  )


}

list_of_matrices_check <- function(x) {
  if(!is.list(x)) {
    return(FALSE)
  }
  if (!is.matrix(x[[1]])) {
    return(FALSE)
  }
  shap <- dim(x[[1]])
  for(y in x[-1]) {
    if(!is.matrix(y) | !all(dim(y) == shap)) {
      return(FALSE)
    }
  }
  return(TRUE)
}

SDN_compatibility_check <- function(Ss, D_matrices, numbers_of_observations) {
  n <- length(numbers_of_observations)
  if (length(Ss) != n | length(D_matrices) != n) {
    return(FALSE)
  }
  if(!all(dim(D_matrices[[1]]) == dim(Ss[[1]]))) {
    return(FALSE)
  }
  return(TRUE)
}

new_gipsmult <- function(list_of_gips_perm, Ss, numbers_of_observations,
                     delta, D_matrices, was_mean_estimated, optimization_info) {
  if (!is.list(list_of_gips_perm) ||
    !inherits(list_of_gips_perm[[1]], "gips_perm") ||
    !list_of_matrices_check(Ss) ||
    !list_of_matrices_check(D_matrices) ||
    !is.numeric(numbers_of_observations) ||
    !SDN_compatibility_check(Ss, D_matrices, numbers_of_observations) ||
    !is.numeric(delta) ||
    !is.logical(was_mean_estimated) ||
    !(is.null(optimization_info) || is.list(optimization_info))) {
    rlang::abort(c("x" = "`gipsmult` object cannot be created from those arguments."))
  }


  structure(list_of_gips_perm,
            Ss = Ss,
            numbers_of_observations = numbers_of_observations,
            delta = delta,
            D_matrices = D_matrices,
            was_mean_estimated = was_mean_estimated,
            optimization_info = optimization_info,
            class = c("gipsmult")
  )
}

validate_gipsmult <- function(g) {
  if (!(inherits(g, "gipsmult"))) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "`g` must be of a `gipsmult` class.",
      "x" = paste0(
        "You provided `g` with `class(g) == (",
        paste(class(g), collapse = ", "), ")`."
      )
    ))
  }

  if (!(length(g) == 1)) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `length(g)` must be `1`.",
      "x" = paste0(
        "You provided `g` with `length(g) == ",
        length(g), "`."
      )
    ))
  }
  if (!is.list(g)) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `g` must be a list.",
      "x" = paste0(
        "You provided `g` with `typeof(g) == '",
        typeof(g), "'."
      )
    ))
  }

  perm <- g[[1]]
  Ss <- attr(g, "Ss")
  numbers_of_observations <- attr(g, "numbers_of_observations")
  delta <- attr(g, "delta")
  D_matrices <- attr(g, "D_matrices")
  was_mean_estimated <- attr(g, "was_mean_estimated")
  optimization_info <- attr(g, "optimization_info")

  if (!inherits(perm, "gips_perm")) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `g[[1]]` must be an object of a `gips_perm` class.",
      "x" = paste0(
        "You provided `g[[1]]` with `class(g[[1]]) == (",
        paste(class(perm), collapse = ", "),
        ")`."
      )
    ))
  } else {
    tryCatch(
      {
        gips::validate_gips_perm(perm)
      },
      error = function(cond) {
        rlang::abort(c("There was a problem identified with provided argument:",
          "i" = "The `g[[1]]` must be an object of a `gips_perm` class.",
          "x" = paste0(
            "You provided `g[[1]]` with `class(g[[1]]) == 'gips_perm'`, but your g[[1]] does not pass `validate_gips_perm(g[[1]])`."
          )
        ))
      }
    )
  }

  check_correctness_of_arguments( # max_iter, return_probabilities and show_progress_bar are to be checked here, but some value has to be passed
    Ss = Ss, numbers_of_observations = numbers_of_observations,
    max_iter = 2, start_perm = perm,
    delta = delta, D_matrices = D_matrices, was_mean_estimated = was_mean_estimated,
    return_probabilities = FALSE, save_all_perms = TRUE, show_progress_bar = FALSE
  )

  if (!(is.null(optimization_info) || is.list(optimization_info))) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `optimization_info` value must be either a `NULL`, or a list.",
      "x" = paste0(
        "You provided `attr(g, 'optimization_info')` with type ",
        typeof(optimization_info), ""
      )
    ))
  }

  if (is.list(optimization_info)) { # Validate the `optimization_info` after the optimization
    legal_fields <- c("original_perm", "acceptance_rate", "log_posteriori_values", "visited_perms", "start_perm", "last_perm", "last_perm_log_posteriori", "iterations_performed", "optimization_algorithm_used", "post_probabilities", "did_converge", "best_perm_log_posteriori", "optimization_time", "whole_optimization_time", "all_n0")

    lacking_fields <- setdiff(legal_fields, names(optimization_info))
    illegal_fields <- setdiff(names(optimization_info), legal_fields)

    abort_text <- character(0)

    if (!(length(lacking_fields) == 0)) {
      abort_text <- c("x" = paste0(
        "Your `attr(g, 'optimization_info')` lacks the following fields: ",
        paste(lacking_fields, collapse = ", "), ""
      ))
    }
    if (!(length(illegal_fields) == 0)) {
      abort_text <- c(abort_text,
        "x" = paste0(
          "Your `attr(g, 'optimization_info')` has the following, unexpected fields: ",
          paste(illegal_fields, collapse = ", "), ""
        )
      )
    }

    # abort the validation
    if (length(abort_text) > 0) {
      rlang::abort(c("There was a problem with the 'optimization_info' attribute.",
        "i" = paste0(
          "After optimization, `attr(g, 'optimization_info')` must be a list of ",
          length(legal_fields), " elements with names: ",
          paste(legal_fields, collapse = ", "), ""
        ),
        "x" = paste0("You have a list of ", length(names(optimization_info)), " elements."),
        abort_text,
        "i" = "Did You accidentally edit `attr(g, 'optimization_info')` by yourself?",
        "i" = "Did You accidentally set one of `attr(g, 'optimization_info')` elements to `NULL` or `NA`?"
      ))
    }

    # All the fields as named as they should be. Check if their content are as expected:
    abort_text <- character(0)
    additional_info <- 0 # for calculation of the number of problems

    # original_perm is not validated :<

    if (!((is.numeric(optimization_info[["acceptance_rate"]]) &&
      (length(optimization_info[["acceptance_rate"]]) == 1) &&
      optimization_info[["acceptance_rate"]] >= 0 &&
      optimization_info[["acceptance_rate"]] <= 1) ||
      optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] == "brute_force")) { # when brute_force, acceptance_rate is NULL
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['acceptance_rate']]` must be a number in range [0, 1].",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['acceptance_rate']] == (",
          paste(optimization_info[["acceptance_rate"]], collapse = ", "),
          ")`."
        )
      )
    }
    if (!(optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force" ||
      is.null(optimization_info[["acceptance_rate"]]))) {
      abort_text <- c(abort_text,
        "i" = "When brute force algorithm was used for optimization, `attr(g, 'optimization_info')[['acceptance_rate']]` must be a `NULL`.",
        "x" = paste0(
          "You have used brute force algorithm, but `attr(g, 'optimization_info')[['acceptance_rate']] == (",
          paste(optimization_info[["acceptance_rate"]], collapse = ", "),
          ")`."
        )
      )
    }
    if (!(is.numeric(optimization_info[["log_posteriori_values"]]))) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['log_posteriori_values']]` must be a vector of numbers.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['log_posteriori_values']] == ",
          typeof(optimization_info[["log_posteriori_values"]]),
          "`."
        )
      )
    }
    if (optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force") {
      if (!(all(is.na(optimization_info[["visited_perms"]])) || (is.list(optimization_info[["visited_perms"]])))) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['visited_perms']]` must be a list or `NA`.",
          "x" = paste0(
            "You have `attr(g, 'optimization_info')[['visited_perms']]` of type ",
            typeof(optimization_info[["visited_perms"]]),
            ""
          )
        )
      } else if (length(optimization_info[["visited_perms"]]) == 0) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['visited_perms']]` must be a list with some elements or an `NA`.",
          "x" = paste0(
            "Your `attr(g, 'optimization_info')[['visited_perms']]` is a list, but of a length 0."
          )
        )
      } else if (!(all(is.na(optimization_info[["visited_perms"]])) || (inherits(optimization_info[["visited_perms"]][[1]], "gips_perm")))) { # It only checks for the first one, because checking for every would be too expensive
        abort_text <- c(abort_text,
          "i" = "Elements of `attr(g, 'optimization_info')[['visited_perms']]` must be of a `gips_perm` class.",
          "x" = paste0(
            "You have `class(attr(g, 'optimization_info')[['visited_perms']][[1]]) == (",
            paste(class(optimization_info[["visited_perms"]][[1]]), collapse = ", "),
            ")`."
          )
        )
      } else if (!(all(is.na(optimization_info[["visited_perms"]])) || (identical(optimization_info[["last_perm"]], optimization_info[["visited_perms"]][[length(optimization_info[["visited_perms"]])]])))) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['last_perm']]` must be the last element of `attr(g, 'optimization_info')[['visited_perms']]` list.",
          "x" = paste0("You have `attr(g, 'optimization_info')[['last_perm']]` different from `attr(g, 'optimization_info')[['visited_perms']][[length(attr(g, 'optimization_info')[['visited_perms']])]]`.")
        )
      }

      if (inherits(optimization_info[["last_perm"]], "gips_perm")) {
        abort_text <- tryCatch(
          {
            gips::validate_gips_perm(optimization_info[["last_perm"]])

            # optimization_info[["last_perm"]] is proper gips_perm object
            last_perm_gips <- gipsmult(Ss, numbers_of_observations,
              delta = delta, D_matrices = D_matrices,
              was_mean_estimated = was_mean_estimated,
              perm = optimization_info[["last_perm"]]
            )

            if (!(abs(optimization_info[["last_perm_log_posteriori"]] - log_posteriori_of_gipsmult(last_perm_gips)) < 0.00000001)) {
              abort_text <- c(abort_text,
                "i" = "`attr(g, 'optimization_info')[['last_perm_log_posteriori']]` must be the log_posteriori of `optimization_info[['last_perm']]`.",
                "x" = paste0(
                  "You have `attr(g, 'optimization_info')[['last_perm_log_posteriori']] == ",
                  optimization_info[["last_perm_log_posteriori"]],
                  "`, but `log_posteriori_of_gips(gips(attr(g, 'Ss'), attr(g, 'numbers_of_observations'), delta=attr(g, 'delta'), D_matrices=attr(g, 'D_matrices'), was_mean_estimated=attr(g, 'was_mean_estimated'), perm=attr(g, 'optimization_info')[['last_perm']])) == ",
                  log_posteriori_of_gipsmult(last_perm_gips), "`."
                )
              )
            }

            abort_text # if optimization_info[["last_perm"]] passes the validation, return original text
          },
          error = function(cond) { # this error can only be thrown in `validate_gips_perm()`, so add an appropriate note to the `abort_text`:
            c(abort_text,
              "i" = "The `attr(g, 'optimization_info')[['last_perm']]` must be an object of a `gips_perm` class.",
              "x" = paste0(
                "You provided `attr(g, 'optimization_info')[['last_perm']]` with `class(attr(g, 'optimization_info')[['last_perm']]) == 'gips_perm'`, but your attr(g, 'optimization_info')[['last_perm']] does not pass `validate_gips_perm(attr(g, 'optimization_info')[['last_perm']])`."
              )
            )
          }
        )
      } else {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['last_perm']]` must be an object of class 'gips_perm'.",
          "x" = paste0(
            "You have `attr(g, 'optimization_info')[['last_perm']]` of class ('",
            paste(class(optimization_info[["last_perm"]]), collapse = "', '"), "')."
          )
        )
      }
    } else { # for brute_force, the visited_perms are of class "cycle"
      if (!(all(is.na(optimization_info[["visited_perms"]])) || (is.list(optimization_info[["visited_perms"]])))) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['visited_perms']]` must be a list.",
          "x" = paste0(
            "You have `attr(g, 'optimization_info')[['visited_perms']]` of type ",
            typeof(optimization_info[["visited_perms"]]),
            ""
          )
        )
      } else if (length(optimization_info[["visited_perms"]]) == 0) {
        abort_text <- c(abort_text,
          "i" = "`attr(g, 'optimization_info')[['visited_perms']]` must be a list with some elements.",
          "x" = paste0(
            "Your `attr(g, 'optimization_info')[['visited_perms']]` is a list, but of a length 0."
          )
        )
      } else if (!(all(is.na(optimization_info[["visited_perms"]])) || (inherits(optimization_info[["visited_perms"]][[1]], "list")))) { # It only checks for the first one, because checking for every would be too expensive
        abort_text <- c(abort_text,
          "i" = "After optimization with brute force algorithm, elements of `attr(g, 'optimization_info')[['visited_perms']]` must be of a `list` class.",
          "x" = paste0(
            "You have `class(attr(g, 'optimization_info')[['visited_perms']][[1]]) == (",
            paste(class(optimization_info[["visited_perms"]][[1]]), collapse = ", "),
            ")`."
          )
        )
      } else if (!(is.null(optimization_info[["last_perm"]]))) {
        abort_text <- c(abort_text,
          "i" = "After optimization with brute force algorithm, `attr(g, 'optimization_info')[['last_perm']]` must be a `NULL`.",
          "x" = paste0("You have `attr(g, 'optimization_info')[['last_perm']]` of type ", typeof(optimization_info[["last_perm"]]), "")
        )
      }
    }

    if (!(all(gips:::is.wholenumber(optimization_info[["iterations_performed"]])))) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['iterations_performed']]` must be a vector of whole numbers.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['iterations_performed']] == (",
          paste(optimization_info[["iterations_performed"]], collapse = ", "),
          ")`."
        )
      )
    } else if (!(sum(optimization_info[["iterations_performed"]]) <= length(optimization_info[["log_posteriori_values"]]))) {
      abort_text <- c(abort_text,
        "i" = "In every iteration at least one value of log_posteriori is calculated.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['iterations_performed']] == ",
          sum(optimization_info[["iterations_performed"]]),
          "`, which is more than `length(attr(g, 'optimization_info')[['log_posteriori_values']]) == ",
          length(optimization_info[["log_posteriori_values"]]), "`."
        )
      )
    }
    if (!all(optimization_info[["optimization_algorithm_used"]] %in% c("Metropolis_Hastings", "hill_climbing", "brute_force"))) { # Even if MH was used, it would produce the text "Metropolis_Hastings"
      abort_text <- c(abort_text,
        "i" = "The available optimization algorithms are 'Metropolis_Hastings', 'hill_climbing' and 'brute_force'.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_algorithm_used']] == (",
          paste(optimization_info[["optimization_algorithm_used"]], collapse = ", "),
          ")`."
        )
      )
    } else if ((all(optimization_info[["optimization_algorithm_used"]] != "Metropolis_Hastings") && # all optimization algorithms
      optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] != "brute_force") && # last optimization algorithm
      !is.null(optimization_info[["post_probabilities"]])) {
      abort_text <- c(abort_text,
        "i" = "`post_probabilities` can only be obtained with 'Metropolis_Hastings' or 'brute_force' optimization method.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_algorithm_used']] == c('",
          paste0(optimization_info[["optimization_algorithm_used"]], collapse = "', '"),
          "')` and the `attr(g, 'optimization_info')[['post_probabilities']]` is not `NULL`, but is of type `",
          typeof(optimization_info[["post_probabilities"]]), "`."
        )
      )
    } else if ((length(optimization_info[["visited_perms"]]) > 0) &&
      !(all(is.na(optimization_info[["visited_perms"]])) || is.null(optimization_info[["post_probabilities"]]) ||
        length(optimization_info[["post_probabilities"]]) <= length(optimization_info[["visited_perms"]]))) {
      abort_text <- c(abort_text,
        "i" = "Every element of `attr(g, 'optimization_info')[['post_probabilities']]` was taken from a visited permutation, so it is in `attr(g, 'optimization_info')[['visited_perms']]`.",
        "x" = paste0(
          "You have `length(attr(g, 'optimization_info')[['visited_perms']]) == ",
          length(optimization_info[["post_probabilities"]]),
          "`, but `length(attr(g, 'optimization_info')[['post_probabilities']]) == ",
          length(optimization_info[["visited_perms"]]),
          "` which are not equal."
        )
      )
    } else if (!(is.null(optimization_info[["post_probabilities"]]) ||
      (all(optimization_info[["post_probabilities"]] <= 1) &&
        all(optimization_info[["post_probabilities"]] >= 0) && # it should be >0, but it sometimes underflow to 0
        (sum(optimization_info[["post_probabilities"]]) < 1.001) && # Allow small error
        (sum(optimization_info[["post_probabilities"]]) > 0.999)))) {
      abort_text <- c(abort_text,
        "i" = "The vector of `attr(g, 'optimization_info')[['post_probabilities']]` must have properties of probability. All elements in range [0, 1] and sums to 1.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['post_probabilities']]` in a range [",
          min(optimization_info[["post_probabilities"]]), ",",
          max(optimization_info[["post_probabilities"]]),
          "] and with the sum ",
          sum(optimization_info[["post_probabilities"]]),
          ""
        )
      )
    }
    if ((!(optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] %in% c("hill_climbing", "brute_force"))) && # The last optimization_algorithm_used has to be hill_climbing or brute_force to make the convergence
      !is.null(optimization_info[["did_converge"]])) {
      abort_text <- c(abort_text,
        "i" = "`did_converge` can only be obtained with 'hill_climbing' or 'brute_force' optimization method.",
        "x" = paste0(
          "The last optimization method You used was `attr(g, 'optimization_info')[['optimization_algorithm_used']][length(attr(g, 'optimization_info')[['optimization_algorithm_used']])] == ",
          optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])],
          "` and the `attr(g, 'optimization_info')[['did_converge']]` is not `NULL`, but is of type ",
          typeof(optimization_info[["did_converge"]]), ""
        )
      )
    } else if ((optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] == "hill_climbing") &&
      !is.logical(optimization_info[["did_converge"]])) {
      abort_text <- c(abort_text,
        "i" = "When 'hill_climbing' optimization method, the `did_converge` must be `TRUE` or `FALSE`.",
        "x" = paste0(
          "The last optimization method You used was `attr(g, 'optimization_info')[['optimization_algorithm_used']][length(attr(g, 'optimization_info')[['optimization_algorithm_used']])] == ",
          optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])],
          "` and the `attr(g, 'optimization_info')[['did_converge']]` is not of type logical, but it is of type ",
          typeof(optimization_info[["did_converge"]]), ""
        )
      )
    } else if ((optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])] == "hill_climbing") &&
      is.na(optimization_info[["did_converge"]])) {
      abort_text <- c(abort_text,
        "i" = "When 'hill_climbing' optimization method, the `did_converge` must be `TRUE` or `FALSE`.",
        "x" = paste0(
          "The last optimization method You used was `attr(g, 'optimization_info')[['optimization_algorithm_used']][length(optimization_info[['optimization_algorithm_used']])] == ",
          optimization_info[["optimization_algorithm_used"]][length(optimization_info[["optimization_algorithm_used"]])],
          "` and the `attr(g, 'optimization_info')[['did_converge']]` is of type logical, but it is a NA."
        )
      )
    }
    best_perm_gips <- gipsmult(Ss, numbers_of_observations, delta = delta, D_matrices = D_matrices, was_mean_estimated = was_mean_estimated, perm = perm) # this perm is g[[1]]
    if (!(abs(optimization_info[["best_perm_log_posteriori"]] - log_posteriori_of_gipsmult(best_perm_gips)) < 0.00000001)) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['best_perm_log_posteriori']]` must be the log_posteriori of the base object, `g[[1]]`.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['best_perm_log_posteriori']] == ",
          optimization_info[["best_perm_log_posteriori"]],
          "`, but `log_posteriori_of_gips(gips(attr(g, 'Ss'), attr(g, 'numbers_of_observations'), delta=attr(g, 'delta'), D_matrices=attr(g, 'D_matrices'), was_mean_estimated=attr(g, 'was_mean_estimated'), perm=g[[1]])) == ",
          log_posteriori_of_gipsmult(best_perm_gips), "`."
        )
      )
    }
    if (any(is.na(optimization_info[["optimization_time"]]))) {
      additional_info <- additional_info + 2
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['optimization_time']]` is initially set to `NA`, but that state of the gips object should not be available to the user.",
        "x" = "You have `is.na(attr(g, 'optimization_info')[['optimization_time']]) == TRUE`.",
        "i" = "Did You use the inner optimizers like `gips:::Metropolis_Hastings()` or `gips:::hill_climbing()` in stead of the exported function `gips::find_MAP()`?",
        "i" = "Did You modify the `find_MAP()` function?"
      )
    } else if (!inherits(optimization_info[["optimization_time"]], "difftime")) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['optimization_time']]` has to be of a class 'difftime'.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_time']]` of a class (",
          paste0(class(optimization_info[["optimization_time"]]), collapse = ", "), ")."
        )
      )
    } else if (any(optimization_info[["optimization_time"]] < 0)) { # allow underflow of time float to 0
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['optimization_time']]` has to be a non negative time difference.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['optimization_time']] == ",
          optimization_info[["optimization_time"]], "`."
        )
      )
    }
    if (is.na(optimization_info[["whole_optimization_time"]])) {
      additional_info <- additional_info + 2
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['whole_optimization_time']]` is initially set to `NA`, but that state of the gips object should not be available to the user.",
        "x" = "You have `is.na(attr(g, 'optimization_info')[['whole_optimization_time']]) == TRUE`.",
        "i" = "Did You use the inner optimizers like `gips:::Metropolis_Hastings()` or `gips:::hill_climbing()` in stead of the exported function `gips::find_MAP()`?",
        "i" = "Did You modify the `find_MAP()` function?"
      )
    } else if (!inherits(optimization_info[["whole_optimization_time"]], "difftime")) {
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['whole_optimization_time']]` has to be of a class 'difftime'.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['whole_optimization_time']]` of a class (",
          paste0(class(optimization_info[["whole_optimization_time"]]), collapse = ", "), ")."
        )
      )
    } else if (optimization_info[["whole_optimization_time"]] < 0) { # allow underflow of time float to 0
      abort_text <- c(abort_text,
        "i" = "`attr(g, 'optimization_info')[['whole_optimization_time']]` has to be a non negative time difference.",
        "x" = paste0(
          "You have `attr(g, 'optimization_info')[['whole_optimization_time']] == ",
          optimization_info[["whole_optimization_time"]], "`."
        )
      )
    }

    # TODO(Validate that there is the same number of algorithms used and length of other things dependent on it)


    if (length(abort_text) > 0) {
      if (!gips:::is.wholenumber((length(abort_text) - additional_info) / 2)) {
        rlang::inform(paste0(
          "You found a small bug in gips package. We calculated there was ",
          (length(abort_text) - additional_info) / 2,
          " problems, but it is not a whole number. Please inform us about that bug by opening an ISSUE on https://github.com/PrzeChoj/gips/issues"
        ))
      }
      abort_text <- c(
        paste0(
          "There were ", (length(abort_text) - additional_info) / 2,
          " problems identified with `attr(g, 'optimization_info')`:"
        ),
        abort_text
      )

      if (length(abort_text) > 11) {
        abort_text <- c(abort_text[1:11],
          "x" = paste0("... and ", (length(abort_text) - 1) / 2 - 5, " more problems")
        )
      }

      abort_text <- c(abort_text,
        "i" = "Did You accidentally edit `attr(g, 'optimization_info')` by yourself?",
        ">" = "If You think You've found a bug in a package, please open an ISSUE on https://github.com/PrzeChoj/gips/issues"
      )

      rlang::abort(abort_text)
    }
  }

  g
}

check_correctness_of_arguments <- function(Ss, numbers_of_observations, max_iter,
                                           start_perm, delta, D_matrices, was_mean_estimated,
                                           return_probabilities, save_all_perms, show_progress_bar) {
  if (!is.list(Ss)) {
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `Ss` must be a list.",
      "x" = paste0(
        "You provided `Ss` with `typeof(Ss) == '",
        typeof(Ss), "'."
      )
    ))
  }
  for (S in Ss){
    s_check(S)
    if (ncol(S) != ncol(Ss[[1]])){
      rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `Ss` must be a list of matrices with the same shape",
      "x" = paste0(
        "There are at least 2 different shapes in the matrices you provided"
      )
    ))
    }
  }

  if (!is.numeric(numbers_of_observations)){
    rlang::abort(c("There was a problem identified with provided argument:",
      "i" = "The `numbers_of_observation` must be a numeric type.",
      "x" = paste0(
        "You provided `numbers_of_observations` with `typeof(numbers_of_observations) == '",
        typeof(numbers_of_observations), "' which is not numeric."
      )
    ))
  }
  if(!exists('abort_text')){
    abort_text <- character(0)
  }
  additional_info <- 0
  for (i in numbers_of_observations){
    noo_check(i)
  }
  if (!(is.infinite(max_iter) || gips:::is.wholenumber(max_iter))) {
    abort_text <- c(abort_text,
      "i" = "`max_iter` must be either infinite (for hill_climbing optimizer) or a whole number.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    )
  } else if (max_iter < 2) {
    abort_text <- c(abort_text,
      "i" = "`max_iter` must be at least 2.",
      "x" = paste0("You provided `max_iter == ", max_iter, "`.")
    )
  }
  if (!(permutations::is.cycle(start_perm) || inherits(start_perm, "gips_perm"))) {
    abort_text <- c(abort_text,
      "i" = "`start_perm` must be the output of `gips_perm()` function, or of a `cycle` class form `permutations` package.", # this is not true, but it is close enough
      "x" = paste0(
        "You provided `start_perm` with `class(start_perm) == (",
        paste(class(start_perm), collapse = ", "),
        ")`."
      )
    )
  } else if (!(permutations::is.cycle(start_perm) || attr(start_perm, "size") == ncol(Ss[[1]]))) {
    abort_text <- c(abort_text,
      "i" = "`start_perm` must have the `size` attribute equal to the shape of a square matrices `Ss`",
      "x" = paste0(
        "You provided `start_perm` with `size == ",
        attr(start_perm, "size"),
        "`, but the `Ss` matrices You provided has ",
        ncol(Ss[[1]]), " columns."
      )
    )
  }
  if (is.null(delta)) {
    abort_text <- c(abort_text,
      "i" = "`delta` must not be `NULL`.",
      "x" = "Your provided `delta` is a `NULL`."
    )
  } else if (delta <= 1) { # See documentation of internal function `G_function()` in `calculate_gamma_function.R`
    abort_text <- c(abort_text,
      "i" = "`delta` must be strictly bigger than 1.",
      "x" = paste0("You provided `delta == ", delta, "`.")
    )
  }
  if (!(is.null(D_matrices) || is.list(D_matrices))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrices` must either be `NULL` or a list.",
      "x" = paste0(
        "You provided `D_matrices` with type ",
        typeof(D_matrices), ""
      )
    )
  } else if (!(is.null(D_matrices) || all(vapply(D_matrices, function(x) is.matrix(x), logical(1))))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrices` must either be `NULL` or a list of square matrices with equal dimensions.",
      "x" = "You provided `D_matrices` as a list, but not all its elements are matrices"
    )
  } else if (!(is.null(D_matrices) || all(vapply(D_matrices, function (x) ncol(x) == nrow(x), logical(1))))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrices` must either be `NULL` or a list of square matrices with equal dimensions.",
      "x" = "You provided `D_matrices` as a list of matrices, but there are some that are not square"
    )
  } else if (!(is.null(D_matrices) || all(vapply(D_matrices[-1], function (x) ncol(x) == ncol(D_matrices[[1]]), logical(1))))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrices` must either be `NULL` or a list of square matrices with equal dimensions.",
      "x" = "You provided `D_matrices` as a list of matrices, but they don't have equal dimensions"
    )
  } else if (!(is.null(D_matrices) || all(vapply(Ss, function (x) ncol(x) == ncol(D_matrices[[1]]), logical(1))))) {
    abort_text <- c(abort_text,
      "i" = "`Ss` must be a list of square matrices with the same shape as  `D_matrices`.",
      "x" = "There is a mismatch between dimensions of Ss and D_matrices"
    )
  } else if (any(vapply(D_matrices, function(m) any(is.nan(m)), logical(1)))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrices` must not contain any `NaN`s.",
      "x" = "You provided `D_matrices` with `NaN`s!"
    )
  } else if (any(vapply(D_matrices, function(m) any(is.infinite(m)), logical(1)))) {
    abort_text <- c(abort_text,
      "i" = "`D_matrices` must not contain any infinite values.",
      "x" = "You provided `D_matrices` with infinite values!"
    )
  }
  if (!is.logical(was_mean_estimated)) {
    abort_text <- c(abort_text,
      "i" = "`was_mean_estimated` must be a logic value (`TRUE` or `FALSE`).",
      "x" = paste0(
        "You provided `was_mean_estimated` with type ",
        typeof(was_mean_estimated), ""
      )
    )
  } else if (is.na(was_mean_estimated)) {
    abort_text <- c(abort_text,
      "i" = "`was_mean_estimated` must be a logic value (`TRUE` or `FALSE`).",
      "x" = "You provided `was_mean_estimated` as an `NA`."
    )
  }
  if (!is.logical(return_probabilities)) {
    abort_text <- c(abort_text,
      "i" = "`return_probabilities` must be a logic value (`TRUE` or `FALSE`).",
      "x" = paste0(
        "You provided `return_probabilities` with type ",
        typeof(return_probabilities), ""
      )
    )
  } else if (is.na(return_probabilities)) {
    abort_text <- c(abort_text,
      "i" = "`return_probabilities` must be a logic value (`TRUE` or `FALSE`).",
      "x" = "You provided `return_probabilities` as an `NA`."
    )
  }
  if (!is.logical(save_all_perms)) {
    abort_text <- c(abort_text,
      "i" = "`save_all_perms` must be a logic value (`TRUE` or `FALSE`).",
      "x" = paste0(
        "You provided `save_all_perms` with type ",
        typeof(save_all_perms), ""
      )
    )
  } else if (is.na(save_all_perms)) {
    abort_text <- c(abort_text,
      "i" = "`save_all_perms` must be a logic value (`TRUE` or `FALSE`).",
      "x" = "You provided `save_all_perms` as an `NA`."
    )
  }
  if (!is.logical(show_progress_bar)) {
    abort_text <- c(abort_text,
      "i" = "`show_progress_bar` must be a logic value (`TRUE` or `FALSE`).",
      "x" = paste0(
        "You provided `show_progress_bar` with type ",
        typeof(show_progress_bar), ""
      )
    )
  } else if (is.na(show_progress_bar)) {
    abort_text <- c(abort_text,
      "i" = "`show_progress_bar` must be a logic value (`TRUE` or `FALSE`).",
      "x" = "You provided `show_progress_bar` as an `NA`."
    )
  }
  if (length(abort_text) > 0) {
    abort_text <- c(
      paste0(
        "There were ", (length(abort_text) - additional_info) / 2,
        " problems identified with the provided arguments:"
      ),
      abort_text
    )

    if (length(abort_text) > 11) {
      abort_text <- c(
        abort_text[1:11],
        paste0("... and ", (length(abort_text) - 1) / 2 - 5, " more problems")
      )
    }

    abort_text <- c(
      abort_text,
      ">" = "If You think You've found a bug in a package, please open an ISSUE on https://github.com/PrzeChoj/gips/issues"
    )

    rlang::abort(abort_text)
  }

  if (return_probabilities && !save_all_perms) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "For calculations of probabilities, all perms have to be available after the optimization process.",
      "x" = "You provided `return_probabilities == TRUE` and `save_all_perms == FALSE`!",
      "i" = "Did You want to set `save_all_perms = TRUE`?",
      "i" = "Did You want to set `return_probabilities = FALSE`?",
      "!" = paste0(
        "Remember that setting `return_probabilities == TRUE` can be computationally costly",
        ifelse(show_progress_bar, " and second progress bar will be shown.", "")
      )
    ))
  }
}

s_check <- function(S){
  if (!is.matrix(S)) {
    rlang::abort(c("There was a problem identified with provided S argument:",
      "i" = "`S` must be a matrix.",
      "x" = paste0(
        "You provided `S` with type ",
        typeof(S), ""
      )
    ))
  }
  abort_text <- character(0)
  additional_info <- 0 # for calculation of the number of problems
  if (ncol(S) != nrow(S)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a square matrix.",
      "x" = paste0(
        "You provided `S` as a matrix, but with different sizes: ",
        ncol(S), " and ", nrow(S), ""
      )
    )
  } else if (!is.numeric(S)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a numeric matrix.",
      "x" = paste0(
        "You provided `S` as a matrix, but with non-numeric values. Your provided type ",
        typeof(S), ""
      )
    )
  } else if (!all(abs(S - t(S)) < 0.000001)) { # this would mean the matrix is not symmetric
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be a symmetric matrix.",
      "x" = "You provided `S` as a matrix, but a non-symmetric one.",
      "i" = "Is your matrix approximately symmetric? Maybe try setting `S <- (S+t(S))/2`?"
    )
    additional_info <- additional_info + 1 # for calculation of the number of problems
  } else if (!gips:::is.positive.semi.definite.matrix(S, tolerance = 1e-06)) {
    abort_text <- c(abort_text,
      "i" = "`S` matrix must be positive semi-definite matrix.",
      "x" = "You provided `S` as a symmetric matrix, but a non-positive-semi-definite one."
    )
  }
}
noo_check <- function(number_of_observations){
  if (is.null(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must not be `NULL`.",
      "x" = "Your provided `number_of_observations` is `NULL`."
    )
  } else if (number_of_observations < 1) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must be at least 1.",
      "x" = paste0(
        "You provided `number_of_observations == ",
        number_of_observations, "`."
      )
    )
  } else if (!gips:::is.wholenumber(number_of_observations)) {
    abort_text <- c(abort_text,
      "i" = "`number_of_observations` must be a whole number.",
      "x" = paste0(
        "You provided `number_of_observations == ",
        number_of_observations, "`."
      )
    )
  }
}

#' @exportS3Method
print.gipsmult <- function(x, digits = 3, compare_to_original = TRUE,
                       log_value = FALSE, oneline = FALSE, ...) {
  validate_gipsmult(x)

  printing_text <- paste0("The permutation ", as.character(x[[1]]))

  if (is.null(attr(x, "optimization_info"))) { # it is unoptimized gips object
    log_posteriori <- log_posteriori_of_gipsmult(x)
    if (is.nan(log_posteriori) || is.infinite(log_posteriori)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(log_posteriori), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrices. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))
    }

    if (compare_to_original) {
      x_id <- gipsmult(
        Ss = attr(x, "Ss"),
        numbers_of_observations = attr(x, "numbers_of_observations"),
        delta = attr(x, "delta"), D_matrices = attr(x, "D_matrices"),
        was_mean_estimated = attr(x, "was_mean_estimated"), perm = ""
      )
      log_posteriori_id <- log_posteriori_of_gipsmult(x_id)

      printing_text <- c(
        printing_text,
        paste0(
          "is ", gips:::convert_log_diff_to_str(log_posteriori - log_posteriori_id, digits),
          " times more likely than the () permutation"
        )
      )
    }
  } else { # it is optimized gips object
    log_posteriori <- attr(x, "optimization_info")[["best_perm_log_posteriori"]]
    x_original <- gipsmult(
      Ss = attr(x, "Ss"),
      numbers_of_observations = attr(x, "numbers_of_observations"),
      delta = attr(x, "delta"), D_matrices = attr(x, "D_matrices"),
      was_mean_estimated = attr(x, "was_mean_estimated"), perm = attr(x, "optimization_info")[["original_perm"]]
    )
    log_posteriori_original <- log_posteriori_of_gipsmult(x_original)

    if (is.nan(log_posteriori) || is.infinite(log_posteriori)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::warn(c("gips is yet unable to process this S matrix, and produced a NaN or Inf value while trying.",
        "x" = paste0("The posteriori value of ", ifelse(is.nan(log_posteriori), "NaN", "Inf"), " occured!"),
        "i" = "We think it can only happen for ncol(S) > 500 or for huge D_matrices. If it is not the case for You, please get in touch with us on ISSUE#5."
      ))
    }

    printing_text <- c(printing_text, paste0(
      "was found after ",
      length(attr(x, "optimization_info")[["log_posteriori_values"]]),
      " posteriori calculations"
    ))

    if (compare_to_original) {
      printing_text <- c(printing_text, paste0(
        "is ", gips:::convert_log_diff_to_str(log_posteriori - log_posteriori_original, digits),
        " times more likely than the ",
        as.character(x_original[[1]]), " permutation"
      ))
    }
  }

  if (log_value) {
    printing_text <- c(
      printing_text,
      paste0(
        "has log posteriori ",
        round(log_posteriori, digits = digits)
      )
    )
  }

  # The first line will end with ":", all following lines will end with ";".
  cat(
    paste0(c(
      printing_text[1],
      paste0(printing_text[-1],
        collapse = ifelse(oneline, "; ", ";\n - ")
      )
    ), collapse = ifelse(oneline, ": ", ":\n - ")),
    ".\n",
    sep = "", ...
  )

  invisible(NULL)
}

#' @exportS3Method
plot.gipsmult <- function(x, type = NA,
                      logarithmic_y = TRUE, logarithmic_x = FALSE,
                      color = NULL,
                      title_text = "Convergence plot",
                      xlabel = NULL, ylabel = NULL,
                      show_legend = TRUE,
                      ylim = NULL, xlim = NULL, ...) {
  # checking the correctness of the arguments:
  if (!requireNamespace("graphics", quietly = TRUE)) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "Package 'graphics' must be installed to use this function.",
      "x" = "Package 'graphics' seems to be unavailable."
    ))
  }

  validate_gipsmult(x)

  if (length(type) != 1) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`type` must be an character vector of length 1.",
      "x" = paste0("You provided `type` with length ", length(type), " which is wrong!")
    ))
  }
  if (is.na(type)) {
    type <- ifelse(is.null(attr(x, "optimization_info")),
      "heatmap",
      "both"
    )

    rlang::inform(c("You used the default value of the 'type' argument in `plot()` for gips object.",
      "i" = paste0(
        "The `type = NA` was automatically changed to `type = '",
        type, "'`."
      )
    ))
  }

  if (type == "MLE") {
    type <- "heatmap"
  }

  if (!(type %in% c("heatmap", "block_heatmap", "all", "best", "both", "n0"))) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`type` must be one of: c('heatmap', 'MLE', 'block_heatmap', 'all', 'best', 'both', 'n0').",
      "x" = paste0("You provided `type == ", type, "`."),
      "i" = "Did You misspell the 'type' argument?"
    ))
  }

  if (type != "block_heatmap" && type != "heatmap" &&
    is.null(attr(x, "optimization_info"))) {
    rlang::abort(
      c(
        "There was a problem identified with the provided arguments:",
        "i" = "For non-optimized `gips` objects only the `type = 'heatmap', 'MLE' or 'block_heatmap'` can be used.",
        "x" = paste0(
          "You did not optimized `x` and provided `type = '",
          type, "'`."
        ),
        "i" = paste0(
          "Did You want to call `x <- find_MAP(g)` and then `plot(x, type = '",
          type, "')`?"
        ),
        "i" = "Did You want to use `type = 'heatmap'`?"
      )
    )
  }

  # plotting:
  if (type == "heatmap" || type == "block_heatmap") {
    rlang::check_installed(c("dplyr", "tidyr", "tibble", "ggplot2"),
      reason = "to use `plot.gips()` with `type %in% c('heatmap', 'MLE', 'block_heatmap')`; without those packages, the `stats::heatmap()` will be used"
    )
    if (type == "block_heatmap") {
      my_projected_matrices <- get_diagonalized_matrices_for_heatmap(x)
    } else {
      my_projected_matrices <- mapply(gips:::project_matrix,attr(x, "Ss"), MoreArgs = list(x[[1]]), SIMPLIFY = FALSE)
    }
    n_Sd3 <- as.integer(ceiling(length(my_projected_matrices)/3))
    if (rlang::is_installed(c("dplyr", "tidyr", "tibble", "ggplot2", "patchwork"))) {
      plots <- lapply(my_projected_matrices, function(mat){
        plot_single_gg(mat, x[[1]])
      })
      combined_plots <- patchwork::wrap_plots(plots, ncol=3)
      return(combined_plots)
    } else { # use the basic plot in R, package `graphics`
      par(mfrow = c(n_Sd3, 3))
      for (matrix in my_projected_matrices) {
        plot_single_stats(matrix)
      }
    }
  }
  if (type %in% c("all", "best", "both")) {
    if (is.null(ylabel)) {
      ylabel <- ifelse(logarithmic_y,
        "log posteriori",
        "posteriori"
      )
    }
    if (is.null(xlabel)) {
      xlabel <- ifelse(logarithmic_x,
        "log10 of number of function calls",
        "number of function calls"
      )
    }
    if (is.null(color)) {
      if (type == "both") {
        color <- c("red", "blue")
      } else {
        color <- "red"
      }
    }
    if (logarithmic_y) {
      y_values_from <- attr(x, "optimization_info")[["log_posteriori_values"]] # values of log_posteriori are logarithmic by default
    } else {
      y_values_from <- exp(attr(x, "optimization_info")[["log_posteriori_values"]])
    }

    y_values_max <- cummax(y_values_from)
    y_values_all <- y_values_from

    num_of_steps <- length(y_values_max)

    if (is.null(xlim)) {
      xlim <- c(1, num_of_steps)
    }

    if (is.null(ylim)) {
      ylim_plot <- c(min(y_values_from), y_values_max[num_of_steps])
      if (type == "best") {
        ylim_plot[1] <- y_values_from[1] # for the "best" type this is the smallest point of the graph
      }
    } else {
      ylim_plot <- ylim
    }

    # make the plot stairs-like
    x_points <- c(1, rep(2:num_of_steps, each = 2))

    if (logarithmic_x) {
      x_points <- log10(x_points)
      xlim <- log10(xlim)
    }

    graphics::plot.new()
    graphics::plot.window(xlim, ylim_plot)

    if (type != "best") {
      # make the plot stairs-like
      y_points <- c(
        rep(y_values_all[1:(length(y_values_all) - 1)], each = 2),
        y_values_all[length(y_values_all)]
      )

      graphics::lines.default(x_points, y_points,
        type = "l", lwd = 3,
        col = color[1], # the first color
        ...
      )
    }
    if (type != "all") {
      # make the plot stairs-like
      y_points <- c(
        rep(y_values_max[1:(length(y_values_max) - 1)], each = 2),
        y_values_max[length(y_values_max)]
      )

      graphics::lines.default(x_points, y_points,
        lwd = 3, lty = 1,
        col = color[length(color)], # the last color
        ...
      )
    }

    graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
    graphics::axis(1, ...)
    graphics::axis(2, ...)
    graphics::box(...)

    if (show_legend) {
      if (type == "both") {
        legend_text <- c(
          "All calculated a posteriori",
          "Maximum a posteriori calculated"
        )
        lty <- c(1, 1)
        lwd <- c(3, 3)
      } else if (type == "all") {
        legend_text <- c("All calculated function values")
        lty <- 1
        lwd <- 3
      } else if (type == "best") {
        legend_text <- c("Maximum function values calculated")
        lty <- 1
        lwd <- 3
      }

      graphics::legend("bottomright",
        inset = .002,
        legend = legend_text,
        col = color,
        lty = lty, lwd = lwd,
        cex = 0.7, box.lty = 0
      )
    }
  }
  if (type == "n0") {
    if (is.null(ylabel)) {
      ylabel <- ifelse(logarithmic_y,
        "log n0",
        "n0"
      )
    }
    if (is.null(xlabel)) {
      xlabel <- ifelse(logarithmic_x,
        "log10 of number of function calls",
        "number of function calls"
      )
    }
    if (is.null(color)) {
      color <- "red"
    }

    if (logarithmic_y) {
      y_values <- log(attr(x, "optimization_info")[["all_n0"]])
    } else {
      y_values <- attr(x, "optimization_info")[["all_n0"]]
    }

    num_of_steps <- length(y_values)

    if (is.null(xlim)) {
      xlim <- c(1, num_of_steps)
    }

    if (is.null(ylim)) {
      ylim_plot <- c(0, max(y_values))
    } else {
      ylim_plot <- ylim
    }

    # make the plot stairs-like
    x_points <- c(1, rep(2:num_of_steps, each = 2))

    if (logarithmic_x) {
      x_points <- log10(x_points)
      xlim <- log10(xlim)
    }

    graphics::plot.new()
    graphics::plot.window(xlim, ylim_plot)

    # make the plot stairs-like
    y_points <- c(
      rep(y_values[1:(length(y_values) - 1)], each = 2),
      y_values[length(y_values)]
    )

    graphics::lines.default(x_points, y_points,
      type = "l", lwd = 3,
      col = color[1], # the first color
      ...
    )

    graphics::title(main = title_text, xlab = xlabel, ylab = ylabel, ...)
    graphics::axis(1, ...)
    graphics::axis(2, ...)
    graphics::box(...)

    if (show_legend) {
      legend_text <- "all perms n0"
      lty <- c(1, 1)
      lwd <- c(3, 3)

      graphics::legend("topright",
                       inset = .002,
                       legend = legend_text,
                       col = color,
                       lty = lty, lwd = lwd,
                       cex = 0.7, box.lty = 0
      )
    }
  }


}

plot_single_gg <- function(my_projected_matrix, perm) {
  p <- ncol(my_projected_matrix)

  if (is.null(colnames(my_projected_matrix))) {
    colnames(my_projected_matrix) <- paste0(seq(1, p))
  }
  if (is.null(rownames(my_projected_matrix))) {
    rownames(my_projected_matrix) <- paste0(seq(1, p))
  }

  my_rownames <- rownames(my_projected_matrix)
  my_colnames <- colnames(my_projected_matrix)
  rownames(my_projected_matrix) <- as.character(1:p)
  colnames(my_projected_matrix) <- as.character(1:p)

  # With this line, the R CMD check's "no visible binding for global variable" warning will not occur:
  col_id <- covariance <- row_id <- NULL

  # Life would be easier with pipes (%>%)
  my_transformed_matrix <- tibble::rownames_to_column(
    as.data.frame(my_projected_matrix),
    "row_id"
  )
  my_transformed_matrix <- tidyr::pivot_longer(my_transformed_matrix,
    -c(row_id),
    names_to = "col_id",
    values_to = "covariance"
  )
  my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
    col_id = as.numeric(col_id)
  )
  my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
    row_id = as.numeric(row_id)
  )
  g_plot <- ggplot2::ggplot(
    my_transformed_matrix,
    ggplot2::aes(x = col_id, y = row_id, fill = covariance)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(na.value = "white") +
    ggplot2::scale_x_continuous(breaks = 1:p, labels = my_rownames) +
    ggplot2::scale_y_reverse(breaks = 1:p, labels = my_colnames) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("Estimated covariance matrix\nprojected on permutation ", perm),
      x = "", y = ""
    )

  return(g_plot)
}

plot_single_stats <- function(my_projected_matrix) {
  if (is.null(color)) { # Setting col = NA or col = NULL turns off the whole plot.
        stats::heatmap(my_projected_matrix,
          symm = TRUE,
          Rowv = NA, Colv = NA, ...
        )
      } else {
        stats::heatmap(my_projected_matrix,
          symm = TRUE,
          Rowv = NA, Colv = NA, col = color, ...
        )
      }
}
get_diagonalized_matrices_for_heatmap <- function(x) {
  Ss <- attr(x, "Ss")
  lapply(Ss, gips:::get_diagonalized_matrix_for_heatmap)
}

get_probabilities_from_gipsmult <- function(g) {
  validate_gipsmult(g)

    if (is.null(attr(g, "optimization_info"))) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`gips` objects has to be optimized with `find_MAP(return_probabilities=TRUE)` to use `get_probabilities_from_gips()` function.",
      "x" = "You did not optimized `g`.",
      "i" = "Did You use the wrong `g` as an argument for this function?",
      "i" = "Did You forget to optimize `g`?"
    ))
  }

  if (is.null(attr(g, "optimization_info")[["post_probabilities"]])) {
    rlang::inform(c(
      "You called `get_probabilities_from_gips(g)` on the `gips` object that does not have saved probabilities.",
      "x" = "`NULL` will be returned",
      "i" = "Did You use the wrong `g` as an argument for this function?",
      "i" = "Did You forget to optimize with `g <- find_MAP(return_probabilities = TRUE)`?",
      "i" = "Did You unintentionally use `g <- forget_perms(g)`?"
    ))
  }

  attr(g, "optimization_info")[["post_probabilities"]]
}





