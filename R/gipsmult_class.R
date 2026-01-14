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
  if (inherits(perm, "gipsmult")) {
      perm <- perm[[1]]
    }
  if (!inherits(perm, c("gips_perm", "permutation"))) {
    perm <- permutations::permutation(perm)
  }
  if (inherits(perm, "gips_perm")) {
    gips_perm_object <- perm
  } else {
    gips_perm_object <- gips::gips_perm(perm, nrow(Ss[[1]]))
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

#' @describeIn gipsmult Constructor. It is only intended for low-level use.
#'
#' @param list_of_gips_perm A list with a single element of
#'     a `gips_perm` class. The base object for the `gipsmult` object.
#' @param optimization_info For internal use only. `NULL` or the list with
#'     information about the optimization process.
#'
#' @returns `new_gipsmult()` returns an object of
#'     a `gipsmult` class without the safety checks.
#'
#' @export
new_gipsmult <- function(list_of_gips_perm, Ss, numbers_of_observations,
                     delta, D_matrices, was_mean_estimated, optimization_info) {
  if (!is.list(list_of_gips_perm) ||
    !inherits(list_of_gips_perm[[1]], "gips_perm") ||
    !list_of_matrices_check(Ss) ||
    !list_of_matrices_check(D_matrices) ||
    !noo_check(numbers_of_observations) ||
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




#' Printing `gipsmult` object
#'
#' Printing function for a `gipsmult` class.
#'
#' @param x An object of a `gipsmult` class.
#' @param digits The number of digits after the comma
#'     for a posteriori to be presented. It can be negative.
#'     By default, `Inf`. It is passed to [base::round()].
#' @param compare_to_original A logical. Whether to print how many
#'     times more likely is the current permutation compared to:
#' * the identity permutation `()` (for unoptimized `gipsmult` object);
#' * the starting permutation (for optimized `gipsmult` object).
#' @param log_value A logical. Whether to print the logarithmic value.
#'     Default to `FALSE`.
#' @param oneline A logical. Whether to print in
#'     one or multiple lines. Default to `FALSE`.
#' @param ... The additional arguments passed to [base::cat()].
#'
#' @seealso
#' * [find_MAP()] - The function that makes
#'     an optimized `gipsmult` object out of the unoptimized one.
#'
#' @returns Returns an invisible `NULL`.
#' @export
#'
#' @examples
#' Ss <- list(matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE), matrix(c(2, 1, 3, 7), nrow = 2, byrow = TRUE))
#' noo <- c(10,13)
#' g <- gipsmult(S, noo, perm = "(12)")
#' print(g, digits = 4, oneline = TRUE)
print.gipsmult <- function(x, digits = 3, compare_to_original = TRUE,
                       log_value = FALSE, oneline = FALSE, ...) {

  printing_text <- paste0("The permutation ", as.character(x[[1]]))

  if (is.null(attr(x, "optimization_info"))) { # it is unoptimized gipsmult object
    log_posteriori <- log_posteriori_of_gipsmult(x)
    if (is.nan(log_posteriori) || is.infinite(log_posteriori)) {
      # See ISSUE#5; We hope the implementation of log calculations have stopped this problem.
      rlang::warn(c("gipsmult is yet unable to process these Ss matrices, and produced a NaN or Inf value while trying.",
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
  } else { # it is optimized gipsmult object
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
      rlang::warn(c("gipsmult is yet unable to process these Ss matrices, and produced a NaN or Inf value while trying.",
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

#' Plot optimized matrix or optimization `gipsmult` object
#'
#' Plot heatmaps of the MAP covariance matrices estimator
#' or the convergence of the optimization method.
#' The plot depends on the `type` argument.
#'
#' @param x Object of a `gipsmult` class.
#' @param type A character vector of length 1. One of
#'     `c("heatmap", "MLE", "best", "all", "both", "n0", "block_heatmap")`:
#'   * `"heatmap"`, `"MLE"` - Plots a heatmap of the Maximum Likelihood
#'       Estimator of the covariance matrix given the permutation.
#'       That is, the `S` matrix inside the `gipsmult` object
#'       projected on the permutation in the `gipsmult` object.
#'   * `"best"` - Plots the line of the biggest a posteriori found over time.
#'   * `"all"` - Plots the line of a posteriori for all visited states.
#'   * `"both"` - Plots both lines from "all" and "best".
#'   * `"n0"` - Plots the line of `n0`s that were spotted during optimization
#'       (only for "MH" optimization).
#'   * `"block_heatmap"` - Plots a heatmap of diagonally block representation of `S`.
#'       Non-block entries (equal to 0) are white for better clarity.
#'
#' The default value is `NA`, which will be changed to "heatmap" for
#'     non-optimized `gipsmult` objects and to "both" for optimized ones.
#'     Using the default produces a warning.
#'     All other arguments are ignored for
#'     the `type = "heatmap"`, `type = "MLE"`, or `type = "block_heatmap"`.
#' @param logarithmic_y,logarithmic_x A boolean.
#'     Sets the axis of the plot in logarithmic scale.
#' @param color Vector of colors to be used to plot lines.
#' @param title_text Text to be in the title of the plot.
#' @param xlabel Text to be on the bottom of the plot.
#' @param ylabel Text to be on the left of the plot.
#' @param show_legend A boolean. Whether or not to show a legend.
#' @param ylim Limits of the y axis. When `NULL`,
#'     the minimum, and maximum of the [log_posteriori_of_gipsmult()] are taken.
#' @param xlim Limits of the x axis. When `NULL`,
#'     the whole optimization process is shown.
#' @param ... Additional arguments passed to
#'     other various elements of the plot.
#'
#' @returns When `type` is one of `"best"`, `"all"`, `"both"` or `"n0"`,
#'     returns an invisible `NULL`.
#'     When `type` is one of `"heatmap"`, `"MLE"` or `"block_heatmap"`,
#'     returns an object of class `ggplot`.
#'
#' @seealso
#' * [find_MAP()] - Usually, the `plot.gipsmult()`
#'     is called on the output of `find_MAP()`.
#' * [gipsmult()] - The constructor of a `gipsmult` class.
#'     The `gipsmult` object is used as the `x` parameter.
#'
#' @export
#'
#' @examples
#' require("MASS") # for mvrnorm()
#'
#' perm_size <- 6
#' mu1 <- runif(6, -10, 10)
#' mu2 <- runif(6, -10, 10) # Assume we don't know the means
#' sigma1 <- matrix(
#'   data = c(
#'     1.0, 0.8, 0.6, 0.4, 0.6, 0.8,
#'     0.8, 1.0, 0.8, 0.6, 0.4, 0.6,
#'     0.6, 0.8, 1.0, 0.8, 0.6, 0.4,
#'     0.4, 0.6, 0.8, 1.0, 0.8, 0.6,
#'     0.6, 0.4, 0.6, 0.8, 1.0, 0.8,
#'     0.8, 0.6, 0.4, 0.6, 0.8, 1.0
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' )
#' sigma2 <- matrix(
#'   data = c(
#'     1.0, 0.5, 0.2, 0.0, 0.2, 0.5,
#'     0.5, 1.0, 0.5, 0.2, 0.0, 0.2,
#'     0.2, 0.5, 1.0, 0.5, 0.2, 0.0,
#'     0.0, 0.2, 0.5, 1.0, 0.5, 0.2,
#'     0.2, 0.0, 0.2, 0.5, 1.0, 0.5,
#'     0.5, 0.2, 0.0, 0.2, 0.5, 1.0
#'   ),
#'   nrow = perm_size, byrow = TRUE
#' )
#' # sigma1 and sigma2 are matrices invariant under permutation (1,2,3,4,5,6)
#' numbers_of_observations <- c(21,37)
#' Z1 <- MASS::mvrnorm(numbers_of_observations[1], mu = mu1, Sigma = sigma1)
#' Z2 <- MASS::mvrnorm(numbers_of_observations[2], mu = mu2, Sigma = sigma2)
#' S1 <- cov(Z1)
#' S2 <- cov(Z2) # Assume we have to estimate the mean
#'
#' g <- gipsmult(list(S1,S2), numbers_of_observations)
#' if (require("graphics")) {
#'   plot(g, type = "MLE")
#' }
#'
#' g_map <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "hill_climbing")
#' if (require("graphics")) {
#'   plot(g_map, type = "both", logarithmic_x = TRUE)
#' }
#'
#' if (require("graphics")) {
#'   plot(g_map, type = "MLE")
#' }
#' # Now, the output is (most likely) different because the permutation
#'   # `g_map[[1]]` is (most likely) not an identity permutation.
#'
#' g_map_MH <- find_MAP(g, max_iter = 30, show_progress_bar = FALSE, optimizer = "MH")
#' if (require("graphics")) {
#'   plot(g_map_MH, type = "n0")
#' }
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

    rlang::inform(c("You used the default value of the 'type' argument in `plot()` for gipsmult object.",
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
        "i" = "For non-optimized `gipsmult` objects only the `type = 'heatmap', 'MLE' or 'block_heatmap'` can be used.",
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
      reason = "to use `plot.gipsmult()` with `type %in% c('heatmap', 'MLE', 'block_heatmap')`; without those packages, the `stats::heatmap()` will be used"
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


#' Extract probabilities for `gipsmult` object optimized with `return_probabilities = TRUE`
#'
#' After the `gipsmult` object was optimized with
#' the `find_MAP(return_probabilities = TRUE)` function, then
#' those calculated probabilities can be extracted with this function.
#'
#' @param g An object of class `gipsmult`.
#'     A result of a `find_MAP(return_probabilities = TRUE)`.
#'
#' @returns Returns a numeric vector, calculated values of probabilities.
#' Names contain permutations this probabilities represent.
#' For `gipsmult` object optimized with `find_MAP(return_probabilities = FALSE)`,
#' it returns a `NULL` object.
#' It is sorted according to the probability.
#'
#' @export
#'
#' @seealso
#' * [find_MAP()] - The `get_probabilities_from_gipsmult()`
#'     is called on the output of
#'     `find_MAP(return_probabilities = TRUE, save_all_perms = TRUE)`.
#'
#' @examples
#' Ss <- list(matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE), matrix(c(2, 1, 3, 7), nrow = 2, byrow = TRUE))
#' noo <- c(10,13)
#' g <- gipsmult(Ss, noo)
#' g_map <- find_MAP(g,
#'   optimizer = "BF", show_progress_bar = FALSE,
#'   return_probabilities = TRUE, save_all_perms = TRUE
#' )
#'
#' get_probabilities_from_gipsmult(g_map)
get_probabilities_from_gipsmult <- function(g) {

    if (is.null(attr(g, "optimization_info"))) {
    rlang::abort(c("There was a problem identified with the provided arguments:",
      "i" = "`gipsmult` objects has to be optimized with `find_MAP(return_probabilities=TRUE)` to use `get_probabilities_from_gipsmult()` function.",
      "x" = "You did not optimized `g`.",
      "i" = "Did You use the wrong `g` as an argument for this function?",
      "i" = "Did You forget to optimize `g`?"
    ))
  }

  if (is.null(attr(g, "optimization_info")[["post_probabilities"]])) {
    rlang::inform(c(
      "You called `get_probabilities_from_gipsmult(g)` on the `gipsmult` object that does not have saved probabilities.",
      "x" = "`NULL` will be returned",
      "i" = "Did You use the wrong `g` as an argument for this function?",
      "i" = "Did You forget to optimize with `g <- find_MAP(return_probabilities = TRUE)`?",
      "i" = "Did You unintentionally use `g <- forget_perms(g)`?"
    ))
  }

  attr(g, "optimization_info")[["post_probabilities"]]
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





