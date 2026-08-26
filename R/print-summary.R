#' Print a fitted gipsDA model
#'
#' @param x A fitted gipsDA model.
#' @param ... Further arguments passed to printing methods.
#'
#' @return Invisibly returns `x`.
#' @keywords internal
#' @exportS3Method
print.gipslda <- function(x, ...) {
  .print_gipsda_model(x, model_name = "gipslda", ...)
  cat("\nCoefficients of linear discriminants:\n")
  print(x$scaling, ...)

  if (!is.null(x$svd)) {
    svd <- x$svd
    if (!is.null(dimnames(x$scaling))) {
      names(svd) <- dimnames(x$scaling)[[2L]]
    }
    if (length(svd) > 1L) {
      cat("\nProportion of trace:\n")
      print(round(svd^2 / sum(svd^2), 4L), ...)
    }
  }

  invisible(x)
}

#' @rdname print.gipslda
#' @exportS3Method
print.gipsqda <- function(x, ...) {
  .print_gipsda_model(x, model_name = "gipsqda", ...)

  if (!is.null(x$ldet)) {
    cat("\nLog determinants of projected covariance matrices:\n")
    print(x$ldet, ...)
  }

  invisible(x)
}

#' @rdname print.gipslda
#' @exportS3Method
print.gipsmultqda <- function(x, ...) {
  .print_gipsda_model(x, model_name = "gipsmultqda", ...)

  if (!is.null(x$ldet)) {
    cat("\nLog determinants of projected covariance matrices:\n")
    print(x$ldet, ...)
  }

  invisible(x)
}

.print_gipsda_model <- function(x, model_name, ...) {
  if (!is.null(cl <- x$call)) {
    names(cl)[2L] <- ""
    cat("Call:\n")
    dput(cl, control = NULL)
  }

  cat("\nModel:", model_name, "\n")
  cat("Number of observations:", x$N, "\n")
  cat("Number of groups:", length(x$prior), "\n")
  cat("Number of predictors:", ncol(x$means), "\n")

  if (!is.null(x$fit_info)) {
    cat("\nFitting options:\n")
    print(x$fit_info, ...)
  }

  cat("\nPrior probabilities of groups:\n")
  print(x$prior, ...)

  cat("\nClass counts:\n")
  print(x$counts, ...)

  cat("\nGroup means:\n")
  print(x$means, ...)

  cat("\nPermutations with their estimated probabilities:\n")
  .print_optimization_info(x$optimization_info, ...)
}

#' Summarize a fitted gipsDA model
#'
#' @param object A fitted gipsDA model.
#' @param ... Further arguments passed to or from methods.
#'
#' @return An object of class `"summary.gipsda"`.
#' @keywords internal
#' @exportS3Method
summary.gipslda <- function(object, ...) {
  structure(
    list(
      model = "gipslda",
      call = object$call,
      n = object$N,
      p = ncol(object$means),
      groups = object$lev,
      counts = object$counts,
      prior = object$prior,
      means = object$means,
      fit_info = object$fit_info,
      scaling = object$scaling,
      svd = object$svd,
      proportion_trace = .gipslda_trace_proportion(object),
      optimization_info = object$optimization_info
    ),
    class = c("summary.gipslda", "summary.gipsda")
  )
}

#' @rdname summary.gipslda
#' @exportS3Method
summary.gipsqda <- function(object, ...) {
  structure(
    list(
      model = "gipsqda",
      call = object$call,
      n = object$N,
      p = ncol(object$means),
      groups = object$lev,
      counts = object$counts,
      prior = object$prior,
      means = object$means,
      fit_info = object$fit_info,
      scaling_dim = dim(object$scaling),
      ldet = object$ldet,
      optimization_info = object$optimization_info
    ),
    class = c("summary.gipsqda", "summary.gipsda")
  )
}

#' @rdname summary.gipslda
#' @exportS3Method
summary.gipsmultqda <- function(object, ...) {
  structure(
    list(
      model = "gipsmultqda",
      call = object$call,
      n = object$N,
      p = ncol(object$means),
      groups = object$lev,
      counts = object$counts,
      prior = object$prior,
      means = object$means,
      fit_info = object$fit_info,
      scaling_dim = dim(object$scaling),
      ldet = object$ldet,
      optimization_info = object$optimization_info
    ),
    class = c("summary.gipsmultqda", "summary.gipsda")
  )
}

.gipslda_trace_proportion <- function(object) {
  if (is.null(object$svd) || length(object$svd) == 0L) {
    return(NULL)
  }

  prop <- object$svd^2 / sum(object$svd^2)

  if (!is.null(dimnames(object$scaling))) {
    names(prop) <- dimnames(object$scaling)[[2L]]
  }

  prop
}

#' Print a gipsDA model summary
#'
#' @param x A summary object produced by `summary.gipslda()`,
#'   `summary.gipsqda()`, or `summary.gipsmultqda()`.
#' @param ... Further arguments passed to printing methods.
#'
#' @return Invisibly returns `x`.
#' @keywords internal
#' @exportS3Method
print.summary.gipsda <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    names(cl)[2L] <- ""
    cat("Call:\n")
    dput(cl, control = NULL)
  }

  cat("\nModel:", x$model, "\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of groups:", length(x$groups), "\n")
  cat("Number of predictors:", x$p, "\n")

  if (!is.null(x$fit_info)) {
    cat("\nFitting options:\n")
    print(x$fit_info, ...)
  }

  cat("\nClass counts:\n")
  print(x$counts, ...)

  cat("\nPrior probabilities of groups:\n")
  print(x$prior, ...)

  cat("\nGroup means:\n")
  print(x$means, ...)

  if (!is.null(x$proportion_trace)) {
    cat("\nProportion of trace:\n")
    print(round(x$proportion_trace, 4L), ...)
  }

  if (!is.null(x$ldet)) {
    cat("\nLog determinants of projected covariance matrices:\n")
    print(x$ldet, ...)
  }

  if (!is.null(x$scaling_dim)) {
    cat("\nScaling array dimensions:\n")
    print(x$scaling_dim, ...)
  }

  cat("\nPermutation optimization information:\n")
  print(x$optimization_info, ...)

  invisible(x)
}

.print_optimization_info <- function(x, ...) {
  if (is.null(x)) {
    cat("NULL\n")
    return(invisible(NULL))
  }

  if (is.list(x) && !inherits(x, "data.frame") && !inherits(x, "gips")) {
    for (nm in names(x)) {
      cat("\nGroup:", nm, "\n")
      print(x[[nm]], ...)
    }
  } else {
    print(x, ...)
  }

  invisible(NULL)
}
