#' Quadratic Discriminant Analysis with multiple gips-projected covariances
#'
#' Quadratic Discriminant Analysis (QDA) in which each class covariance matrix
#' is projected using the \emph{gipsmult} framework, allowing for structured
#' permutation symmetry across multiple covariance matrices.
#'
#' This function is a modification of \code{\link[MASS]{qda}} in which the
#' class-specific covariance matrices are jointly projected to improve
#' numerical stability and exploit shared symmetry assumptions.
#'
#' @name gipsmultqda
#' @aliases
#'   gipsmultqda gipsmultqda.default gipsmultqda.formula
#'   gipsmultqda.data.frame gipsmultqda.matrix
#'   predict.gipsmultqda print.gipsmultqda model.frame.gipsmultqda
#'
#' @usage
#' gipsmultqda(x, ...)
#'
#' \method{gipsmultqda}{formula}(formula, data, ..., subset, na.action)
#'
#' \method{gipsmultqda}{default}(x, grouping, prior = proportions,
#'   nu = 5, MAP = TRUE, optimizer = NULL, max_iter = NULL, ...)
#'
#' \method{gipsmultqda}{data.frame}(x, ...)
#'
#' \method{gipsmultqda}{matrix}(x, grouping, ..., subset, na.action)
#'
#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...}.
#'   The response is the grouping factor and the right-hand side specifies
#'   the (non-factor) discriminators.
#' @param data An optional data frame, list or environment from which variables
#'   specified in \code{formula} are preferentially taken.
#' @param x (required if no formula is given as the principal argument)
#'   a matrix or data frame containing the explanatory variables.
#' @param grouping A factor specifying the class for each observation.
#' @param prior Prior probabilities of class membership. Must sum to one.
#' @param nu Degrees of freedom parameter used internally during covariance
#'   projection.
#' @param MAP Logical; if \code{TRUE}, a maximum a posteriori covariance
#'   projection is used.
#' @param optimizer Character string specifying the optimization method used
#'   for covariance projection. If \code{NULL}, a default choice is made
#'   based on the problem dimension.
#' @param max_iter Maximum number of iterations for stochastic optimizers.
#' @param subset An index vector specifying the cases to be used in the training
#'   sample. (NOTE: must be named.)
#' @param na.action A function specifying the action to be taken if \code{NA}s
#'   are found.
#' @param ... Arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{"gipsmultqda"} containing:
#' \itemize{
#'   \item \code{prior}: prior probabilities of the groups
#'   \item \code{counts}: number of observations per group
#'   \item \code{means}: group means
#'   \item \code{scaling}: array of group-specific scaling matrices derived
#'     from the projected covariance matrices
#'   \item \code{ldet}: log-determinants of the projected covariance matrices
#'   \item \code{lev}: class labels
#'   \item \code{N}: total number of observations
#'   \item \code{optimization_info}: information returned by the covariance
#'     projection optimizer
#'   \item \code{call}: the matched call
#' }
#'
#' @details
#' In contrast to classical QDA, which estimates each class covariance matrix
#' independently, \code{gipsmultqda} performs a joint projection of all class
#' covariance matrices using the \emph{gipsmult} framework. This allows the
#' incorporation of shared permutation symmetries and can improve classification
#' performance in high-dimensional or small-sample regimes.
#'
#' Several classification rules are available via
#' \code{\link{predict.gipsmultqda}}, including plug-in, predictive, debiased,
#' and leave-one-out cross-validation.
#'
#' @note
#' This function is not a drop-in replacement for \code{\link[MASS]{qda}}.
#' The covariance estimation, returned object, and classification rules
#' differ substantially.
#'
#' The theoretical background and details of the covariance projection are
#' documented in the \code{gipsmult} package.
#'
#' @seealso
#' \code{\link[MASS]{qda}}, \code{\link{predict.gipsmultqda}},
#' \code{\link{gipsqda}}, \code{\link{gipslda}}
#'
#' @examples
#' tr <- sample(1:50, 25)
#' train <- rbind(iris3[tr, , 1], iris3[tr, , 2], iris3[tr, , 3])
#' test <- rbind(iris3[-tr, , 1], iris3[-tr, , 2], iris3[-tr, , 3])
#' cl <- factor(c(rep("s", 25), rep("c", 25), rep("v", 25)))
#' z <- gipsmultqda(train, cl)
#' predict(z, test)$class
#'
#' @keywords multivariate classification
#'
#' @export
gipsmultqda <- function(x, ...) UseMethod("gipsmultqda")

#' @exportS3Method
gipsmultqda.formula <- function(formula, data, ..., subset, na.action) {
  m <- match.call(expand.dots = FALSE)
  m$... <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  grouping <- model.response(m)
  x <- model.matrix(Terms, m)
  xvars <- as.character(attr(Terms, "variables"))[-1L]
  if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)
  if (xint > 0) x <- x[, -xint, drop = FALSE]
  res <- gipsmultqda.default(x, grouping, ...)
  res$terms <- Terms
  cl <- match.call()
  cl[[1L]] <- as.name("gipsmultqda")
  res$call <- cl
  res$contrasts <- attr(x, "contrasts")
  res$xlevels <- .getXlevels(Terms, m)
  res$na.action <- attr(m, "na.action")
  res
}

#' @exportS3Method
gipsmultqda.data.frame <- function(x, ...) {
  res <- gipsmultqda(structure(data.matrix(x), class = "matrix"), ...)
  cl <- match.call()
  cl[[1L]] <- as.name("gipsmultqda")
  res$call <- cl
  res
}

#' @exportS3Method
gipsmultqda.matrix <- function(x, grouping, ..., subset, na.action) {
  if (!missing(subset)) {
    x <- x[subset, , drop = FALSE]
    grouping <- grouping[subset]
  }
  if (!missing(na.action)) {
    dfr <- na.action(structure(list(g = grouping, x = x),
      class = "data.frame"
    ))
    grouping <- dfr$g
    x <- dfr$x
  }
  #    res <- NextMethod("gipsmultqda")
  res <- gipsmultqda.default(x, grouping, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("gipsmultqda")
  res$call <- cl
  res
}

#' @exportS3Method
gipsmultqda.default <-
  function(x, grouping, prior = proportions, nu = 5, MAP = TRUE, optimizer = NULL, max_iter = NULL, ...) {
    if (is.null(dim(x))) stop("'x' is not a matrix")
    x <- as.matrix(x)
    if (any(!is.finite(x))) {
      stop("infinite, NA or NaN values in 'x'")
    }
    n <- nrow(x)
    p <- ncol(x)
    if (n != length(grouping)) {
      stop("nrow(x) and length(grouping) are different")
    }
    g <- as.factor(grouping)
    lev <- levels(g)
    counts <- as.vector(table(g))
    names(counts) <- lev
    # if(any(counts < p+1)) stop("some group is too small for 'gipsmultqda'")
    proportions <- counts / length(g)
    ng <- length(proportions)
    # allow for supplied prior
    if (any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if (length(prior) != ng) stop("'prior' is of incorrect length")
    names(prior) <- lev

    if (is.null(optimizer)) {
      if (p < 10) {
        optimizer <- "BF"
      } else {
        optimizer <- "MH"
      }
    }
    if (optimizer == "MH" & is.null(max_iter)) {
      max_iter <- 100
      warning("MH optimizer set but 'max_iter' argument is unspecified \n
        Setting max_iter = 100 by default")
    }

    # means by group (rows) and variable (columns)
    group.means <- tapply(x, list(rep(g, ncol(x)), col(x)), mean)
    scaling <- array(dim = c(p, p, ng))
    ldet <- numeric(ng)
    ####################################################################################
    cXs <- list()

    for (i in 1L:ng) {
      x_i <- (x[unclass(g) == i, ])
      cXs[[i]] <- cov(x_i)
      group.means[i, ] <- colMeans(x_i)
    }

    # for (i in 1L:ng){
    #     cX <- MASS::cov.mve(x[unclass(g) == i, ])
    #     cXs[[i]] <- cX$cov
    #     group.means[i,] <- cX$center
    # }
    pr_cov_opt_info <- project_covs(cXs, counts, MAP, optimizer, max_iter)
    cov_proj <- pr_cov_opt_info$covs
    cov_proj <- lapply(cov_proj, function(mat) desingularize(mat, 0.05))
    optimization_info <- pr_cov_opt_info$opt_info
    for (i in 1L:ng) {
      cov_proj_sng <- cov_proj[[i]]
      sX <- svd(cov_proj_sng, nu = 0)
      scaling[, , i] <- sX$v %*% diag(sqrt(1 / sX$d), , p)
      ldet[i] <- sum(log(sX$d))
    }
    ####################################################################################
    if (is.null(dimnames(x))) {
      dimnames(scaling) <- list(NULL, as.character(1L:p), lev)
    } else {
      dimnames(scaling) <- list(colnames(x), as.character(1L:p), lev)
      dimnames(group.means)[[2L]] <- colnames(x)
    }
    cl <- match.call()
    cl[[1L]] <- as.name("gipsmultqda")
    res <- list(
      prior = prior, counts = counts, means = group.means,
      scaling = scaling, ldet = ldet, lev = lev, N = n, call = cl, optimization_info = optimization_info
    )
    class(res) <- "gipsmultqda"
    res
  }

#' @exportS3Method
predict.gipsmultqda <- function(object, newdata, prior = object$prior,
                                method = c(
                                  "plug-in", "predictive", "debiased",
                                  "looCV"
                                ), ...) {
  if (!inherits(object, "gipsmultqda")) stop("object not of class \"gipsmultqda\"")
  method <- match.arg(method)
  if (method == "looCV" && !missing(newdata)) {
    stop("cannot have leave-one-out CV with 'newdata'")
  }
  if (is.null(mt <- object$call$method)) mt <- "moment"
  if (method == "looCV" && !(mt == "moment" || mt == "mle")) {
    stop(gettext("cannot use leave-one-out CV with method %s",
      sQuote(mt),
      domain = NA
    ))
  }
  ngroup <- length(object$prior)
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if (length(prior) != ngroup) stop("'prior' is of incorrect length")
  }
  if (!is.null(Terms <- object$terms)) {
    # formula fit
    if (missing(newdata)) {
      newdata <- model.frame(object)
    } else {
      newdata <- model.frame(as.formula(delete.response(Terms)),
        newdata,
        na.action = function(x) x,
        xlev = object$xlevels
      )
    }
    x <- model.matrix(delete.response(Terms), newdata,
      contrasts.arg = object$contrasts
    )
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    if (xint > 0) x <- x[, -xint, drop = FALSE]
    if (method == "looCV") g <- model.response(newdata)
  } else { #
    # matrix or data-frame fit
    if (missing(newdata)) {
      if (!is.null(sub <- object$call$subset)) {
        newdata <-
          eval.parent(parse(text = paste(
            deparse(object$call$x,
              backtick = TRUE
            ),
            "[", deparse(sub, backtick = TRUE), ",]"
          )))
        g <- eval.parent(parse(text = paste(
          deparse(object$call[[3L]],
            backtick = TRUE
          ),
          "[", deparse(sub, backtick = TRUE), "]"
        )))
      } else {
        newdata <- eval.parent(object$call$x)
        g <- eval.parent(object$call[[3L]])
      }
      if (!is.null(nas <- object$call$na.action)) {
        df <- data.frame(g = g, X = newdata)
        df <- eval(call(nas, df))
        g <- df$g
        newdata <- df$X
      }
      g <- as.factor(g)
    }
    if (is.null(dim(newdata))) {
      dim(newdata) <- c(1, length(newdata))
    } # a row vector
    x <- as.matrix(newdata) # to cope with dataframes
  }
  p <- ncol(object$means)
  if (ncol(x) != p) stop("wrong number of variables")
  if (length(colnames(x)) > 0L &&
    any(colnames(x) != dimnames(object$means)[[2L]])) {
    warning("variable names in 'newdata' do not match those in 'object'")
  }
  dist <- matrix(0, nrow = nrow(x), ncol = ngroup)
  if (method == "plug-in") {
    for (i in 1L:ngroup) {
      dev <- ((x - matrix(object$means[i, ], nrow(x),
        ncol(x),
        byrow = TRUE
      )) %*% object$scaling[, , i])
      dist[, i] <- 0.5 * rowSums(dev^2) + 0.5 * object$ldet[i] - log(prior[i])
    }
    #        dist <- exp( -(dist - min(dist, na.rm=T)))
    dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
  } else if (method == "looCV") {
    n <- nrow(x)
    NG <- 1
    if (mt == "mle") NG <- 0
    ldet <- matrix(0, n, ngroup)
    for (i in 1L:ngroup) {
      dev <- ((x - matrix(object$means[i, ], nrow(x), p, byrow = TRUE))
      %*% object$scaling[, , i])
      dist[, i] <- rowSums(dev^2)
      ldet[, i] <- object$ldet[i]
    }
    nc <- object$counts[g]
    ind <- cbind(1L:n, g)
    fac <- 1 - nc / (nc - 1) / (nc - NG) * dist[ind]
    fac[] <- pmax(fac, 1e-10) # possibly degenerate dsn
    ldet[ind] <- log(fac) + p * log((nc - NG) / (nc - 1 - NG)) + ldet[ind]
    dist[ind] <- dist[ind] * (nc^2 / (nc - 1)^2) * (nc - 1 - NG) / (nc - NG) / fac
    dist <- 0.5 * dist + 0.5 * ldet -
      matrix(log(prior), n, ngroup, byrow = TRUE)
    dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
  } else if (method == "debiased") {
    for (i in 1L:ngroup) {
      nk <- object$counts[i]
      Bm <- p * log((nk - 1) / 2) - sum(digamma(0.5 * (nk - 1L:ngroup)))
      dev <- ((x - matrix(object$means[i, ],
        nrow = nrow(x),
        ncol = ncol(x), byrow = TRUE
      )) %*% object$scaling[, , i])
      dist[, i] <- 0.5 * (1 - (p - 1) / (nk - 1)) * rowSums(dev^2) +
        0.5 * object$ldet[i] - log(prior[i]) + 0.5 * Bm - p / (2 * nk)
    }
    dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
  } else {
    for (i in 1L:ngroup) {
      nk <- object$counts[i]
      dev <- ((x - matrix(object$means[i, ],
        nrow = nrow(x),
        ncol = ncol(x), byrow = TRUE
      ))
      %*% object$scaling[, , i])
      dev <- 1 + rowSums(dev^2) / (nk + 1)
      dist[, i] <- prior[i] * exp(-object$ldet[i] / 2) *
        dev^(-nk / 2) * (1 + nk)^(-p / 2)
    }
  }
  posterior <- dist / drop(dist %*% rep(1, ngroup))
  cl <- factor(max.col(posterior),
    levels = seq_along(object$lev),
    labels = object$lev
  )
  dimnames(posterior) <- list(rownames(x), object$lev)
  list(class = cl, posterior = posterior)
}

#' @exportS3Method
print.gipsmultqda <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    names(cl)[2L] <- ""
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\nPrior probabilities of groups:\n")
  print(x$prior, ...)
  cat("\nGroup means:\n")
  print(x$means, ...)
  cat("\nPermutations with their estimated probabilities:\n")
  print(x$optimization_info)
  invisible(x)
}

#' @exportS3Method
model.frame.gipsmultqda <- model.frame.gipslda
