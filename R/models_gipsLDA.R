# file MASS/R/lda.R
# copyright (C) 1994-2023 W. N. Venables and B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
gipslda <- function(x, ...) UseMethod("gipslda")


gipslda.formula <- function(formula, data, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1L]] <- quote(stats::model.frame)
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    if(xint > 0L) x <- x[, -xint, drop = FALSE]
    res <- gipslda.default(x, grouping, ...)
    res$terms <- Terms
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1L]] <- as.name("gipslda")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- .getXlevels(Terms, m)
    res$na.action <- attr(m, "na.action")
    res
}

gipslda.data.frame <- function(x, ...)
{
    res <- gipslda(structure(data.matrix(x), class = "matrix"), ...)
    cl <- match.call()
    cl[[1L]] <- as.name("gipslda")
    res$call <- cl
    res
}


gipslda.matrix <- function(x, grouping, ..., subset, na.action)
{
    if(!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if(!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x),
                                   class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
#    res <- NextMethod("gipslda")
    res <- gipslda.default(x, grouping, ...)
    cl <- match.call()
    cl[[1L]] <- as.name("gipslda")
    res$call <- cl
    res
}
gipslda.default <-
function(x, grouping, prior = proportions, tol = 1.0e-4,
       nu = 5, weighted_avg = FALSE, MAP = TRUE, optimizer = NULL, max_iter = NULL, ...) {
    if(is.null(dim(x))) stop("'x' is not a matrix")
    x <- as.matrix(x)
    if(any(!is.finite(x)))
      stop("infinite, NA or NaN values in 'x'")
    n <- nrow(x)
    p <- ncol(x)
    if(n != length(grouping))
      stop("nrow(x) and length(grouping) are different")
    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))
    if(!missing(prior)) {
      if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
      if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
      prior <- prior[counts > 0L]
    }
    if(any(counts == 0L)) {
      empty <- lev[counts == 0L]
      warning(sprintf(ngettext(length(empty),
                               "group %s is empty",
                               "groups %s are empty"),
                      paste(empty, collapse = " ")), domain = NA)
      lev1 <- lev[counts > 0L]
      g <- factor(g, levels = lev1)
      counts <- as.vector(table(g))
    }
    if(is.null(optimizer)) {
      if (p < 10) {
          optimizer = "BF"
      }
      else {
          optimizer = "MH"
      }
    }
    if (optimizer == "MH" & is.null(max_iter)) {
        max_iter <- 100
        warning("MH optimizer set but 'max_iter' argument is unspecified \n
        Setting max_iter = 100 by default")
    }

    proportions <- counts/n
    ng <- length(proportions)
    names(prior) <- names(counts) <- lev1
    ## drop attributes to avoid e.g. matrix() methods
    group.means <- tapply(c(x), list(rep(g, p), col(x)), mean)
    f1 <- sqrt(diag(var(x - group.means[g,  ])))
    if(any(f1 < tol)) {
      const <- format((1L:p)[f1 < tol])
      stop(sprintf(ngettext(length(const),
                   "variable %s appears to be constant within groups",
                   "variables %s appear to be constant within groups"),
                   paste(const, collapse = " ")),
           domain = NA)
    }
    # scale columns to unit variance before checking for collinearity
    scaling <- diag(1/f1, , p)
    ####################################################################################
    if (weighted_avg) {
    # compute each class covariance matrix S_g
    cov_list <- tapply(
        1:nrow(x),
        g,
        function(idx) cov(x[idx,, drop=FALSE])
    )

    # weighted pooled covariance: S = (1/n) * sum_g n_g * S_g
    cov_adj <- Reduce(`+`, Map(function(Sg, ng) ng * Sg, cov_list, counts)) / n
    }

    else {
      # adjust to "unbiased" scaling of covariance matrix
      cov_raw <- cov((x - group.means[g,  ]) %*% scaling)
      cov_adj <- n/(n - ng) * cov_raw
    }

    #project covariance using gips
    pr_cov_opt_info <- project_covs(list(cov_adj), n, MAP, optimizer, max_iter)
    cov_proj <- pr_cov_opt_info$covs[[1]]
    cov_proj <- shrink_offdiag_until_det(cov_proj, 0.05, 0.95)
    optimization_info <- pr_cov_opt_info$opt_info
    ####################################################################################
    sX <- svd(cov_proj, nu = 0L)
    rank <- sum(sX$d > tol^2)
    if(rank == 0L) stop("rank = 0: variables are numerically constant")
    if(rank < p) warning("variables are collinear")
    scaling <- scaling %*% sX$v[, 1L:rank] %*%
      diag(sqrt(1/sX$d[1L:rank]),,rank)

    xbar <- colSums(prior %*% group.means)
    fac <- 1/(ng - 1)
    X <- sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE) %*% scaling
    X.s <- svd(X, nu = 0L)
    rank <- sum(X.s$d > tol * X.s$d[1L])
    if(rank == 0L) stop("group means are numerically identical")
    scaling <- scaling %*% X.s$v[, 1L:rank]
    if(is.null(dimnames(x)))
      dimnames(scaling) <- list(NULL, paste("LD", 1L:rank, sep = ""))
    else {
      dimnames(scaling) <- list(colnames(x), paste("LD", 1L:rank, sep = ""))
      dimnames(group.means)[[2L]] <- colnames(x)
    }
    cl <- match.call()
    cl[[1L]] <- as.name("gipslda")
    cl$weighted_avg <- NULL
    cl$MAP <- NULL
    cl$optimizer <- NULL
    cl$max_iter <- NULL
    structure(list(prior = prior, counts = counts, means = group.means,
                 scaling = scaling, lev = lev, svd = X.s$d[1L:rank],
                 N = n, optimization_info = optimization_info, call = cl),
            class = "gipslda")
}


predict.gipslda <- function(object, newdata, prior = object$prior, dimen,
			method = c("plug-in", "predictive", "debiased"), ...)
{
    if(!inherits(object, "gipslda")) stop("object not of class \"gipslda\"")
    if(!is.null(Terms <- object$terms)) { # formula fit
        Terms <- delete.response(Terms)
        if(missing(newdata)) newdata <- model.frame(object)
        else {
            newdata <- model.frame(Terms, newdata, na.action=na.pass,
                                   xlev = object$xlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, newdata)
        }
        x <- model.matrix(Terms, newdata, contrasts.arg = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch = 0L)
        if(xint > 0L) x <- x[, -xint, drop = FALSE]
    } else { # matrix or data-frame fit
        if(missing(newdata)) {
            if(!is.null(sub <- object$call$subset))
                newdata <-
                    eval.parent(parse(text = paste(deparse(object$call$x,
                                      backtick = TRUE),
                                      "[", deparse(sub, backtick = TRUE),",]")))
            else newdata <- eval.parent(object$call$x)
            if(!is.null(nas <- object$call$na.action))
                newdata <- eval(call(nas, newdata))
        }
        if(is.null(dim(newdata)))
            dim(newdata) <- c(1L, length(newdata))  # a row vector
        x <- as.matrix(newdata)		# to cope with dataframes
    }

    if(ncol(x) != ncol(object$means)) stop("wrong number of variables")
    if(length(colnames(x)) > 0L &&
      any(colnames(x) != dimnames(object$means)[[2L]]))
         warning("variable names in 'newdata' do not match those in 'object'")
    ng <- length(object$prior)
    if(!missing(prior)) {
        if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
        if(length(prior) != ng) stop("'prior' is of incorrect length")
    }
    ## remove overall means to keep distances small
    means <- colSums(prior*object$means)
    scaling <- object$scaling
    x <- scale(x, center = means, scale = FALSE) %*% scaling
    dm <- scale(object$means, center = means, scale = FALSE) %*% scaling
    method <- match.arg(method)
    dimen <- if(missing(dimen)) length(object$svd) else min(dimen, length(object$svd))
    N <- object$N
    if(method == "plug-in") {
        dm <- dm[, 1L:dimen, drop = FALSE]
        dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x),
                       length(prior), byrow = TRUE) - x[, 1L:dimen, drop=FALSE] %*% t(dm)
        dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
    } else if (method == "debiased") {
        dm <- dm[, 1L:dimen, drop=FALSE]
        dist <- matrix(0.5 * rowSums(dm^2), nrow(x), ng, byrow = TRUE) -
            x[, 1L:dimen, drop=FALSE] %*% t(dm)
        dist <- (N - ng - dimen - 1)/(N - ng) * dist -
            matrix(log(prior) - dimen/object$counts , nrow(x), ng, byrow=TRUE)
        dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
    } else {                            # predictive
        dist <- matrix(0, nrow = nrow(x), ncol = ng)
        p <- ncol(object$means)
        # adjust to ML estimates of covariances
        X <- x * sqrt(N/(N-ng))
        for(i in 1L:ng) {
            nk <- object$counts[i]
            dev <- scale(X, center = dm[i, ], scale = FALSE)
            dev <- 1 + rowSums(dev^2) * nk/(N*(nk+1))
            dist[, i] <- prior[i] * (nk/(nk+1))^(p/2) * dev^(-(N - ng + 1)/2)
        }
    }
    posterior <- dist / drop(dist %*% rep(1, ng))
    nm <- names(object$prior)
    cl <- factor(nm[max.col(posterior)], levels = object$lev)
    dimnames(posterior) <- list(rownames(x), nm)
    list(class = cl, posterior = posterior, x = x[, 1L:dimen, drop = FALSE])
}

print.gipslda <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        names(cl)[2L] <- ""
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nPrior probabilities of groups:\n")
    print(x$prior, ...)
    cat("\nGroup means:\n")
    print(x$means, ...)
    cat("\nCoefficients of linear discriminants:\n")
    print(x$scaling, ...)
    svd <- x$svd
    names(svd) <- dimnames(x$scaling)[[2L]]
    if(length(svd) > 1L) {
        cat("\nProportion of trace:\n")
        print(round(svd^2/sum(svd^2), 4L), ...)
    }
    cat("\nPermutations with their estimated probabilities:\n")
    print(x$optimization_info)
    invisible(x)
}

plot.gipslda <- function(x, panel = panel.gipslda, ..., cex = 0.7,
                     dimen, abbrev = FALSE,
                     xlab = "LD1", ylab = "LD2")
{
    panel.gipslda <- function(x, y, ...) text(x, y, as.character(g), cex = cex, ...)
    if(!is.null(Terms <- x$terms)) { # formula fit
        data <- model.frame(x)
        X <- model.matrix(delete.response(Terms), data)
        g <- model.response(data)
        xint <- match("(Intercept)", colnames(X), nomatch = 0L)
        if(xint > 0L) X <- X[, -xint, drop=FALSE]
    } else { # matrix or data-frame fit
        xname <- x$call$x
        gname <- x$call[[3L]]
        if(!is.null(sub <- x$call$subset)) {
            X <- eval.parent(parse(text=paste(deparse(xname, backtick=TRUE),
                                   "[", deparse(sub, backtick=TRUE),",]")))
            g <- eval.parent(parse(text=paste(deparse(gname, backtick=TRUE),
                                   "[", deparse(sub, backtick=TRUE),"]")))
        } else {
            X <- eval.parent(xname)
            g <- eval.parent(gname)
        }
        if(!is.null(nas <- x$call$na.action)) {
            df <- data.frame(g = g, X = X)
            df <- eval(call(nas, df))
            g <- df$g
            X <- df$X
        }
    }
    if(abbrev) levels(g) <- abbreviate(levels(g), abbrev)
    means <- colMeans(x$means)
    X <- scale(X, center=means, scale=FALSE) %*% x$scaling
    if(!missing(dimen) && dimen < ncol(X)) X <- X[, 1L:dimen, drop = FALSE]
    if(ncol(X) > 2L) {
        pairs(X, panel = panel, ...)
    } else if(ncol(X) == 2L)  {
        eqscplot(X[, 1L:2L], xlab = xlab, ylab = ylab, type = "n", ...)
        panel(X[, 1L], X[, 2L], ...)
    } else ldahist(X[, 1L], g, xlab = xlab, ...)
    invisible(NULL)
}

ldahist <-
function(data, g, nbins = 25, h, x0 = -h/1000, breaks,
	 xlim = range(breaks), ymax = 0, width,
         type = c("histogram", "density", "both"), sep = (type != "density"),
         col = 5L,
	 xlab = deparse(substitute(data)), bty = "n", ...)
{
    xlab
    type <- match.arg(type)
    data <- data[!is.na(data)]
    g <- g[!is.na(data)]
    counts <- table(g)
    groups <- names(counts)[counts > 0L]
    if(missing(breaks)) {
        if(missing(h)) h <- diff(pretty(data, nbins))[1L]
        first <- floor((min(data) - x0)/h)
        last <- ceiling((max(data) - x0)/h)
        breaks <- x0 + h * c(first:last)
    }
    if(type == "histogram" || type == "both") {
        if(any(diff(breaks) <= 0)) stop("'breaks' must be strictly increasing")
        if(min(data) < min(breaks) || max(data) > max(breaks))
            stop("'breaks' do not cover the data")
        est <- vector("list", length(groups))
        names(est) <- groups
        for (grp in groups){
            bin <- cut(data[g == grp], breaks, include.lowest = TRUE)
            est1 <- tabulate(bin, length(levels(bin)))
            est1 <- est1/(diff(breaks) * length(data[g == grp]))
            ymax <- max(ymax, est1)
            est[[grp]] <- est1
        }
    }
    if(type == "density" || type == "both"){
        xd <- vector("list", length(groups))
        for (grp in groups){
            if(missing(width)) width <- width.SJ(data[g == grp])
            xd1 <- density(data[g == grp], n = 200L, width = width,
                           from = xlim[1L], to = xlim[2L])
            ymax <- max(ymax, xd1$y)
            xd[[grp]] <- xd1
        }
    }
    dev.hold(); on.exit(dev.flush())
    if(!sep)
        plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = "", bty = bty)
    else {
        oldpar <- par(mfrow = c(length(groups), 1L))
        on.exit(par(oldpar), add = TRUE)
    }
    for (grp in groups) {
        if(sep) plot(xlim, c(0, ymax), type = "n",
                     xlab = paste("group", grp), ylab = "", bty = bty)
        if(type == "histogram" || type == "both") {
            n <- length(breaks)
            rect(breaks[-n], 0, breaks[-1L], est[[grp]], col = col, ...)
        }
        if(type == "density" || type  ==  "both") lines(xd[[grp]])
    }
    invisible()
}

pairs.gipslda <- function(x, labels = colnames(x), panel = panel.gipslda,
                      dimen, abbrev = FALSE, ..., cex = 0.7,
                      type = c("std", "trellis"))
{
    panel.gipslda <- function(x,y, ...) text(x, y, as.character(g), cex = cex, ...)
    type <- match.arg(type)
    if(!is.null(Terms <- x$terms)) { # formula fit
        data <- model.frame(x)
        X <- model.matrix(delete.response(Terms), data)
        g <- model.response(data)
        xint <- match("(Intercept)", colnames(X), nomatch = 0L)
        if(xint > 0L) X <- X[, -xint, drop = FALSE]
    } else { # matrix or data-frame fit
        xname <- x$call$x
        gname <- x$call[[3L]]
        if(!is.null(sub <- x$call$subset)) {
            X <- eval.parent(parse(text=paste(deparse(xname, backtick=TRUE),
                                   "[", deparse(sub, backtick=TRUE),",]")))
            g <- eval.parent(parse(text=paste(deparse(gname, backtick=TRUE),
                                   "[", deparse(sub, backtick=TRUE),"]")))
        } else {
            X <- eval.parent(xname)
            g <- eval.parent(gname)
        }
        if(!is.null(nas <- x$call$na.action)) {
            df <- data.frame(g = g, X = X)
            df <- eval(call(nas, df))
            g <- df$g
            X <- df$X
        }
    }
    g <- as.factor(g)
    if(abbrev) levels(g) <- abbreviate(levels(g), abbrev)
    means <- colMeans(x$means)
    X <- scale(X, center = means, scale = FALSE) %*% x$scaling
    if(!missing(dimen) && dimen < ncol(X)) X <- X[, 1L:dimen]
    if(type == "std") pairs(X, panel = panel, ...)
    else {
        print(lattice::splom(~X, groups = g, panel = lattice::panel.superpose,
                             key = list(
                             text = list(levels(g)),
                             points = lattice::Rows(lattice::trellis.par.get("superpose.symbol"),
                             seq_along(levels(g))),
                             columns = min(5L, length(levels(g)))
                             )
                    ))
    }
    invisible(NULL)
}

model.frame.gipslda <- function(formula, ...)
{
    oc <- formula$call
    oc$prior <- oc$tol <- oc$method <- oc$CV <- oc$nu <- NULL
    oc[[1L]] <- quote(stats::model.frame)
    if(length(dots <- list(...))) {
        nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
        oc[names(nargs)] <- nargs
    }
    if (is.null(env <- environment(formula$terms))) env <- parent.frame()
    eval(oc, env)
}

coef.gipslda <- function(object, ...) object$scaling
