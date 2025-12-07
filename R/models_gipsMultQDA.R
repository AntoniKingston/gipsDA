
gipsmultqda <- function(x, ...) UseMethod("gipsmultqda")

gipsmultqda.formula <- function(formula, data, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1L]] <- quote(stats::model.frame)
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1L]
    if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
    xint <- match("(Intercept)", colnames(x), nomatch=0L)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
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

gipsmultqda.data.frame <- function(x, ...)
{
    res <- gipsmultqda(structure(data.matrix(x), class="matrix"), ...)
    cl <- match.call()
    cl[[1L]] <- as.name("gipsmultqda")
    res$call <- cl
    res
}


gipsmultqda.matrix <- function(x, grouping, ..., subset, na.action)
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
#    res <- NextMethod("gipsmultqda")
    res <- gipsmultqda.default(x, grouping, ...)
    cl <- match.call()
    cl[[1L]] <- as.name("gipsmultqda")
    res$call <- cl
    res
}

gipsmultqda.default <-
  function(x, grouping, prior = proportions, nu = 5, MAP = TRUE, optimizer = NULL, max_iter = NULL, ...)
{
    if(is.null(dim(x))) stop("'x' is not a matrix")
    x <- as.matrix(x)
    if(any(!is.finite(x)))
        stop("infinite, NA or NaN values in 'x'")
    n <- nrow(x)
    p <- ncol(x)
    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")
    g <- as.factor(grouping)
    lev <- levels(g)
    counts <- as.vector(table(g))
    names(counts) <- lev
    if(any(counts < p+1)) stop("some group is too small for 'gipsmultqda'")
    proportions <- counts/length(g)
    ng <- length(proportions)
# allow for supplied prior
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if(length(prior) != ng) stop("'prior' is of incorrect length")
    names(prior) <- lev

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

# means by group (rows) and variable (columns)
    group.means <- tapply(x, list(rep(g, ncol(x)), col(x)), mean)
    scaling <- array(dim=c(p,p,ng))
    ldet <- numeric(ng)
####################################################################################
    cXs <- list()

    for (i in 1L:ng){
        cX <- MASS::cov.mve(x[unclass(g) == i, ])
        cXs[[i]] <- cX$cov
        group.means[i,] <- cX$center
    }
    pr_cov_opt_info <- project_covs(cXs, counts, MAP, optimizer, max_iter)
    cov_proj <- pr_cov_opt_info$covs
    optimization_info <- pr_cov_opt_info$opt_info
    for (i in 1L:ng) {
      cov_proj <- cXs[[i]]
      sX <- svd(cov_proj, nu=0)
      scaling[, , i] <- sX$v %*% diag(sqrt(1/sX$d),,p)
      ldet[i] <- sum(log(sX$d))
    }
####################################################################################
    if(is.null(dimnames(x)))
        dimnames(scaling) <- list(NULL, as.character(1L:p), lev)
    else {
        dimnames(scaling) <- list(colnames(x), as.character(1L:p), lev)
        dimnames(group.means)[[2L]] <- colnames(x)
    }
    cl <- match.call()
    cl[[1L]] <- as.name("gipsmultqda")
    res <- list(prior = prior, counts = counts, means = group.means,
                scaling = scaling, ldet = ldet, lev = lev, N = n, call = cl, optimization_info = optimization_info)
    class(res) <- "gipsmultqda"
    res
}

predict.gipsmultqda <- function(object, newdata, prior = object$prior,
			method = c("plug-in", "predictive", "debiased",
                          "looCV"), ...)
{
    if(!inherits(object, "gipsmultqda")) stop("object not of class \"gipsmultqda\"")
    method <- match.arg(method)
    if(method == "looCV" && !missing(newdata))
        stop("cannot have leave-one-out CV with 'newdata'")
    if(is.null(mt <- object$call$method)) mt <- "moment"
    if(method == "looCV" && !(mt == "moment" || mt == "mle"))
        stop(gettext("cannot use leave-one-out CV with method %s",
                     sQuote(mt), domain = NA))
    ngroup <- length(object$prior)
    if(!missing(prior)) {
        if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
        if(length(prior) != ngroup) stop("'prior' is of incorrect length")
    }
    if(!is.null(Terms <- object$terms)) {
    # formula fit
        if(missing(newdata)) newdata <- model.frame(object)
        else {
            newdata <- model.frame(as.formula(delete.response(Terms)),
                                   newdata, na.action=function(x) x,
                                   xlev = object$xlevels)
        }
        x <- model.matrix(delete.response(Terms), newdata,
                          contrasts.arg = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch=0L)
        if(xint > 0) x <- x[, -xint, drop=FALSE]
        if(method == "looCV") g <- model.response(newdata)
    } else { #
    # matrix or data-frame fit
        if(missing(newdata)) {
            if(!is.null(sub <- object$call$subset)) {
                newdata <-
                    eval.parent(parse(text=paste(deparse(object$call$x,
                                      backtick=TRUE),
                                      "[", deparse(sub, backtick=TRUE),",]")))
                g <- eval.parent(parse(text=paste(deparse(object$call[[3L]],
                                       backtick=TRUE),
                                       "[", deparse(sub, backtick=TRUE),"]")))
            } else {
                newdata <- eval.parent(object$call$x)
                g <- eval.parent(object$call[[3L]])
            }
            if(!is.null(nas <- object$call$na.action)) {
                df <- data.frame(g = g, X = newdata)
                df <- eval(call(nas, df))
                g <- df$g
                newdata <- df$X
            }
            g <- as.factor(g)
        }
        if(is.null(dim(newdata)))
            dim(newdata) <- c(1, length(newdata))  # a row vector
        x <- as.matrix(newdata)		# to cope with dataframes
    }
    p <- ncol(object$means)
    if(ncol(x) != p) stop("wrong number of variables")
    if(length(colnames(x)) > 0L &&
       any(colnames(x) != dimnames(object$means)[[2L]]))
        warning("variable names in 'newdata' do not match those in 'object'")
    dist <- matrix(0, nrow = nrow(x), ncol = ngroup)
    if(method == "plug-in") {
        for(i in 1L:ngroup) {
            dev <- ((x - matrix(object$means[i,  ], nrow(x),
                                ncol(x), byrow = TRUE)) %*% object$scaling[,,i])
            dist[, i] <- 0.5 * rowSums(dev^2) + 0.5 * object$ldet[i] - log(prior[i])
        }
#        dist <- exp( -(dist - min(dist, na.rm=T)))
        dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
    } else if(method == "looCV") {
        n <- nrow(x)
        NG <- 1
        if(mt == "mle") NG <- 0
        ldet <- matrix(0, n, ngroup)
        for(i in 1L:ngroup) {
            dev <- ((x - matrix(object$means[i,  ], nrow(x), p, byrow = TRUE))
                    %*% object$scaling[,,i])
            dist[, i] <- rowSums(dev^2)
            ldet[, i] <- object$ldet[i]
        }
        nc <- object$counts[g]
        ind <- cbind(1L:n, g)
        fac <- 1 - nc/(nc-1)/(nc-NG) * dist[ind]
        fac[] <- pmax(fac, 1e-10)  # possibly degenerate dsn
        ldet[ind] <- log(fac) + p * log((nc-NG)/(nc-1-NG)) + ldet[ind]
        dist[ind] <- dist[ind] * (nc^2/(nc-1)^2) * (nc-1-NG)/(nc-NG) / fac
        dist <- 0.5 * dist + 0.5 * ldet -
            matrix(log(prior), n, ngroup, byrow=TRUE)
        dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
    } else if(method == "debiased") {
        for(i in 1L:ngroup) {
            nk <- object$counts[i]
            Bm <- p * log((nk-1)/2) - sum(digamma(0.5 * (nk - 1L:ngroup)))
            dev <- ((x - matrix(object$means[i,  ], nrow = nrow(x),
                                ncol = ncol(x), byrow = TRUE)) %*% object$scaling[,,i])
            dist[, i] <- 0.5 * (1 - (p-1)/(nk-1)) * rowSums(dev^2) +
                0.5 * object$ldet[i] - log(prior[i]) + 0.5 * Bm - p/(2*nk)
        }
        dist <- exp( -(dist - apply(dist, 1L, min, na.rm=TRUE)))
    } else {
        for(i in 1L:ngroup) {
            nk <- object$counts[i]
            dev <- ((x - matrix(object$means[i,  ], nrow = nrow(x),
                                ncol = ncol(x), byrow = TRUE))
                    %*% object$scaling[,,i])
            dev <- 1 + rowSums(dev^2)/(nk+1)
            dist[, i] <- prior[i] * exp(-object$ldet[i]/2) *
                dev^(-nk/2) * (1 + nk)^(-p/2)
        }
    }
    posterior <- dist/drop(dist %*% rep(1, ngroup))
    cl <- factor(max.col(posterior), levels=seq_along(object$lev),
                 labels=object$lev)
    dimnames(posterior) <- list(rownames(x), object$lev)
    list(class = cl, posterior = posterior)
}

print.gipsmultqda <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        names(cl)[2L] <- ""
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    cat("\nPrior probabilities of groups:\n")
    print(x$prior, ...)
    cat("\nGroup means:\n")
    print(x$means, ...)
    cat("\nPermutations with their estimated probabilities:\n")
    print(x$optimization_info)
    invisible(x)
}

model.frame.gipsmultqda <-  model.frame.gipslda
