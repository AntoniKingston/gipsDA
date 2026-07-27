# This test independently reconstructs the fitted parameters of all three
# classifiers. The reference calculations below deliberately use only base R,
# stats, and the public gips API; no gipsDA helper or internal function is used.

manual_r3_data <- function() {
  # Generate three heteroscedastic Gaussian populations in R^3. Unequal class
  # sizes exercise the empirical-prior calculation and the sample-size vector
  # used by the joint covariance projection.
  set.seed(20260727)
  sizes <- c(A = 9L, B = 10L, C = 11L)
  means <- rbind(
    A = c(-1.5, 0.2, 0.8),
    B = c(0.6, 1.7, -0.5),
    C = c(1.4, -1.0, 1.2)
  )
  covariances <- list(
    A = matrix(c(
      1.20, 0.35, 0.10,
      0.35, 0.80, -0.15,
      0.10, -0.15, 0.65
    ), 3L, 3L, byrow = TRUE),
    B = matrix(c(
      0.75, -0.20, 0.18,
      -0.20, 1.35, 0.30,
      0.18, 0.30, 0.90
    ), 3L, 3L, byrow = TRUE),
    C = matrix(c(
      1.05, 0.25, -0.22,
      0.25, 0.70, 0.12,
      -0.22, 0.12, 1.25
    ), 3L, 3L, byrow = TRUE)
  )

  samples <- lapply(names(sizes), function(group) {
    z <- matrix(stats::rnorm(sizes[[group]] * 3L), ncol = 3L)
    sweep(z %*% chol(covariances[[group]]), 2L, means[group, ], `+`)
  })
  x <- do.call(rbind, samples)
  colnames(x) <- c("x1", "x2", "x3")

  list(
    x = x,
    grouping = factor(
      rep(names(sizes), times = sizes),
      levels = names(sizes)
    )
  )
}

manual_desingularize <- function(A, target = 0.05) {
  # If lambda is the eigenvalue nearest zero, the package regularization is
  # A* = (A + s I) / (1 + s), where s = (target - lambda)/(1 - target).
  # The operation is skipped when |lambda| is already at least `target`.
  symmetric <- isTRUE(all.equal(A, t(A)))
  eigenvalues <- eigen(A, symmetric = symmetric, only.values = TRUE)$values
  lambda <- eigenvalues[which.min(abs(eigenvalues))]
  if (abs(lambda) >= target) {
    return(A)
  }
  s <- (target - lambda) / (1 - target)
  (A + diag(s, nrow(A))) / (1 + s)
}

manual_map_projection <- function(empirical_covariances, sample_sizes) {
  # For one covariance gips expects a matrix; for the joint model it expects
  # the complete list. `was_mean_estimated = TRUE` accounts for estimating
  # class means. BF enumerates all R^3 permutations and is deterministic.
  gips_input <- if (length(empirical_covariances) == 1L) {
    empirical_covariances[[1L]]
  } else {
    empirical_covariances
  }
  search <- gips::gips(
    gips_input,
    sample_sizes,
    was_mean_estimated = TRUE
  )
  search <- gips::find_MAP(
    search,
    optimizer = "BF",
    show_progress_bar = FALSE
  )
  permutation <- search[[1L]]

  list(
    covariances = lapply(
      empirical_covariances,
      gips::project_matrix,
      permutation
    ),
    permutation = permutation
  )
}

manual_group_statistics <- function(x, grouping) {
  groups <- levels(grouping)
  counts <- table(grouping)
  counts <- as.vector(counts)
  names(counts) <- groups
  means <- t(vapply(
    groups,
    function(group) colMeans(x[grouping == group, , drop = FALSE]),
    numeric(ncol(x))
  ))
  dimnames(means) <- list(groups, colnames(x))

  list(
    groups = groups,
    counts = counts,
    prior = counts / nrow(x),
    means = means
  )
}

manual_lda_parameters <- function(x, grouping, tol = 1e-4) {
  stats <- manual_group_statistics(x, grouping)
  n <- nrow(x)
  p <- ncol(x)
  number_of_groups <- length(stats$groups)

  # Let M_g be the class means and E_i = x_i - M_{g_i}. First standardize
  # columns by D = diag(1 / sqrt(var(E))). The unbiased pooled within-class
  # correlation estimate used for projection is
  #
  #   S_W = n/(n-G) cov(E D).
  centered <- x - stats$means[grouping, , drop = FALSE]
  within_sd <- sqrt(diag(stats::var(centered)))
  standardizer <- diag(1 / within_sd, nrow = p, ncol = p)
  pooled <- n / (n - number_of_groups) *
    stats::cov(centered %*% standardizer)

  projection <- manual_map_projection(list(pooled), n)
  projected <- manual_desingularize(projection$covariances[[1L]])

  # If S_W* = V diag(d) V', then D V diag(d^(-1/2)) whitens the
  # within-class covariance. A second SVD of the weighted, centered class
  # means supplies the Fisher discriminant directions.
  within_svd <- svd(projected, nu = 0L)
  within_rank <- sum(within_svd$d > tol^2)
  whitening <- standardizer %*%
    within_svd$v[, seq_len(within_rank), drop = FALSE] %*%
    diag(sqrt(1 / within_svd$d[seq_len(within_rank)]), within_rank)

  grand_mean <- colSums(stats$prior * stats$means)
  between <- sqrt(n * stats$prior / (number_of_groups - 1L)) *
    scale(stats$means, center = grand_mean, scale = FALSE) %*%
    whitening
  between_svd <- svd(between, nu = 0L)
  between_rank <- sum(between_svd$d > tol * between_svd$d[1L])
  scaling <- whitening %*%
    between_svd$v[, seq_len(between_rank), drop = FALSE]
  dimnames(scaling) <- list(
    colnames(x),
    paste0("LD", seq_len(between_rank))
  )

  c(stats, list(
    scaling = scaling,
    svd = between_svd$d[seq_len(between_rank)],
    permutation = projection$permutation,
    N = n
  ))
}

manual_qda_parameters <- function(x, grouping, joint = FALSE) {
  stats <- manual_group_statistics(x, grouping)
  n <- nrow(x)
  p <- ncol(x)
  number_of_groups <- length(stats$groups)
  empirical <- lapply(stats$groups, function(group) {
    stats::cov(x[grouping == group, , drop = FALSE])
  })

  if (joint) {
    # gipsmultqda performs one search from (S_1, ..., S_G), weighted by the
    # class sample sizes, and applies the resulting permutation to every S_g.
    projection <- manual_map_projection(empirical, stats$counts)
    projected <- projection$covariances
    permutations <- list(projection$permutation)
  } else {
    # gipsqda searches independently for each S_g. Its implementation supplies
    # total n to every single-covariance gips fit and stores the final class's
    # optimizer result in `optimization_info`.
    projections <- lapply(empirical, function(S) {
      manual_map_projection(list(S), n)
    })
    projected <- lapply(projections, function(result) result$covariances[[1L]])
    permutations <- lapply(projections, function(result) result$permutation)
  }
  projected <- lapply(projected, manual_desingularize)

  # For each projected covariance S_g* = V_g diag(d_g) V_g', QDA stores
  # L_g = V_g diag(d_g^(-1/2)) and log|S_g*| = sum(log(d_g)).
  # L_g L_g' is the sign-invariant precision matrix used for prediction.
  scaling <- array(dim = c(p, p, number_of_groups))
  ldet <- numeric(number_of_groups)
  for (i in seq_len(number_of_groups)) {
    covariance_svd <- svd(projected[[i]], nu = 0L)
    scaling[, , i] <- covariance_svd$v %*%
      diag(sqrt(1 / covariance_svd$d), p)
    ldet[i] <- sum(log(covariance_svd$d))
  }
  dimnames(scaling) <- list(
    colnames(x),
    as.character(seq_len(p)),
    stats$groups
  )

  c(stats, list(
    scaling = scaling,
    ldet = ldet,
    permutation = permutations[[length(permutations)]],
    N = n
  ))
}

expect_common_manual_parameters <- function(fit, reference, tolerance = 1e-10) {
  expect_equal(fit$prior, reference$prior, tolerance = tolerance)
  expect_equal(fit$counts, reference$counts, tolerance = tolerance)
  expect_equal(fit$means, reference$means, tolerance = tolerance)
  expect_identical(fit$lev, reference$groups)
  expect_identical(fit$N, reference$N)
  expect_identical(
    as.character(fit$optimization_info),
    as.character(reference$permutation)
  )
}

test_that("R3 model parameters equal independent base R and gips calculations", {
  data <- manual_r3_data()

  lda_reference <- manual_lda_parameters(data$x, data$grouping)
  qda_reference <- manual_qda_parameters(data$x, data$grouping, joint = FALSE)
  mult_reference <- manual_qda_parameters(data$x, data$grouping, joint = TRUE)

  lda_fit <- gipsDA::gipslda(
    data$x,
    data$grouping,
    MAP = TRUE,
    optimizer = "BF"
  )
  qda_fit <- gipsDA::gipsqda(
    data$x,
    data$grouping,
    MAP = TRUE,
    optimizer = "BF"
  )
  mult_fit <- gipsDA::gipsmultqda(
    data$x,
    data$grouping,
    MAP = TRUE,
    optimizer = "BF"
  )

  expect_common_manual_parameters(lda_fit, lda_reference)
  expect_common_manual_parameters(qda_fit, qda_reference)
  expect_common_manual_parameters(mult_fit, mult_reference)

  # SVD vectors are defined only up to column sign. Align signs before checking
  # LDA coefficients; for QDA compare L_g L_g', the precision matrix actually
  # used by the discriminant rule.
  signs <- sign(colSums(lda_fit$scaling * lda_reference$scaling))
  signs[signs == 0] <- 1
  expect_equal(
    sweep(lda_fit$scaling, 2L, signs, `*`),
    lda_reference$scaling,
    tolerance = 1e-10
  )
  expect_equal(lda_fit$svd, lda_reference$svd, tolerance = 1e-10)

  for (group in seq_along(qda_reference$groups)) {
    expect_equal(
      tcrossprod(qda_fit$scaling[, , group]),
      tcrossprod(qda_reference$scaling[, , group]),
      tolerance = 1e-10
    )
    expect_equal(
      tcrossprod(mult_fit$scaling[, , group]),
      tcrossprod(mult_reference$scaling[, , group]),
      tolerance = 1e-10
    )
  }
  expect_equal(qda_fit$ldet, qda_reference$ldet, tolerance = 1e-10)
  expect_equal(mult_fit$ldet, mult_reference$ldet, tolerance = 1e-10)
})
