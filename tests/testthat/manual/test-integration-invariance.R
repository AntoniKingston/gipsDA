invariance_fitters <- list(
  gipslda = gipsDA::gipslda,
  gipsqda = gipsDA::gipsqda,
  gipsmultqda = gipsDA::gipsmultqda
)

fit_invariance_model <- function(fitter, x, grouping, ...) {
  suppressWarnings(fitter(x, grouping, optimizer = "BF", ...))
}

posterior_for <- function(object, newdata) {
  suppressWarnings(predict(object, newdata)$posterior)
}

test_that("fits are invariant to row order and equivariant to class labels", {
  fixture <- make_binary_fixture(p = 3, n_per_class = 6)
  x <- fixture$x
  grouping <- fixture$grouping
  evaluation_x <- x[c(2, 5, 8, 11), , drop = FALSE]
  row_order <- c(seq(2, nrow(x), by = 2), seq(1, nrow(x), by = 2))

  renamed_values <- c(
    setosa = "zeta",
    versicolor = "alpha"
  )[as.character(grouping)]
  renamed_grouping <- factor(renamed_values, levels = c("alpha", "zeta"))
  renamed_to_original <- c(alpha = "versicolor", zeta = "setosa")

  for (fitter in invariance_fitters) {
    original <- fit_invariance_model(fitter, x, grouping)
    reordered <- fit_invariance_model(
      fitter,
      x[row_order, , drop = FALSE],
      grouping[row_order]
    )
    renamed <- fit_invariance_model(fitter, x, renamed_grouping)

    original_posterior <- posterior_for(original, evaluation_x)
    reordered_posterior <- posterior_for(reordered, evaluation_x)
    renamed_posterior <- posterior_for(renamed, evaluation_x)

    expect_equal(reordered$means, original$means, tolerance = 1e-12)
    expect_equal(reordered$prior, original$prior, tolerance = 1e-12)
    expect_equal(reordered_posterior, original_posterior, tolerance = 1e-9)

    semantic_order <- names(original$prior)
    renamed_columns <- names(renamed_to_original)[
      match(semantic_order, renamed_to_original)
    ]
    expect_equal(
      unname(renamed$means[renamed_columns, , drop = FALSE]),
      unname(original$means),
      tolerance = 1e-12
    )
    expect_equal(
      unname(renamed$prior[renamed_columns]),
      unname(original$prior),
      tolerance = 1e-12
    )
    expect_equal(
      unname(renamed_posterior[, renamed_columns, drop = FALSE]),
      unname(original_posterior),
      tolerance = 1e-9
    )
    expect_identical(
      unname(renamed_to_original[as.character(predict(renamed, evaluation_x)$class)]),
      as.character(predict(original, evaluation_x)$class)
    )
  }
})

test_that("fits are equivariant to a consistent feature permutation", {
  fixture <- make_binary_fixture(p = 3, n_per_class = 6)
  feature_order <- c(3, 1, 2)
  evaluation_rows <- c(1, 4, 7, 10, 12)
  evaluation_x <- fixture$x[evaluation_rows, , drop = FALSE]

  for (fitter in invariance_fitters) {
    original <- fit_invariance_model(fitter, fixture$x, fixture$grouping)
    permuted <- fit_invariance_model(
      fitter,
      fixture$x[, feature_order, drop = FALSE],
      fixture$grouping
    )

    expect_equal(
      unname(permuted$means),
      unname(original$means[, feature_order, drop = FALSE]),
      tolerance = 1e-12
    )
    expect_equal(
      unname(posterior_for(
        permuted,
        evaluation_x[, feature_order, drop = FALSE]
      )),
      unname(posterior_for(original, evaluation_x)),
      tolerance = 1e-8
    )
    expect_identical(
      as.character(predict(
        permuted,
        evaluation_x[, feature_order, drop = FALSE]
      )$class),
      as.character(predict(original, evaluation_x)$class)
    )
  }
})

test_that("uniform observation duplication preserves valid model quantities", {
  fixture <- make_binary_fixture(p = 3, n_per_class = 6)
  duplicated_rows <- rep(seq_len(nrow(fixture$x)), each = 2L)

  for (fitter in invariance_fitters) {
    original <- fit_invariance_model(fitter, fixture$x, fixture$grouping)
    duplicated <- fit_invariance_model(
      fitter,
      fixture$x[duplicated_rows, , drop = FALSE],
      fixture$grouping[duplicated_rows]
    )

    expect_equal(duplicated$means, original$means, tolerance = 1e-12)
    expect_equal(duplicated$prior, original$prior, tolerance = 1e-12)
    expect_equal(duplicated$counts, 2 * original$counts)
    expect_equal(duplicated$N, 2 * original$N)
  }
})

test_that("MAP projection is invariant under its selected permutation", {
  empirical_covariance <- matrix(
    c(
      2.0, 0.4, 0.1,
      0.4, 1.5, 0.3,
      0.1, 0.3, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  )

  result <- gipsDA:::project_covs(
    list(empirical_covariance),
    ns_obs = 12,
    MAP = TRUE,
    optimizer = "BF",
    max_iter = NULL
  )
  projected <- result$covs[[1L]]

  expect_equal(
    gips::project_matrix(projected, result$permutation),
    projected,
    tolerance = 1e-12
  )
})

test_that("MAP FALSE projection is the posterior-weighted projection", {
  empirical_covariance <- matrix(
    c(
      2.0, 0.4, 0.1,
      0.4, 1.5, 0.3,
      0.1, 0.3, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  )

  result <- gipsDA:::project_covs(
    list(empirical_covariance),
    ns_obs = 12,
    MAP = FALSE,
    optimizer = "BF",
    max_iter = NULL
  )
  probabilities <- result$opt_info
  independently_reconstructed <- Reduce(
    `+`,
    Map(
      function(probability, permutation) {
        probability * gips::project_matrix(empirical_covariance, permutation)
      },
      as.numeric(probabilities),
      names(probabilities)
    )
  ) / sum(probabilities)

  expect_equal(
    result$covs[[1L]],
    independently_reconstructed,
    tolerance = 1e-12
  )
})
