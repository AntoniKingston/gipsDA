qda_fitters <- list(
  gipsqda = list(
    generic = gipsDA::gipsqda,
    default = gipsDA:::gipsqda.default,
    class = "gipsqda"
  ),
  gipsmultqda = list(
    generic = gipsDA::gipsmultqda,
    default = gipsDA:::gipsmultqda.default,
    class = "gipsmultqda"
  )
)

fit_with <- function(spec, x, grouping, ...) {
  suppressWarnings(spec$generic(x, grouping, optimizer = "BF", ...))
}

expect_same_qda_estimates <- function(x, y) {
  expect_equal(x$prior, y$prior)
  expect_equal(x$counts, y$counts)
  expect_equal(x$means, y$means)
  expect_equal(x$scaling, y$scaling, tolerance = 1e-10)
  expect_equal(x$ldet, y$ldet, tolerance = 1e-10)
}

test_that("QDA generics dispatch over matrix, data frame, formula, and default", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)

  for (name in names(qda_fitters)) {
    spec <- qda_fitters[[name]]

    matrix_fit <- fit_with(spec, fixture$x, fixture$grouping)
    data_frame_fit <- fit_with(spec, as.data.frame(fixture$x), fixture$grouping)
    formula_fit <- suppressWarnings(spec$generic(
      class ~ x1 + x2,
      data = fixture$data,
      optimizer = "BF"
    ))
    default_fit <- suppressWarnings(spec$default(
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    ))

    expect_valid_qda_fit(matrix_fit, spec$class, 18, 2, 3)
    expect_valid_qda_fit(data_frame_fit, spec$class, 18, 2, 3)
    expect_valid_qda_fit(formula_fit, spec$class, 18, 2, 3)
    expect_valid_qda_fit(default_fit, spec$class, 18, 2, 3)
    expect_equal(matrix_fit$means, data_frame_fit$means)
    expect_equal(matrix_fit$means, formula_fit$means)
    expect_equal(matrix_fit$means, default_fit$means)
    expect_s3_class(formula_fit$terms, "terms")
    expect_identical(as.character(formula_fit$call[[1L]]), name)
    expect_identical(as.character(data_frame_fit$call[[1L]]), name)
  }
})

test_that("subset and formula na.action are honored", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  keep <- c(1:6, 9:14, 17:22)
  data_na <- fixture$data
  data_na[c(2, 10, 18), "x1"] <- NA_real_

  for (spec in qda_fitters) {
    subset_fit <- fit_with(
      spec,
      fixture$x,
      fixture$grouping,
      subset = keep
    )
    formula_fit <- suppressWarnings(spec$generic(
      class ~ x1 + x2,
      data = data_na,
      subset = keep,
      na.action = stats::na.omit,
      optimizer = "BF"
    ))

    expect_valid_qda_fit(subset_fit, spec$class, length(keep), 2, 3)
    expect_equal(unname(subset_fit$counts), rep(6, 3))
    expect_valid_qda_fit(formula_fit, spec$class, 15, 2, 3)
    expect_equal(unname(formula_fit$counts), rep(5, 3))
    expect_s3_class(formula_fit$na.action, "omit")
  }
})

for (name in names(qda_fitters)) {
  test_that(paste(name, "matrix na.action omits only incomplete rows"), {
    fixture <- make_classification_fixture(p = 2, n_per_class = 8)
    matrix_na <- fixture$x
    matrix_na[c(2, 10, 18), 1] <- NA_real_
    spec <- qda_fitters[[name]]

    fit <- fit_with(
      spec,
      matrix_na,
      fixture$grouping,
      na.action = stats::na.omit
    )

    expect_valid_qda_fit(fit, spec$class, 21, 2, 3)
    expect_equal(unname(fit$counts), rep(7, 3))
  })
}

test_that("explicit and omitted priors have the intended values", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  keep <- c(1:8, 9:14, 17:20)
  expected_proportions <- c(setosa = 8 / 18, versicolor = 6 / 18, virginica = 4 / 18)
  supplied <- c(setosa = 0.2, versicolor = 0.3, virginica = 0.5)

  for (spec in qda_fitters) {
    omitted <- fit_with(
      spec,
      fixture$x,
      fixture$grouping,
      subset = keep
    )
    explicit <- fit_with(
      spec,
      fixture$x,
      fixture$grouping,
      subset = keep,
      prior = supplied
    )

    expect_equal(omitted$prior, expected_proportions)
    expect_equal(explicit$prior, supplied)
  }
})

test_that("MAP can be toggled and nu currently leaves estimates unchanged", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  for (spec in qda_fitters) {
    map_fit <- fit_with(spec, fixture$x, fixture$grouping, MAP = TRUE, nu = 2)
    averaged_fit <- fit_with(spec, fixture$x, fixture$grouping, MAP = FALSE, nu = 2)
    low_nu <- fit_with(spec, fixture$x, fixture$grouping, nu = 1)
    high_nu <- fit_with(spec, fixture$x, fixture$grouping, nu = 100)

    expect_valid_qda_fit(map_fit, spec$class, 12, 2, 2)
    expect_valid_qda_fit(averaged_fit, spec$class, 12, 2, 2)
    expect_same_qda_estimates(low_nu, high_nu)
  }
})

test_that("MH supplies a default max_iter only when it is omitted", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 4)
  seen_max_iter <- list()

  fake_project_covs <- function(emp_covs, ns_obs, MAP, optimizer, max_iter, ...) {
    seen_max_iter[[length(seen_max_iter) + 1L]] <<- max_iter
    list(covs = emp_covs, opt_info = "mock optimization")
  }
  testthat::local_mocked_bindings(
    project_covs = fake_project_covs,
    .package = "gipsDA"
  )

  for (spec in qda_fitters) {
    expect_warning(
      spec$generic(fixture$x, fixture$grouping, optimizer = "MH"),
      "Setting max_iter = 100"
    )
    expect_no_warning(
      spec$generic(
        fixture$x,
        fixture$grouping,
        optimizer = "MH",
        max_iter = 7
      )
    )
  }

  expect_equal(unlist(seen_max_iter), c(100, 100, 7, 7, 100, 7))
})

test_that("QDA fit validation rejects malformed inputs and priors", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 5)

  for (spec in qda_fitters) {
    expect_error(
      spec$generic(as.numeric(fixture$x), fixture$grouping, optimizer = "BF"),
      "'x' is not a matrix",
      fixed = TRUE
    )

    bad_x <- fixture$x
    bad_x[1, 1] <- NA_real_
    expect_error(
      spec$generic(bad_x, fixture$grouping, optimizer = "BF"),
      "infinite, NA or NaN values"
    )
    expect_error(
      spec$generic(fixture$x, fixture$grouping[-1], optimizer = "BF"),
      "nrow\\(x\\) and length\\(grouping\\) are different"
    )
    expect_error(
      spec$generic(
        fixture$x,
        fixture$grouping,
        prior = c(-0.1, 1.1),
        optimizer = "BF"
      ),
      "invalid 'prior'",
      fixed = TRUE
    )
    expect_error(
      spec$generic(
        fixture$x,
        fixture$grouping,
        prior = c(0.2, 0.2),
        optimizer = "BF"
      ),
      "invalid 'prior'",
      fixed = TRUE
    )
    expect_error(
      spec$generic(
        fixture$x,
        fixture$grouping,
        prior = c(0.2, 0.3, 0.5),
        optimizer = "BF"
      ),
      "'prior' is of incorrect length",
      fixed = TRUE
    )
  }
})

test_that("QDA fits expose stable dimensions, labels, and finite estimates", {
  fixture <- make_classification_fixture(p = 3, n_per_class = 6)

  for (spec in qda_fitters) {
    fit <- fit_with(spec, fixture$x, fixture$grouping)

    expect_valid_qda_fit(fit, spec$class, 18, 3, 3)
    expect_named(fit$counts, levels(fixture$grouping))
    expect_named(fit$prior, levels(fixture$grouping))
    expect_equal(rownames(fit$means), levels(fixture$grouping))
    expect_equal(colnames(fit$means), colnames(fixture$x))
    expect_equal(dimnames(fit$scaling)[[1L]], colnames(fixture$x))
    expect_equal(dimnames(fit$scaling)[[3L]], levels(fixture$grouping))
    expect_true(all(is.finite(fit$means)))
    expect_true(all(is.finite(fit$ldet)))
    expect_false(is.null(fit$optimization_info))
  }
})
