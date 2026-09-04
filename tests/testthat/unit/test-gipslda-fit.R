test_that("gipslda dispatches from matrix, data frame, formula, and default", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)

  matrix_fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")
  frame_fit <- gipslda(
    fixture$data[c("x1", "x2")],
    fixture$grouping,
    optimizer = "BF"
  )
  formula_fit <- gipslda(class ~ ., fixture$data, optimizer = "BF")
  default_fit <- gipsDA:::gipslda.default(
    unclass(fixture$x),
    fixture$grouping,
    optimizer = "BF"
  )

  for (fit in list(matrix_fit, frame_fit, formula_fit, default_fit)) {
    expect_valid_lda_fit(fit, n = 18, p = 2, groups = 3)
  }
  expect_equal(unname(frame_fit$means), unname(matrix_fit$means))
  expect_equal(unname(formula_fit$means), unname(matrix_fit$means))
  expect_equal(unname(default_fit$means), unname(matrix_fit$means))
  expect_s3_class(formula_fit$terms, "terms")
  expect_length(formula_fit$xlevels, 0)
})

test_that("matrix and formula interfaces honor subset and na.action", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  keep <- c(1:6, 9:14, 17:22)

  matrix_fit <- gipslda(
    fixture$x,
    fixture$grouping,
    subset = keep,
    optimizer = "BF"
  )
  expect_equal(matrix_fit$N, length(keep))
  expect_equal(unname(matrix_fit$counts), c(6, 6, 6))

  data_na <- fixture$data
  data_na$x1[c(2, 10, 18)] <- NA_real_
  formula_fit <- gipslda(
    class ~ .,
    data_na,
    subset = keep,
    na.action = stats::na.omit,
    optimizer = "BF"
  )
  expect_equal(formula_fit$N, length(keep) - 3)
  expect_s3_class(formula_fit$na.action, "omit")
  expect_equal(as.integer(formula_fit$na.action), c(2, 8, 14))

  x_na <- fixture$x
  x_na[c(2, 10, 18), 1] <- NA_real_
  matrix_na_fit <- gipslda(
    x_na,
    fixture$grouping,
    na.action = stats::na.omit,
    optimizer = "BF"
  )
  expect_equal(matrix_na_fit$N, nrow(x_na) - 3)
  expect_equal(unname(matrix_na_fit$counts), c(7, 7, 7))
})

test_that("explicit and omitted priors have their documented behavior", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)
  explicit <- c(setosa = 0.2, versicolor = 0.3, virginica = 0.5)
  unequal_rows <- c(1:4, 7:11, 13:18)
  unequal_x <- fixture$x[unequal_rows, , drop = FALSE]
  unequal_grouping <- droplevels(fixture$grouping[unequal_rows])

  explicit_fit <- gipslda(
    fixture$x,
    fixture$grouping,
    prior = explicit,
    optimizer = "BF"
  )
  omitted_fit <- gipslda(unequal_x, unequal_grouping, optimizer = "BF")

  expect_equal(explicit_fit$prior, explicit)
  # This is the intended public contract for an omitted prior.
  expect_equal(
    omitted_fit$prior,
    prop.table(table(unequal_grouping)),
    ignore_attr = TRUE
  )
  expect_equal(names(omitted_fit$prior), levels(unequal_grouping))
})

test_that("weighted covariance and MAP modes produce usable fits", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)

  pooled <- gipslda(
    fixture$x,
    fixture$grouping,
    weighted_avg = FALSE,
    MAP = TRUE,
    optimizer = "BF"
  )
  weighted <- gipslda(
    fixture$x,
    fixture$grouping,
    weighted_avg = TRUE,
    MAP = TRUE,
    optimizer = "BF"
  )
  averaged <- gipslda(
    fixture$x,
    fixture$grouping,
    MAP = FALSE,
    optimizer = "BF"
  )

  for (fit in list(pooled, weighted, averaged)) {
    expect_valid_lda_fit(fit, n = 18, p = 2, groups = 3)
  }
  expect_false(is.null(pooled$optimization_info))
  expect_false(is.null(averaged$optimization_info))
})

test_that("MH without max_iter warns and still returns a fit", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  expect_warning(
    fit <- gipslda(
      fixture$x,
      fixture$grouping,
      optimizer = "MH"
    ),
    "max_iter.*unspecified"
  )
  expect_valid_lda_fit(fit, n = 12, p = 2, groups = 2)
})

test_that("fit output has stable structural invariants", {
  fixture <- make_classification_fixture(p = 3, n_per_class = 6)
  fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")

  expect_named(
    fit,
    c(
      "prior", "counts", "means", "scaling", "lev", "svd", "N",
      "optimization_info", "selected_map_permutation", "fit_info", "call"
    )
  )

  expect_named(fit$fit_info, c("MAP", "optimizer", "max_iter", "weighted_avg"))
  expect_true(is.logical(fit$fit_info$MAP))
  expect_identical(fit$fit_info$optimizer, "BF")
  expect_null(fit$fit_info$max_iter)
  expect_false(fit$fit_info$weighted_avg)

  expect_valid_lda_fit(fit, n = 18, p = 3, groups = 3)
  expect_equal(names(fit$counts), levels(fixture$grouping))
  expect_equal(rownames(fit$means), levels(fixture$grouping))
  expect_equal(rownames(fit$scaling), colnames(fixture$x))
  expect_equal(colnames(fit$scaling), paste0("LD", seq_along(fit$svd)))
  expect_lte(ncol(fit$scaling), length(levels(fixture$grouping)) - 1)
  expect_true(all(fit$svd > 0))
})

test_that("fit validates inputs and degenerate data", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  expect_error(
    gipsDA:::gipslda.default(
      seq_len(nrow(fixture$x)),
      fixture$grouping,
      optimizer = "BF"
    ),
    "'x' is not a matrix",
    fixed = TRUE
  )
  bad <- fixture$x
  bad[1, 1] <- NA_real_
  expect_error(
    gipslda(bad, fixture$grouping, optimizer = "BF"),
    "infinite, NA or NaN"
  )
  bad[1, 1] <- Inf
  expect_error(
    gipslda(bad, fixture$grouping, optimizer = "BF"),
    "infinite, NA or NaN"
  )
  expect_error(
    gipslda(fixture$x, fixture$grouping[-1], optimizer = "BF"),
    "nrow\\(x\\).*length\\(grouping\\)"
  )
  expect_error(
    gipslda(
      fixture$x,
      fixture$grouping,
      prior = c(0.5, -0.5),
      optimizer = "BF"
    ),
    "invalid 'prior'",
    fixed = TRUE
  )
  expect_error(
    gipslda(
      fixture$x,
      fixture$grouping,
      prior = c(0.4, 0.4),
      optimizer = "BF"
    ),
    "invalid 'prior'",
    fixed = TRUE
  )
  expect_error(
    gipslda(
      fixture$x,
      fixture$grouping,
      prior = c(0.2, 0.3, 0.5),
      optimizer = "BF"
    ),
    "incorrect length"
  )

  constant <- cbind(x1 = rep(c(0, 1), each = 6), x2 = rep(1, 12))
  expect_error(
    gipslda(constant, fixture$grouping, optimizer = "BF"),
    "constant within groups"
  )

  identical_means <- rbind(
    cbind(x1 = 1:6, x2 = c(1, 4, 2, 6, 3, 5)),
    cbind(x1 = 1:6, x2 = c(1, 4, 2, 6, 3, 5))
  )
  expect_error(
    gipslda(identical_means, fixture$grouping, optimizer = "BF"),
    "group means are numerically identical"
  )
})

test_that("empty grouping levels are dropped with a warning", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)
  grouping <- factor(fixture$grouping, levels = c(levels(fixture$grouping), "empty"))

  warnings <- character()
  fit <- withCallingHandlers(
    gipslda(
      fixture$x,
      grouping,
      prior = c(0.45, 0.55, 0),
      optimizer = "BF"
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("group empty is empty", warnings, fixed = TRUE)))
  expect_equal(fit$prior, c(setosa = 0.45, versicolor = 0.55))
  expect_equal(fit$lev, c("setosa", "versicolor", "empty"))
  expect_equal(names(fit$counts), c("setosa", "versicolor"))
})
