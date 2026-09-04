test_that("print methods report key fitted model content", {
  binary_fixture <- make_binary_fixture(p = 2, n_per_class = 6)
  multi_fixture <- make_classification_fixture(p = 2, n_per_class = 6)

  lda_fit <- suppressWarnings(gipslda(
    multi_fixture$x,
    multi_fixture$grouping,
    optimizer = "BF"
  ))

  qda_fit <- suppressWarnings(gipsqda(
    binary_fixture$x,
    binary_fixture$grouping,
    optimizer = "BF"
  ))

  mult_fit <- suppressWarnings(gipsmultqda(
    binary_fixture$x,
    binary_fixture$grouping,
    optimizer = "BF"
  ))

  fits <- list(
    gipslda = lda_fit,
    gipsqda = qda_fit,
    gipsmultqda = mult_fit
  )

  for (model_name in names(fits)) {
    fit <- fits[[model_name]]

    output <- capture.output(return_value <- withVisible(print(fit)))

    expect_false(return_value$visible)
    expect_identical(return_value$value, fit)

    expect_true(any(grepl("Call:", output, fixed = TRUE)))
    expect_true(any(grepl(paste("Model:", model_name), output, fixed = TRUE)))
    expect_true(any(grepl("Number of observations:", output, fixed = TRUE)))
    expect_true(any(grepl("Number of groups:", output, fixed = TRUE)))
    expect_true(any(grepl("Number of predictors:", output, fixed = TRUE)))
    expect_true(any(grepl("Fitting options:", output, fixed = TRUE)))
    expect_true(any(grepl("Prior probabilities of groups:", output, fixed = TRUE)))
    expect_true(any(grepl("Class counts:", output, fixed = TRUE)))
    expect_true(any(grepl("Group means:", output, fixed = TRUE)))
    expect_true(any(grepl("Selected MAP permutation:", output, fixed = TRUE)))
    expect_true(any(grepl("Posterior probabilities of retained permutations:", output, fixed = TRUE)))

    if (inherits(fit, "gipslda")) {
      expect_true(any(grepl("Coefficients of linear discriminants:", output, fixed = TRUE)))
      expect_true(any(grepl("Proportion of trace:", output, fixed = TRUE)))
    }

    if (inherits(fit, "gipsqda")) {
      expect_true(any(grepl("Group: setosa", output, fixed = TRUE)))
      expect_true(any(grepl("Group: versicolor", output, fixed = TRUE)))
      expect_true(any(grepl("Log determinants of projected covariance matrices:", output, fixed = TRUE)))
    }

    if (inherits(fit, "gipsmultqda")) {
      expect_true(any(grepl("Log determinants of projected covariance matrices:", output, fixed = TRUE)))
    }
  }
})

test_that("summary methods return structured summary objects", {
  binary_fixture <- make_binary_fixture(p = 2, n_per_class = 6)
  multi_fixture <- make_classification_fixture(p = 2, n_per_class = 6)

  lda_fit <- suppressWarnings(gipslda(
    multi_fixture$x,
    multi_fixture$grouping,
    optimizer = "BF"
  ))

  qda_fit <- suppressWarnings(gipsqda(
    binary_fixture$x,
    binary_fixture$grouping,
    optimizer = "BF"
  ))

  mult_fit <- suppressWarnings(gipsmultqda(
    binary_fixture$x,
    binary_fixture$grouping,
    optimizer = "BF"
  ))

  lda_summary <- summary(lda_fit)
  qda_summary <- summary(qda_fit)
  mult_summary <- summary(mult_fit)

  expect_s3_class(lda_summary, "summary.gipslda")
  expect_s3_class(qda_summary, "summary.gipsqda")
  expect_s3_class(mult_summary, "summary.gipsmultqda")

  summaries <- list(
    lda_summary,
    qda_summary,
    mult_summary
  )

  for (summary_object in summaries) {
    expect_s3_class(summary_object, "summary.gipsda")

    expect_named(
      summary_object,
      c(
        "model", "call", "n", "p", "groups", "counts", "prior",
        "means", "fit_info", "optimization_info",
        setdiff(
          names(summary_object),
          c(
            "model", "call", "n", "p", "groups", "counts", "prior",
            "means", "fit_info", "optimization_info"
          )
        )
      ),
      ignore.order = TRUE
    )

    output <- capture.output(return_value <- withVisible(print(summary_object)))

    expect_false(return_value$visible)
    expect_identical(return_value$value, summary_object)

    expect_true(any(grepl("Call:", output, fixed = TRUE)))
    expect_true(any(grepl("Model:", output, fixed = TRUE)))
    expect_true(any(grepl("Number of observations:", output, fixed = TRUE)))
    expect_true(any(grepl("Number of groups:", output, fixed = TRUE)))
    expect_true(any(grepl("Number of predictors:", output, fixed = TRUE)))
    expect_true(any(grepl("Fitting options:", output, fixed = TRUE)))
    expect_true(any(grepl("Class counts:", output, fixed = TRUE)))
    expect_true(any(grepl("Prior probabilities of groups:", output, fixed = TRUE)))
    expect_true(any(grepl("Group means:", output, fixed = TRUE)))
  }

  expect_equal(lda_summary$model, "gipslda")
  expect_equal(qda_summary$model, "gipsqda")
  expect_equal(mult_summary$model, "gipsmultqda")

  expect_false(is.null(lda_summary$proportion_trace))
  expect_false(is.null(qda_summary$ldet))
  expect_false(is.null(mult_summary$ldet))
  expect_false(is.null(qda_summary$scaling_dim))
  expect_false(is.null(mult_summary$scaling_dim))
})

test_that("print helpers handle missing optional optimization and svd content", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  fit <- suppressWarnings(gipslda(
    fixture$x,
    fixture$grouping,
    optimizer = "BF"
  ))

  fit$optimization_info <- NULL
  fit$svd <- NULL

  output <- capture.output(print(fit))

  expect_true(any(grepl("Selected MAP permutation:", output, fixed = TRUE)))
  expect_true(any(grepl("Posterior probabilities of retained permutations:", output, fixed = TRUE)))
  expect_true(any(output == "NULL"))

  summary_object <- summary(fit)
  expect_null(summary_object$proportion_trace)
})
