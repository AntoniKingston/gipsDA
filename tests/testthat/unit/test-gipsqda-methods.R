qda_method_specs <- list(
  gipsqda = list(
    generic = gipsDA::gipsqda,
    predict = gipsDA:::predict.gipsqda,
    class = "gipsqda"
  ),
  gipsmultqda = list(
    generic = gipsDA::gipsmultqda,
    predict = gipsDA:::predict.gipsmultqda,
    class = "gipsmultqda"
  )
)

fit_qda_method <- function(spec, ...) {
  suppressWarnings(do.call(spec$generic, list(...)))
}

expect_valid_qda_prediction <- function(prediction, n, groups) {
  expect_named(prediction, c("class", "posterior"))
  expect_s3_class(prediction$class, "factor")
  expect_length(prediction$class, n)
  expect_valid_posterior(prediction$posterior, n, length(groups))
  expect_equal(levels(prediction$class), groups)
  expect_equal(colnames(prediction$posterior), groups)
}

test_that("all QDA prediction rules return normalized posteriors", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 7)
  methods <- c("plug-in", "predictive", "debiased")

  for (spec in qda_method_specs) {
    fit <- suppressWarnings(spec$generic(
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    ))

    for (method in methods) {
      prediction <- predict(
        fit,
        fixture$x[1:5, , drop = FALSE],
        method = method
      )
      expect_valid_qda_prediction(prediction, 5, levels(fixture$grouping))
    }

    loo <- predict(fit, method = "looCV")
    expect_valid_qda_prediction(
      loo,
      nrow(fixture$x),
      levels(fixture$grouping)
    )
  }
})

test_that("formula QDA predicts from supplied model frames", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 6)

  for (spec in qda_method_specs) {
    fit <- fit_qda_method(
      spec,
      class ~ x1 + x2,
      data = fixture$data,
      optimizer = "BF"
    )
    supplied <- predict(fit, fixture$data[1:4, , drop = FALSE])

    expect_valid_qda_prediction(supplied, 4, levels(fixture$grouping))
  }
})

for (name in names(qda_method_specs)) {
  local({
    current_name <- name

    test_that(paste(current_name, "formula prediction reconstructs missing newdata"), {
      fixture <- make_classification_fixture(p = 2, n_per_class = 6)
      spec <- qda_method_specs[[current_name]]
      fit <- suppressWarnings(spec$generic(
        class ~ x1 + x2,
        data = fixture$data,
        optimizer = "BF"
      ))

      prediction <- predict(fit)
      expect_valid_qda_prediction(
        prediction,
        nrow(fixture$data),
        levels(fixture$grouping)
      )
      expect_equal(rownames(prediction$posterior), rownames(fixture$data))
    })

    test_that(paste(current_name, "formula prediction supports looCV"), {
      fixture <- make_classification_fixture(p = 2, n_per_class = 6)
      spec <- qda_method_specs[[current_name]]
      fit <- suppressWarnings(spec$generic(
        class ~ x1 + x2,
        data = fixture$data,
        optimizer = "BF"
      ))

      prediction <- predict(fit, method = "looCV")
      expect_valid_qda_prediction(
        prediction,
        nrow(fixture$data),
        levels(fixture$grouping)
      )
    })
  })
}

test_that("matrix QDA handles missing, matrix, data-frame, and single-row newdata", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  for (spec in qda_method_specs) {
    fit <- suppressWarnings(spec$generic(
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    ))

    missing_newdata <- predict(fit)
    matrix_newdata <- predict(fit, fixture$x[1:3, , drop = FALSE])
    frame_newdata <- predict(fit, as.data.frame(fixture$x[1:3, , drop = FALSE]))
    row_newdata <- predict(fit, fixture$x[1, ])

    expect_valid_qda_prediction(
      missing_newdata,
      nrow(fixture$x),
      levels(fixture$grouping)
    )
    expect_valid_qda_prediction(
      matrix_newdata,
      3,
      levels(fixture$grouping)
    )
    expect_valid_qda_prediction(
      frame_newdata,
      3,
      levels(fixture$grouping)
    )
    expect_valid_qda_prediction(
      row_newdata,
      1,
      levels(fixture$grouping)
    )
    expect_equal(matrix_newdata$posterior, frame_newdata$posterior)
  }
})

test_that("prediction honors supplied priors and validates bad priors", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  for (spec in qda_method_specs) {
    fit <- fit_qda_method(
      spec,
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    )
    balanced <- predict(fit, fixture$x[1:3, , drop = FALSE])
    shifted <- predict(
      fit,
      fixture$x[1:3, , drop = FALSE],
      prior = c(0.99, 0.01)
    )

    expect_false(isTRUE(all.equal(balanced$posterior, shifted$posterior)))
    expect_error(
      predict(
        fit,
        fixture$x[1:2, , drop = FALSE],
        prior = c(-0.1, 1.1)
      ),
      "invalid 'prior'",
      fixed = TRUE
    )
    expect_error(
      predict(
        fit,
        fixture$x[1:2, , drop = FALSE],
        prior = c(0.2, 0.2)
      ),
      "invalid 'prior'",
      fixed = TRUE
    )
    expect_error(
      predict(
        fit,
        fixture$x[1:2, , drop = FALSE],
        prior = c(0.2, 0.3, 0.5)
      ),
      "'prior' is of incorrect length",
      fixed = TRUE
    )
  }
})

test_that("prediction validates object classes, methods, and newdata shape", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  for (spec in qda_method_specs) {
    fit <- fit_qda_method(
      spec,
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    )

    expect_error(
      spec$predict(unclass(fit), fixture$x[1:2, , drop = FALSE]),
      paste0("object not of class \"", spec$class, "\""),
      fixed = TRUE
    )
    expect_error(
      predict(fit, fixture$x[1:2, , drop = FALSE], method = "unknown"),
      "'arg' should be one of"
    )
    expect_error(
      predict(fit, fixture$x[1:2, 1, drop = FALSE]),
      "wrong number of variables",
      fixed = TRUE
    )
    expect_error(
      predict(fit, matrix(letters[1:4], ncol = 2)),
      "non-numeric argument|requires numeric"
    )
    expect_error(
      predict(fit, fixture$x[1:2, , drop = FALSE], method = "looCV"),
      "cannot have leave-one-out CV with 'newdata'",
      fixed = TRUE
    )

    incompatible <- fit
    incompatible$call$method <- "robust"
    expect_error(
      predict(incompatible, method = "looCV"),
      "cannot use leave-one-out CV with method"
    )
  }
})

test_that("prediction warns when variable names disagree", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)
  renamed <- fixture$x[1:2, , drop = FALSE]
  colnames(renamed) <- c("wrong_1", "wrong_2")

  for (spec in qda_method_specs) {
    fit <- fit_qda_method(
      spec,
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    )
    expect_warning(
      prediction <- predict(fit, renamed),
      "variable names in 'newdata' do not match"
    )
    expect_valid_qda_prediction(
      prediction,
      2,
      levels(fixture$grouping)
    )
  }
})

test_that("formula prediction rejects missing predictors and handles coercion", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 6)

  for (spec in qda_method_specs) {
    fit <- fit_qda_method(
      spec,
      class ~ x1 + x2,
      data = fixture$data,
      optimizer = "BF"
    )

    expect_error(
      predict(fit, fixture$data[1:2, "x1", drop = FALSE]),
      "x2"
    )

    wrong_class <- fixture$data[1:2, , drop = FALSE]
    wrong_class$x1 <- as.character(wrong_class$x1)
    coerced <- suppressWarnings(predict(fit, wrong_class))
    expect_valid_qda_prediction(coerced, 2, levels(fixture$grouping))
  }
})

test_that("print methods report model summaries and return invisibly", {
  fixture <- make_binary_fixture(p = 2, n_per_class = 5)

  for (spec in qda_method_specs) {
    fit <- fit_qda_method(
      spec,
      fixture$x,
      fixture$grouping,
      optimizer = "BF"
    )

    expect_output(
      returned <- print(fit),
      "Prior probabilities of groups"
    )
    expect_output(print(fit), "Group means")
    expect_output(print(fit), "Posterior probabilities of retained permutations")
    expect_identical(returned, fit)
  }
})

for (name in names(qda_method_specs)) {
  local({
    current_name <- name

    test_that(paste("model.frame reconstructs", current_name, "formula data"), {
      fixture <- make_binary_fixture(p = 2, n_per_class = 6)
      spec <- qda_method_specs[[current_name]]
      fit <- suppressWarnings(spec$generic(
        class ~ x1 + x2,
        data = fixture$data,
        subset = 2:11,
        optimizer = "BF"
      ))
      reconstructed <- model.frame(fit)

      expect_s3_class(reconstructed, "data.frame")
      expect_equal(names(reconstructed), c("class", "x1", "x2"))
      expect_equal(nrow(reconstructed), 10)
      expect_equal(reconstructed$class, fixture$data$class[2:11])
      expect_equal(
        as.matrix(reconstructed[c("x1", "x2")]),
        as.matrix(fixture$data[2:11, c("x1", "x2")])
      )
    })
  })
}
