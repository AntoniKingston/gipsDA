json_roundtrip <- function(object) {
  jsonlite::unserializeJSON(jsonlite::serializeJSON(object, digits = NA))
}

expect_same_prediction <- function(before, after) {
  expect_identical(after$class, before$class)
  expect_equal(after$posterior, before$posterior, tolerance = 1e-12)
  if (!is.null(before$x)) {
    expect_equal(after$x, before$x, tolerance = 1e-12)
  }
}

test_that("supported model interfaces preserve predictions through JSON", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  new_matrix <- fixture$x[1:6, , drop = FALSE]
  new_frame <- as.data.frame(new_matrix)
  model_funs <- list(
    gipslda = gipsDA::gipslda,
    gipsqda = gipsDA::gipsqda,
    gipsmultqda = gipsDA::gipsmultqda
  )
  interfaces <- list(
    matrix = list(x = fixture$x, newdata = new_matrix),
    data.frame = list(
      x = as.data.frame(fixture$x),
      newdata = new_frame
    )
  )

  for (model_name in names(model_funs)) {
    for (interface_name in names(interfaces)) {
      interface <- interfaces[[interface_name]]
      fit <- suppressWarnings(model_funs[[model_name]](
        interface$x,
        fixture$grouping,
        optimizer = "BF"
      ))
      before <- predict(fit, interface$newdata)
      restored <- json_roundtrip(fit)
      after <- predict(restored, interface$newdata)

      expect_s3_class(restored, model_name)
      expect_identical(restored$lev, fit$lev)
      expect_equal(restored$means, fit$means, tolerance = 1e-12)
      expect_same_prediction(before, after)

      summary_output <- capture.output(returned <- print(restored))
      expect_identical(returned, restored)
      expect_true(any(grepl("^Call:", summary_output)))
      expect_true(any(grepl("Prior probabilities of groups", summary_output)))
      expect_true(any(grepl("Group means", summary_output)))
      expect_true(any(grepl("estimated probabilities", summary_output)))
      if (identical(model_name, "gipslda")) {
        expect_true(any(grepl(
          "Coefficients of linear discriminants",
          summary_output
        )))
      }
    }
  }
})

test_that("the package JSON helpers preserve an LDA prediction workflow", {
  fixture <- make_classification_fixture(p = 2, n_per_class = 8)
  fit <- gipslda(fixture$x, fixture$grouping, optimizer = "BF")
  newdata <- fixture$x[1:6, , drop = FALSE]
  before <- predict(fit, newdata)
  path <- withr::local_tempfile(fileext = ".json")

  gipsDA:::gipsDA_to_json(fit, path)
  restored <- gipsDA:::gipsDA_from_json(path, "gipslda")
  after <- suppressWarnings(predict(restored, newdata))

  expect_s3_class(restored, "gipslda")
  expect_identical(after$class, before$class)
  expect_equal(
    unname(after$posterior),
    unname(before$posterior),
    tolerance = 1e-3
  )
  expect_equal(unname(after$x), unname(before$x), tolerance = 2e-3)
})

test_that("LDA graphics branches write non-empty PDFs with expected dimensions", {
  one_dimension <- make_binary_fixture(p = 1, n_per_class = 8)
  two_dimensions <- make_classification_fixture(p = 2, n_per_class = 8)
  fit_one <- gipslda(
    one_dimension$x,
    one_dimension$grouping,
    optimizer = "BF"
  )
  fit_two <- gipslda(
    two_dimensions$x,
    two_dimensions$grouping,
    optimizer = "BF"
  )

  expect_identical(colnames(fit_one$scaling), "LD1")
  expect_identical(colnames(fit_two$scaling), c("LD1", "LD2"))
  expect_equal(dim(predict(fit_one, one_dimension$x)$x), c(16, 1))
  expect_equal(dim(predict(fit_two, two_dimensions$x)$x), c(24, 2))

  plot_one_pdf <- withr::local_tempfile(fileext = ".pdf")
  withr::with_pdf(
    plot_one_pdf,
    expect_invisible(plot(fit_one, xlab = "LD1"))
  )
  expect_gt(file.info(plot_one_pdf)$size, 0)

  plot_two_pdf <- withr::local_tempfile(fileext = ".pdf")
  withr::with_pdf(
    plot_two_pdf,
    expect_invisible(plot(fit_two, xlab = "LD1", ylab = "LD2"))
  )
  expect_gt(file.info(plot_two_pdf)$size, 0)

  pairs_pdf <- withr::local_tempfile(fileext = ".pdf")
  withr::with_pdf(
    pairs_pdf,
    {
      expect_invisible(pairs(fit_two, labels = c("LD1", "LD2"), type = "std"))
      expect_invisible(pairs(fit_two, type = "trellis"))
    }
  )
  expect_gt(file.info(pairs_pdf)$size, 0)
})
