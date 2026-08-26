test_that("the public API exports the three discriminant-analysis generics", {
  exports <- getNamespaceExports("gipsDA")

  expect_setequal(exports, c("gipslda", "gipsqda", "gipsmultqda"))
  expect_true(all(vapply(
    mget(exports, envir = asNamespace("gipsDA")),
    is.function,
    logical(1)
  )))
})

test_that("documented S3 methods are registered", {
  fit_classes <- c("gipslda", "gipsqda", "gipsmultqda")
  fit_inputs <- c("default", "matrix", "data.frame", "formula")

  for (generic in fit_classes) {
    for (input in fit_inputs) {
      expect_true(is.function(getS3method(generic, input)))
    }
    expect_true(is.function(getS3method("predict", generic)))
    expect_true(is.function(getS3method("print", generic)))
    expect_true(is.function(getS3method("model.frame", generic)))
  }

  expect_true(is.function(getS3method("coef", "gipslda")))
  expect_true(is.function(getS3method("plot", "gipslda")))
  expect_true(is.function(getS3method("pairs", "gipslda")))
})

test_that("default optimizer uses BF for p <= 10 and MH for p > 10", {
  seen_optimizer <- character()

  fake_project_covs <- function(emp_covs, ns_obs, MAP, optimizer, max_iter, ...) {
    seen_optimizer <<- c(seen_optimizer, optimizer)
    list(covs = emp_covs, opt_info = "mock optimization")
  }

  testthat::local_mocked_bindings(
    project_covs = fake_project_covs,
    .package = "gipsDA"
  )

  set.seed(1)
  x10 <- matrix(rnorm(60 * 10), ncol = 10)
  y10 <- factor(rep(c("a", "b", "c"), each = 20))

  x11 <- matrix(rnorm(66 * 11), ncol = 11)
  y11 <- factor(rep(c("a", "b", "c"), each = 22))

  suppressWarnings(gipslda(x10, y10))
  suppressWarnings(gipslda(x11, y11))

  expect_identical(seen_optimizer[1], "BF")
  expect_identical(seen_optimizer[2], "MH")
})
