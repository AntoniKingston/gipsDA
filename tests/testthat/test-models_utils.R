# ==============================================================================
# Global Setup for this file
# ==============================================================================

setup_data <- function() {
  data(iris)
  local_iris <- iris
  local_iris$Species <- as.numeric(local_iris$Species)
  groups <- split(local_iris[, 1:4], local_iris$Species)

  covs <- lapply(groups, cov)
  ns <- sapply(groups, nrow)

  return(list(covs = covs, ns = ns))
}

# ==============================================================================
# GROUP 1: Serialization Tests (JSON Logic)
# ==============================================================================

test_that("recursive_length calculates depth correctly", {
  expect_equal(recursive_length(c(1, 2, 3)), 3)

  nested <- list(a = 1, b = list(c = 2, d = 3), e = 4)
  expect_equal(recursive_length(nested), 4)
  expect_equal(recursive_length(list()), 0)

  func <- function(x) {
    x + 1
  }
  expect_equal(recursive_length(func), 0)
})

test_that("Serialization handles basic types (Matrix, Formula)", {
  original_data <- list(
    mat = matrix(1:4, nrow = 2),
    form = Species ~ Sepal.Length
  )

  serialized <- serialize_for_json(original_data)
  restored <- deserialize_from_json(serialized)

  expect_true(is.matrix(restored$mat))
  expect_equal(restored$mat, original_data$mat)
  expect_equal(format(restored$form), format(original_data$form))
})

test_that("serialize_for_json works with gips_perm class and formulas", {
  fake_perm <- list(1L, 2L, 3L)
  class(fake_perm) <- "gips_perm"

  attr(fake_perm, "size") <- 3L

  result <- serialize_for_json(fake_perm)

  expect_equal(result$`__type`, "gips_perm")
  expect_equal(result$size, 3)

  fake_formula <- list("x", "y")
  class(fake_formula) <- "formula"

  result <- serialize_for_json(fake_formula)

  expect_equal(result$`__type`, "formula")
})

test_that("deserialize_from_json works with gips_perm class and formulas", {
  json_input_perm <- list(
    `__type` = "gips_perm",
    value = "(1)(2)",
    size = 2L
  )

  real_perm_obj <- deserialize_from_json(json_input_perm)
  expect_s3_class(real_perm_obj, "gips_perm")

  json_input_form <- list(
    `__type` = "formula",
    value = "Species ~ ."
  )

  real_form_obj <- deserialize_from_json(json_input_form)
  expect_s3_class(real_form_obj, "formula")
})

test_that("gipsDA_to_json and gipsDA_from_json perform full file round-trip", {
  temp_file <- tempfile(fileext = ".json")
  # Ensure cleanup happens even if test fails
  on.exit(unlink(temp_file))

  complex_obj <- list(
    covs = list(matrix(rnorm(4), 2, 2)),
    model_formula = y ~ x + z,
    meta = list(version = 1.0)
  )
  class(complex_obj) <- "my_gips_DA"

  gipsDA_to_json(complex_obj, temp_file)

  expect_true(file.exists(temp_file))

  restored_obj <- gipsDA_from_json(temp_file, classname = "my_gips_DA")

  expect_s3_class(restored_obj, "my_gips_DA")
  expect_type(restored_obj$covs, "list")
  expect_true(is.matrix(restored_obj$covs[[1]]))
})


# ==============================================================================
# GROUP 2: Projection Logic Tests
# ==============================================================================

test_that("project_matrix_multiperm returns weighted average matrix", {
  emp_cov <- matrix(c(
    4, 2, 0, 1,
    2, 3, 1, 0,
    0, 1, 2, 1,
    1, 0, 1, 5
  ), nrow = 4, byrow = TRUE)
  probs <- c("(1,2,3)" = 0.5, "(2,3,4)" = 0.5)

  res <- project_matrix_multiperm(emp_cov, probs)

  expect_true(is.matrix(res))
  expect_equal(dim(res), dim(emp_cov))
})

test_that("project_covs (MAP=TRUE) returns correct structure using Iris data", {
  d <- setup_data()

  res <- project_covs(d$covs, d$ns, MAP = TRUE, optimizer = "BF", max_iter = 10)

  expect_type(res, "list")
  # expect_named checks names directly
  expect_named(res, c("covs", "opt_info"))

  expect_type(res$covs, "list")
  # expect_length checks list size directly
  expect_length(res$covs, 3)

  expect_equal(dim(res$covs[[1]]), c(4, 4))
  expect_true(is.matrix(res$covs[[1]]))
  expect_false(is.null(res$opt_info))
})

test_that("project_covs (MAP=FALSE) handles probabilities using Iris data", {
  d <- setup_data()

  res <- project_covs(d$covs, d$ns, MAP = FALSE, optimizer = "BF", max_iter = 10)

  expect_type(res, "list")
  expect_named(res, c("covs", "opt_info"))

  # Check if numeric OR list (simple OR logic logic)
  expect_true(is.numeric(res$opt_info) || is.list(res$opt_info))

  expect_length(res$covs, 3)
  expect_equal(dim(res$covs[[2]]), c(4, 4))
})

test_that("project_covs warns if all probabilities are below tolerance", {
  d <- setup_data()

  expect_warning(
    project_covs(d$covs, d$ns, MAP = FALSE, optimizer = "MH", max_iter = 5, tol = 2.0),
    regexp = "There are no perms with estimated probability above threshold"
  )
})
