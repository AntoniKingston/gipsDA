setup_data <- function() {
  data(iris)
  iris$Species <- as.numeric(iris$Species)
  groups <- split(iris[, 1:4], iris$Species)

  covs <- lapply(groups, cov)
  ns <- sapply(groups, nrow)

  return(list(covs = covs, ns = ns))
}

# GROUP 1: Serialization Tests (JSON Logic)

testthat::test_that("recursive_length calculates depth correctly", {
  # Simple atomic vector
  testthat::expect_equal(recursive_length(c(1, 2, 3)), 3)

  # Nested list
  nested <- list(a = 1, b = list(c = 2, d = 3), e = 4)
  testthat::expect_equal(recursive_length(nested), 4)

  # Empty list
  testthat::expect_equal(recursive_length(list()), 0)
})

testthat::test_that("Serialization handles basic types (Matrix, Formula)", {
  # Arrange
  original_data <- list(
    mat = matrix(1:4, nrow = 2),
    form = Species ~ Sepal.Length
  )

  # Act: Serialize then Deserialize
  serialized <- serialize_for_json(original_data)
  restored <- deserialize_from_json(serialized)

  # Assert
  testthat::expect_true(is.matrix(restored$mat))
  testthat::expect_equal(restored$mat, original_data$mat)

  testthat::expect_s3_class(restored$form, "formula")
  testthat::expect_equal(format(restored$form), format(original_data$form))
})

testthat::test_that("gipsDA_to_json and gipsDA_from_json perform full file round-trip", {
  # Arrange
  temp_file <- tempfile(fileext = ".json")

  # Create a complex object mimicking your DA structure
  complex_obj <- list(
    covs = list(matrix(rnorm(4), 2, 2)),
    model_formula = y ~ x + z,
    meta = list(version = 1.0)
  )
  class(complex_obj) <- "my_gips_DA"

  # Act
  gipsDA_to_json(complex_obj, temp_file)

  # Check if file exists
  testthat::expect_true(file.exists(temp_file))

  # Load back
  restored_obj <- gipsDA_from_json(temp_file, classname = "my_gips_DA")

  # Assert
  testthat::expect_s3_class(restored_obj, "my_gips_DA")
  testthat::expect_type(restored_obj$covs, "list")
  testthat::expect_true(is.matrix(restored_obj$covs[[1]]))
  testthat::expect_s3_class(restored_obj$model_formula, "formula")

  # Cleanup
  unlink(temp_file)
})

# ==============================================================================
# GROUP 2: Projection Logic Tests
# ==============================================================================

testthat::test_that("project_matrix_multiperm returns weighted average matrix", {
  # Arrange
  # We assume gips::project_matrix(mat, "perm") works.
  # Using a 4x4 matrix to match Iris dimensions style, though logic is generic.
  emp_cov <- diag(4)

  # Mock probabilities: 50% perm1, 50% perm2 (Identity permutation "()")
  probs <- list("()" = 0.5, "()" = 0.5)

  # Act
  res <- project_matrix_multiperm(emp_cov, probs)

  # Assert
  testthat::expect_true(is.matrix(res))
  testthat::expect_equal(dim(res), dim(emp_cov))

  # Identity projected onto Identity is Identity. Average is still Identity.
  testthat::expect_equal(res, emp_cov)
})

testthat::test_that("project_covs (MAP=TRUE) returns correct structure using Iris data", {
  # Arrange
  d <- setup_data()

  # Act
  # Assuming 'gipsmult' and 'find_MAP' are available.
  res <- project_covs(d$covs, d$ns, MAP = TRUE, optimizer = "hill_climbing", max_iter = 10)

  # Assert
  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("covs", "opt_info"))

  # Check covs
  testthat::expect_type(res$covs, "list")

  # Iris has 3 classes, so we expect 3 covariance matrices
  testthat::expect_length(res$covs, 3)

  # Iris has 4 variables, so matrices should be 4x4
  testthat::expect_equal(dim(res$covs[[1]]), c(4, 4))
  testthat::expect_true(is.matrix(res$covs[[1]]))

  # Check opt_info
  testthat::expect_false(is.null(res$opt_info))
})

testthat::test_that("project_covs (MAP=FALSE) handles probabilities using Iris data", {
  # Arrange
  d <- setup_data()

  # Act
  # This path triggers find_MAP with return_probabilities=TRUE
  # We use 'Metropolis_Hastings' as it is typical for BMA/probabilities
  res <- project_covs(d$covs, d$ns, MAP = FALSE, optimizer = "Metropolis_Hastings", max_iter = 10)

  # Assert
  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("covs", "opt_info"))

  # opt_info should be the probabilities vector/list
  testthat::expect_true(is.numeric(res$opt_info) || is.list(res$opt_info))

  # Check dimensions again to be sure
  testthat::expect_length(res$covs, 3)
  testthat::expect_equal(dim(res$covs[[2]]), c(4, 4))
})

testthat::test_that("project_covs warns if all probabilities are below tolerance", {
  d <- setup_data()

  # Set tolerance > 1.0 to force the warning (since sum of probs is ~1.0)
  testthat::expect_warning(
    project_covs(d$covs, d$ns, MAP = FALSE, optimizer = "Metropolis_Hastings", max_iter = 5, tol = 2.0),
    "There are no perms with estimated probability above threshold"
  )
})