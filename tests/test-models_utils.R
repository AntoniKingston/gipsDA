# tinytest::run_test_dir("tests")

library(tinytest)

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

# Test: recursive_length calculates depth correctly
expect_equal(recursive_length(c(1, 2, 3)), 3)

nested <- list(a = 1, b = list(c = 2, d = 3), e = 4)
expect_equal(recursive_length(nested), 4)
expect_equal(recursive_length(list()), 0)


# Test: Serialization handles basic types (Matrix, Formula)
original_data <- list(
  mat = matrix(1:4, nrow = 2),
  form = Species ~ Sepal.Length
)

serialized <- serialize_for_json(original_data)
restored <- deserialize_from_json(serialized)

expect_true(is.matrix(restored$mat))
expect_equal(restored$mat, original_data$mat)
expect_equal(format(restored$form), format(original_data$form))


# Test: gipsDA_to_json and gipsDA_from_json perform full file round-trip
temp_file <- tempfile(fileext = ".json")

complex_obj <- list(
  covs = list(matrix(rnorm(4), 2, 2)),
  model_formula = y ~ x + z,
  meta = list(version = 1.0)
)
class(complex_obj) <- "my_gips_DA"

gipsDA_to_json(complex_obj, temp_file)

expect_true(file.exists(temp_file))

restored_obj <- gipsDA_from_json(temp_file, classname = "my_gips_DA")

expect_true(inherits(restored_obj, "my_gips_DA"))
# Sprawdzamy czy to lista (zamiast expect_type)
expect_true(is.list(restored_obj$covs))
expect_true(is.matrix(restored_obj$covs[[1]]))

unlink(temp_file)


# ==============================================================================
# GROUP 2: Projection Logic Tests
# ==============================================================================

# Test: project_matrix_multiperm returns weighted average matrix
emp_cov <- matrix(c(4, 2, 0, 1,
                    2, 3, 1, 0,
                    0, 1, 2, 1,
                    1, 0, 1, 5), nrow = 4, byrow = TRUE)
probs <- c("(1,2,3)" = 0.5, "(2,3,4)" = 0.5)

res <- project_matrix_multiperm(emp_cov, probs)

expect_true(is.matrix(res))
expect_equal(dim(res), dim(emp_cov))


# Test: project_covs (MAP=TRUE) returns correct structure using Iris data
d <- setup_data()

res <- project_covs(d$covs, d$ns, MAP = TRUE, optimizer = "BF", max_iter = 10)

expect_true(is.list(res))
# Sprawdzamy nazwy (zamiast expect_named)
expect_equal(names(res), c("covs", "opt_info"))

expect_true(is.list(res$covs))
# Sprawdzamy długość (zamiast expect_length)
expect_equal(length(res$covs), 3)

expect_equal(dim(res$covs[[1]]), c(4, 4))
expect_true(is.matrix(res$covs[[1]]))
expect_false(is.null(res$opt_info))


# Test: project_covs (MAP=FALSE) handles probabilities using Iris data
d <- setup_data()

res <- project_covs(d$covs, d$ns, MAP = FALSE, optimizer = "BF", max_iter = 10)

expect_true(is.list(res))
expect_equal(names(res), c("covs", "opt_info"))

expect_true(is.numeric(res$opt_info) || is.list(res$opt_info))
expect_equal(length(res$covs), 3)
expect_equal(dim(res$covs[[2]]), c(4, 4))


# Test: project_covs warns if all probabilities are below tolerance
d <- setup_data()

expect_warning(
  project_covs(d$covs, d$ns, MAP = FALSE, optimizer = "MH", max_iter = 5, tol = 2.0),
  pattern = "There are no perms with estimated probability above threshold"
)