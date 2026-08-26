test_that("project_matrix_multiperm returns a normalized weighted projection", {
  emp_cov <- matrix(c(
    2, 0.4,
    0.4, 1
  ), nrow = 2, byrow = TRUE)
  probs <- c("(1,2)" = 1, "()" = 3)

  expected <- (
    probs[[1]] * gips::project_matrix(emp_cov, names(probs)[[1]]) +
      probs[[2]] * gips::project_matrix(emp_cov, names(probs)[[2]])
  ) / sum(probs)

  expect_equal(
    gipsDA:::project_matrix_multiperm(emp_cov, probs),
    expected
  )
})

test_that("project_covs MAP mode accepts a matrix and a list", {
  cov_1 <- matrix(c(
    2, 0.4,
    0.4, 1
  ), nrow = 2, byrow = TRUE)
  cov_2 <- matrix(c(
    1.5, 0.2,
    0.2, 0.8
  ), nrow = 2, byrow = TRUE)

  matrix_result <- gipsDA:::project_covs(
    cov_1,
    ns_obs = 20,
    MAP = TRUE,
    optimizer = "BF",
    max_iter = NA
  )

  expect_length(matrix_result$covs, 1)
  expect_s3_class(matrix_result$permutation, "gips_perm")
  expect_type(matrix_result$opt_info, "double")
  expect_false(is.null(names(matrix_result$opt_info)))

  expect_equal(
    matrix_result$covs[[1]],
    gips::project_matrix(cov_1, matrix_result$permutation)
  )

  list_result <- gipsDA:::project_covs(
    list(cov_1, cov_2),
    ns_obs = c(20, 24),
    MAP = TRUE,
    optimizer = "BF",
    max_iter = NA
  )

  expect_length(list_result$covs, 2)
  expect_s3_class(list_result$permutation, "gips_perm")
  expect_type(list_result$opt_info, "double")
  expect_false(is.null(names(list_result$opt_info)))

  expect_equal(
    list_result$covs,
    lapply(list(cov_1, cov_2), gips::project_matrix, list_result$permutation)
  )
})

test_that("project_covs non-MAP mode filters BF probabilities", {
  emp_cov <- matrix(c(
    2, 0.4,
    0.4, 1
  ), nrow = 2, byrow = TRUE)

  all_result <- gipsDA:::project_covs(
    emp_cov,
    ns_obs = 20,
    MAP = FALSE,
    optimizer = "BF",
    max_iter = NA,
    tol = 0
  )

  expect_gt(length(all_result$opt_info), 1)
  expect_equal(sum(all_result$opt_info), 1)
  expect_equal(
    all_result$covs[[1]],
    gipsDA:::project_matrix_multiperm(emp_cov, all_result$opt_info)
  )

  threshold <- mean(range(all_result$opt_info))
  filtered_result <- gipsDA:::project_covs(
    emp_cov,
    ns_obs = 20,
    MAP = FALSE,
    optimizer = "BF",
    max_iter = NA,
    tol = threshold
  )
  expected_probs <- all_result$opt_info[all_result$opt_info > threshold]

  expect_equal(filtered_result$opt_info, expected_probs)
  expect_equal(
    filtered_result$covs[[1]],
    gipsDA:::project_matrix_multiperm(emp_cov, expected_probs)
  )
})

test_that("project_covs non-MAP mode falls back to the MAP projection", {
  emp_cov <- matrix(c(
    2, 0.4,
    0.4, 1
  ), nrow = 2, byrow = TRUE)
  all_result <- gipsDA:::project_covs(
    emp_cov,
    ns_obs = 20,
    MAP = FALSE,
    optimizer = "BF",
    max_iter = NA,
    tol = 0
  )

  expect_warning(
    fallback <- gipsDA:::project_covs(
      emp_cov,
      ns_obs = 20,
      MAP = FALSE,
      optimizer = "BF",
      max_iter = NA,
      tol = 1
    ),
    "no perms with estimated probability above threshold"
  )

  expect_equal(fallback$opt_info, all_result$opt_info[1])
  expect_equal(
    fallback$covs[[1]],
    gips::project_matrix(emp_cov, names(all_result$opt_info)[1])
  )
})

test_that("serialize and deserialize preserve calls, matrices, and nested lists", {
  model_call <- quote(stats::lm(y ~ x, data = dat))
  model_matrix <- matrix(1:6, nrow = 2, dimnames = list(NULL, letters[1:3]))
  value <- list(
    call = model_call,
    matrix = model_matrix,
    nested = list(number = 2.5, text = "value", flag = TRUE)
  )

  encoded_call <- gipsDA:::serialize_for_json(model_call)
  expect_identical(encoded_call$`__type`, "call")
  expect_identical(
    gipsDA:::deserialize_from_json(encoded_call),
    model_call
  )

  encoded_matrix <- gipsDA:::serialize_for_json(model_matrix)
  expect_identical(encoded_matrix$`__type`, "matrix")
  expect_identical(encoded_matrix$nrow, 2L)
  expect_identical(encoded_matrix$ncol, 3L)
  expect_equal(
    gipsDA:::deserialize_from_json(encoded_matrix),
    unname(model_matrix)
  )

  expect_equal(
    gipsDA:::deserialize_from_json(gipsDA:::serialize_for_json(value)),
    list(
      call = model_call,
      matrix = unname(model_matrix),
      nested = value$nested
    )
  )
})

test_that("formula and terms serialization restores formulas", {
  model_formula <- y ~ x + z
  model_terms <- stats::terms(model_formula)

  encoded_formula <- gipsDA:::serialize_for_json(model_formula)
  expect_identical(encoded_formula$`__type`, "formula")
  expect_s3_class(
    gipsDA:::deserialize_from_json(encoded_formula),
    "formula"
  )
  expect_equal(
    gipsDA:::deserialize_from_json(encoded_formula),
    model_formula,
    ignore_attr = TRUE
  )

  encoded_terms <- gipsDA:::serialize_for_json(model_terms)
  expect_identical(encoded_terms$`__type`, "formula")
  expect_s3_class(
    gipsDA:::deserialize_from_json(encoded_terms),
    "formula"
  )
  expect_equal(
    gipsDA:::deserialize_from_json(encoded_terms),
    model_formula,
    ignore_attr = TRUE
  )
})

test_that("gips permutation serialization preserves its declared size", {
  perm <- gips::gips_perm("(1,2)", size = 3)

  encoded <- gipsDA:::serialize_for_json(perm)
  expect_identical(encoded$`__type`, "gips_perm")
  expect_identical(encoded$value, as.character(perm))
  expect_equal(encoded$size, attr(perm, "size"))

  restored <- gipsDA:::deserialize_from_json(encoded)
  expect_s3_class(restored, "gips_perm")
  expect_identical(as.character(restored), as.character(perm))
  expect_equal(attr(restored, "size"), attr(perm, "size"))
})

test_that("named atomic vectors preserve names and storage mode", {
  values <- list(
    integer = c(first = 1L, second = 2L),
    double = c(first = 0.25, second = 0.75),
    logical = c(first = TRUE, second = FALSE),
    character = c(first = "a", second = "b")
  )

  for (value in values) {
    encoded <- gipsDA:::serialize_for_json(value)
    restored <- gipsDA:::deserialize_from_json(encoded)

    expect_identical(encoded$`__type`, "named_atomic")
    expect_identical(restored, value)
  }
})

test_that("JSON file roundtrip restores nested supported types and class", {
  file <- tempfile(fileext = ".json")
  withr::defer(unlink(file))
  perm <- gips::gips_perm("(1,2,3)", size = 3)
  object <- structure(
    list(
      matrix = matrix(c(2, 0.5, 0.5, 1), nrow = 2),
      call = quote(stats::predict(fit, newdata = new_data)),
      perm = perm,
      nested = list(label = "fixture", count = 4)
    ),
    class = "utils_fixture"
  )

  gipsDA:::gipsDA_to_json(object, file)
  restored <- gipsDA:::gipsDA_from_json(file, "utils_fixture")

  expect_s3_class(restored, "utils_fixture")
  expect_equal(restored$matrix, object$matrix)
  expect_identical(restored$call, object$call)
  expect_identical(as.character(restored$perm), as.character(object$perm))
  expect_equal(attr(restored$perm, "size"), attr(object$perm, "size"))
  expect_equal(restored$nested, object$nested)
})

test_that("recursive_length counts atomic leaves recursively", {
  expect_identical(gipsDA:::recursive_length(numeric()), 0L)
  expect_identical(gipsDA:::recursive_length(1:4), 4L)
  expect_identical(
    gipsDA:::recursive_length(list(1:2, list("a", logical(0)), NULL)),
    3L
  )
  expect_identical(gipsDA:::recursive_length(globalenv()), 0L)
})

test_that("desingularize leaves sufficiently regular matrices unchanged", {
  A <- matrix(c(
    2, 0.2,
    0.2, 1
  ), nrow = 2, byrow = TRUE)

  expect_identical(gipsDA:::desingularize(A, target = 0.05), A)
})

test_that("desingularize raises small positive and negative eigenvalues", {
  positive <- diag(c(0.01, 2))
  adjusted_positive <- gipsDA:::desingularize(positive, target = 0.05)
  expect_equal(
    min(abs(eigen(adjusted_positive, symmetric = TRUE)$values)),
    0.05,
    tolerance = 1e-12
  )

  negative <- diag(c(-0.01, 2))
  adjusted_negative <- gipsDA:::desingularize(negative, target = 0.05)
  expect_equal(
    min(abs(eigen(adjusted_negative, symmetric = TRUE)$values)),
    0.05,
    tolerance = 1e-12
  )
})

test_that("desingularize supports nonsymmetric square matrices", {
  A <- matrix(c(
    0, 1,
    0, 2
  ), nrow = 2, byrow = TRUE)

  adjusted <- gipsDA:::desingularize(A, target = 0.1)
  expect_equal(
    min(abs(eigen(adjusted, symmetric = FALSE)$values)),
    0.1,
    tolerance = 1e-12
  )
})

test_that("desingularize reports invalid scaling and nonsquare inputs", {
  expect_error(
    gipsDA:::desingularize(diag(c(0.5, 3)), target = 2),
    "Invalid scaling: 1 \\+ s <= 0",
    fixed = FALSE
  )
  expect_error(
    gipsDA:::desingularize(matrix(1:6, nrow = 2)),
    "non-square matrix"
  )
})
