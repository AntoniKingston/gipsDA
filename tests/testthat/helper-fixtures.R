make_classification_fixture <- function(p = 2L, n_per_class = 12L) {
  stopifnot(p >= 1L, p <= 4L, n_per_class >= 4L, n_per_class <= 50L)

  rows <- seq_len(n_per_class)
  x <- rbind(
    iris3[rows, seq_len(p), "Setosa", drop = FALSE],
    iris3[rows, seq_len(p), "Versicolor", drop = FALSE],
    iris3[rows, seq_len(p), "Virginica", drop = FALSE]
  )
  dim(x) <- c(3L * n_per_class, p)
  colnames(x) <- paste0("x", seq_len(p))
  grouping <- factor(rep(c("setosa", "versicolor", "virginica"),
    each = n_per_class
  ))

  data <- data.frame(class = grouping, x, check.names = FALSE)
  list(x = x, grouping = grouping, data = data)
}

make_binary_fixture <- function(p = 2L, n_per_class = 12L) {
  fixture <- make_classification_fixture(p, n_per_class)
  keep <- fixture$grouping != "virginica"
  fixture$x <- fixture$x[keep, , drop = FALSE]
  fixture$grouping <- droplevels(fixture$grouping[keep])
  fixture$data <- droplevels(fixture$data[keep, , drop = FALSE])
  fixture
}

expect_valid_posterior <- function(posterior, n, groups) {
  if (!is.numeric(groups) || length(groups) != 1L) {
    groups <- length(groups)
  }
  expect_true(is.matrix(posterior))
  expect_equal(nrow(posterior), n)
  expect_equal(ncol(posterior), groups)
  expect_true(all(is.finite(posterior)))
  expect_true(all(posterior >= 0))
  expect_equal(unname(rowSums(posterior)), rep(1, n), tolerance = 1e-8)
}

expect_valid_lda_fit <- function(object, n, p, groups) {
  expect_s3_class(object, "gipslda")
  expect_equal(object$N, n)
  expect_equal(length(object$counts), groups)
  expect_equal(dim(object$means), c(groups, p))
  expect_equal(nrow(object$scaling), p)
  expect_equal(sum(object$prior), 1, tolerance = 1e-10)
  expect_true(all(is.finite(object$scaling)))
}

expect_valid_qda_fit <- function(object, class, n, p, groups) {
  expect_s3_class(object, class)
  expect_equal(object$N, n)
  expect_equal(length(object$counts), groups)
  expect_equal(dim(object$means), c(groups, p))
  expect_equal(dim(object$scaling), c(p, p, groups))
  expect_length(object$ldet, groups)
  expect_equal(sum(object$prior), 1, tolerance = 1e-10)
  expect_true(all(is.finite(object$scaling)))
}

local_plot_device <- function() {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  withr::defer({
    grDevices::dev.off()
    unlink(path)
  }, envir = parent.frame())
  invisible(path)
}
