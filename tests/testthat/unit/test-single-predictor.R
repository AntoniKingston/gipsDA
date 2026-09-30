test_that("one-predictor LDA and QDA store the identity probability", {
  x <- matrix(seq_len(20), ncol = 1, dimnames = list(NULL, "x"))
  grouping <- factor(rep(c("a", "b"), each = 10))

  for (fitter in list(gipslda, gipsqda)) {
    for (optimizer in c("BF", "MH")) {
      fit <- fitter(x, grouping, MAP = TRUE, optimizer = optimizer, max_iter = 10)
      probabilities <- fit$optimization_info
      if (inherits(fit, "gipsqda")) {
        expect_length(probabilities, 2)
        for (prob in probabilities) expect_identical(prob, c("()" = 1))
      } else {
        expect_identical(probabilities, c("()" = 1))
      }
      expect_true(all(is.finite(predict(fit, x)$posterior)))
      expect_error(fitter(x, grouping, MAP = FALSE),
                   "with one predictor requires MAP = TRUE", fixed = TRUE)
    }
    fit <- fitter(x, grouping, MAP = TRUE, store_probabilities = FALSE)
    if (inherits(fit, "gipsqda")) {
      expect_true(all(vapply(fit$optimization_info, is.null, logical(1))))
    } else {
      expect_null(fit$optimization_info)
    }
    data <- data.frame(x = x[, 1], group = grouping)
    fit <- fitter(group ~ x, data = data)
    expect_true(all(is.finite(predict(fit, data)$posterior)))
    expect_error(fitter(group ~ x, data = data, MAP = FALSE),
                 "with one predictor requires MAP = TRUE", fixed = TRUE)
  }
})

test_that("multivariate QDA reports its minimum predictor count", {
  x <- matrix(seq_len(20), ncol = 1)
  grouping <- factor(rep(c("a", "b"), each = 10))
  for (map in c(TRUE, FALSE)) {
    expect_error(gipsmultqda(x, grouping, MAP = map),
                 "gipsmultqda requires at least two predictors", fixed = TRUE)
  }
})

test_that("QDA reports unused grouping levels before covariance estimation", {
  x <- cbind(seq_len(20), sin(seq_len(20)))
  grouping <- factor(rep(c("a", "b"), each = 10),
                     levels = c("a", "b", "unused"))
  for (fitter in list(gipsqda, gipsmultqda)) {
    for (map in c(TRUE, FALSE)) {
      expect_error(fitter(x, grouping, MAP = map),
                   "Unused levels in grouping: unused. Remove them with droplevels(grouping)",
                   fixed = TRUE)
      fit <- fitter(x, droplevels(grouping), MAP = map)
      expect_identical(fit$lev, c("a", "b"))
      expect_true(all(is.finite(predict(fit, x)$posterior)))
    }
  }
})
