library(tinytest)

# --- Data Preparation ---
# We use the standard Iris dataset for testing
data(iris)

# Split into Matrix (X) and Grouping Factor (Y) for testing matrix interfaces
X <- as.matrix(iris[, 1:4])
Y <- iris$Species

# ==============================================================================
# 1. Testing Model Fitting Interfaces (S3 Dispatch)
# ==============================================================================

# Test: Does gipslda work with the formula interface?
fit_formula <- gipslda(Species ~ ., data = iris)

# Check if the returned object has the correct class
expect_true(inherits(fit_formula, "gipslda"))
expect_true(is.list(fit_formula))

# Check if the result contains key components (similar to MASS::lda)
required_components <- c("prior", "means", "scaling", "counts")
expect_true(all(required_components %in% names(fit_formula)))

# Check dimensions of the means (3 classes x 4 variables)
expect_equal(dim(fit_formula$means), c(3, 4))


# Test: Does gipslda work with the matrix/default interface?
fit_matrix <- gipslda(x = X, grouping = Y)

expect_true(inherits(fit_matrix, "gipslda"))

# The results between formula and matrix interface should be consistent
# (e.g., priors and counts should be identical)
expect_equal(fit_matrix$prior, fit_formula$prior)
expect_equal(fit_matrix$counts, fit_formula$counts)


# Test: Does gipslda work with the data.frame interface?
# (Assuming gipslda.data.frame passes data to gipslda.default)
fit_df <- gipslda(x = iris[, 1:4], grouping = Y)
expect_true(inherits(fit_df, "gipslda"))


# Test: Checking specific arguments (MAP, weighted_avg)
# We verify that the function runs without errors when custom arguments are passed
fit_opts <- gipslda(X, Y, MAP = FALSE, weighted_avg = TRUE)
expect_true(inherits(fit_opts, "gipslda"))


# ==============================================================================
# 2. Testing the 'print' method
# ==============================================================================

# Test: Does print.gipslda execute without errors?
# We capture the output to avoid cluttering the console during testing
output <- capture.output(print(fit_formula))

# Expect that it prints at least something (output length > 0)
expect_true(length(output) > 0)


# ==============================================================================
# 3. Testing the 'predict' method
# ==============================================================================

# Test: Basic prediction on training data
pred <- predict(fit_formula, newdata = iris)

expect_true(is.list(pred))
# Check if it returns the standard LDA prediction components
expect_equal(names(pred), c("class", "posterior", "x"))

# Check dimensions matching the input data
expect_equal(length(pred$class), nrow(iris))
expect_equal(nrow(pred$posterior), nrow(iris))
expect_equal(ncol(pred$posterior), 3) # 3 classes

# Check if posterior probabilities sum to 1 (row-wise)
row_sums <- rowSums(pred$posterior)
# We allow for a tiny floating-point tolerance
expect_true(all(abs(row_sums - 1) < 1e-6))


# Test: Prediction on a new, smaller subset of data
new_data <- iris[1:5, ]
pred_small <- predict(fit_formula, newdata = new_data)
expect_equal(length(pred_small$class), 5)