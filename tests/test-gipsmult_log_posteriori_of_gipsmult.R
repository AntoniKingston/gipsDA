library(tinytest)

# --- Helper: Setup Data ---
S1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
Ss <- list(S1)
ns <- c(100L) # Integer

# ==============================================================================
# 1. Calculation Logic Tests
# ==============================================================================

# Test: log_posteriori_of_gipsmult returns a single numeric value
g <- gipsmult(Ss, ns) # By default was_mean_estimated = TRUE

val <- log_posteriori_of_gipsmult(g)

expect_true(is.numeric(val))
expect_equal(length(val), 1)
expect_true(is.finite(val))


# Test: log_posteriori_of_perm (internal function wrapper check)
D_mats <- attr(g, "D_matrices")
delta <- attr(g, "delta")
perm <- g[[1]]

# FIX: Because gipsmult() sets was_mean_estimated = TRUE by default,
# the log_posteriori_of_gipsmult function internally subtracts 1 from ns.
# To match the result manually, we must also subtract 1 here.
val_internal <- log_posteriori_of_perm(perm, Ss, ns - 1L, delta, D_mats)

expect_true(is.numeric(val_internal))
expect_equal(val, val_internal)