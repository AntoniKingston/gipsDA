# =======================================================================
# SCENARIO 2: gipsMultQDA MODEL
# =======================================================================
# - SAME permutation for all classes.
# - DIFFERENT covariance matrices for all classes.
# - DIFFERENT means for all classes.
# =======================================================================

cat("\n\n--- RUNNING SCENARIO 2: gipsMultQDA ---\n")

# Data parameters
p <- 8          # Number of features (data dimension)
n_classes <- 6  # Number of classes
n_per_class <- 50 # Number of observations per class

# Wishart distribution parameters
df <- p + 2     # Degrees of freedom (must be >= p)
psi <- 100.0       # Scaling factor for covariance matrices

# --- Search parameters ---
max_iterations <- 40           # Safety limit for the search loop
target_train_accuracy <- 0.70  # Target accuracy for the training set
target_test_accuracy  <- 0.50  # Target accuracy for the test set
output_filename <- "scenario_2.csv" # Filename for the output data

# Initialize variables to store results
found_divisor <- FALSE
final_divisor <- NA
final_data_to_save <- NULL
final_lda_train_accuracy <- NA
final_lda_test_accuracy <- NA
final_qda_train_accuracy <- NA
final_qda_test_accuracy <- NA

perm <- "(2,3,4)(1,5,7)"

#-----------------------------------------------------------------------
#               ** AUTOMATED SEARCH LOOP FOR SCALING FACTOR **
#-----------------------------------------------------------------------
cat("Starting search for a suitable scaling divisor...\n")
cat(sprintf("Target LDA Accuracy -> Train: %.2f | Test: %.2f \n\n",
            target_train_accuracy, target_test_accuracy))

for (i in 1:max_iterations) {
  
  cat(sprintf("--- Iteration %d: Testing with divisor psi = %d ---\n", i, psi / 2**i))
  
  # ** Generate synthetic data with the new covariance matrix **
  class_means <- matrix(runif(n_classes * p, min = 0, max = 1), nrow = p, ncol = n_classes)
  data_list <- list()
  
  for (k in 1:n_classes) {
    # Generate a different covariance matrix for each class
    A <- matrix(rnorm(p * p), nrow = p, ncol = p)
    Sigma <- t(A) %*% A
    wishart_output <- rWishart(n = 1, df = df, Sigma = Sigma)
    base_cov_matrix <- wishart_output[,,1]
    scaled_cov_matrix <- base_cov_matrix * psi / (2**i)
    cov_matrix <- gips::project_matrix(scaled_cov_matrix, perm)
    
    class_data <- MASS::mvrnorm(n = n_per_class, mu = class_means[, k], Sigma = cov_matrix)
    class_df <- as.data.frame(class_data)
    class_df$Y <- as.factor(k)
    data_list[[k]] <- class_df
  }
  final_data <- do.call(rbind, data_list)
  
  # ** Train and evaluate classification models **
  set.seed(42) # Use a fixed seed for the split for consistency
  train_ratio <- 0.6
  train_indices <- sample(1:nrow(final_data), size = floor(train_ratio * nrow(final_data)))
  train_data <- final_data[train_indices, ]
  test_data <- final_data[-train_indices, ]
  
  # --- QDA Model ---
  qda_model <- MASS::qda(Y ~ ., data = train_data)
  qda_train_predictions <- predict(qda_model, train_data)$class
  qda_train_accuracy <- mean(qda_train_predictions == train_data$Y)
  qda_test_predictions <- predict(qda_model, test_data)$class
  qda_test_accuracy <- mean(qda_test_predictions == test_data$Y)
  
  cat(sprintf("QDA Accuracy -> Train: %.2f | Test: %.2f\n", qda_train_accuracy, qda_test_accuracy))
  
  # --- Check if the targets are met ---
  is_train_acc_met <- qda_train_accuracy > target_train_accuracy
  is_test_acc_met  <- qda_test_accuracy > target_test_accuracy
  
  if (is_train_acc_met && is_test_acc_met) {
    cat("\nSUCCESS: Found a suitable divisor!\n")
    found_divisor <- TRUE
    
    # Store final results
    final_divisor <- 2**i
    final_data_to_save <- final_data
    final_qda_train_accuracy <- qda_train_accuracy
    final_qda_test_accuracy <- qda_test_accuracy
    
    # Train and evaluate LDA for the final report
    lda_model <- MASS::lda(Y ~ ., data = train_data)
    lda_train_predictions <- predict(lda_model, train_data)$class
    final_lda_train_accuracy <- mean(lda_train_predictions == train_data$Y)
    lda_test_predictions <- predict(lda_model, test_data)$class
    final_lda_test_accuracy <- mean(lda_test_predictions == test_data$Y)
    
    break # Exit the loop
  }
  
  cat("---------------------------\n")
}

#-----------------------------------------------------------------------
#                        ** FINAL RESULTS & ACTIONS **
#-----------------------------------------------------------------------
cat("\n\n--- FINAL SIMULATION RESULTS ---\n")

if (found_divisor) {
  cat(sprintf("Optimal 'psi' found: %f\n\n", psi / final_divisor))
  
  # Save the data frame to a CSV file
  cat(sprintf("Saving the generated dataset to '%s'...\n", output_filename))
  write.csv(final_data_to_save, output_filename, row.names = FALSE)
  cat("File saved successfully.\n\n")
  
  cat("LDA Model:\n")
  cat("  - Training Accuracy:", round(final_lda_train_accuracy * 100, 2), "%\n")
  cat("  - Test Accuracy:    ", round(final_lda_test_accuracy * 100, 2), "%\n")
  cat("\n")
  cat("QDA Model:\n")
  cat("  - Training Accuracy:", round(final_qda_train_accuracy * 100, 2), "%\n")
  cat("  - Test Accuracy:    ", round(final_qda_test_accuracy * 100, 2), "%\n")
  
} else {
  cat(sprintf("FAILED: Could not find a suitable divisor within %d iterations.\n", max_iterations))
  cat("Consider increasing `max_iterations` or adjusting the targets\n")
}

cat("--------------------------------\n")