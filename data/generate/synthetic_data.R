# This script is designed to generate synthetic data for a classification problem.
# The data generation process consists of the following steps:

# 1. Sampling the S parameter for the Wishart distribution:
#    - A matrix S is generated, which must be a positive-definite matrix.
#    - This matrix will serve as the scale matrix parameter for the Wishart distribution.

# 2. Sampling covariance matrices:
#    - Covariance matrices are drawn from the Wishart distribution W_p(S, df).

# 3. Permutations of covariance matrices:
#    - Before being used, the sampled covariance matrices are projected onto specific permutations.
#    - The projection method depends on one of the three scenarios described below.

# 4. Sampling means:
#    - The mean vectors for each class are sampled as random points from a p-dimensional unit cube, R^p.

# 5. Modifying class separability (psi parameter):
#    - An additional parameter 'psi' is introduced, which multiplies the covariance matrix.
#    - The purpose of this parameter is to control the degree of linear separability between classes;
#      a larger 'psi' results in greater overlap between the classes.

# 6. Data generation:
#    - The final data for each class is sampled from a multivariate normal distribution
#      using the means and covariance matrices prepared in the previous steps.


# Three scenarios for data generation are considered, differing in the structure
# of parameters for each class:

# Scenario 1: gipsLDA model
# - For each class, data is generated using:
#   - THE SAME permutations.
#   - THE SAME covariance matrices.
#   - DIFFERENT means.
# This is analogous to the assumption in classical Linear Discriminant Analysis (LDA),
# which assumes a common covariance matrix for all classes.

# Scenario 2: gipsMultQDA model
# - For each class, data is generated using:
#   - THE SAME permutations.
#   - DIFFERENT covariance matrices.
#   - DIFFERENT means.
# In this scenario, each class has its own covariance matrix, but the
# dependency structure (defined by the permutations) remains the same across classes.

# Scenario 3: gipsQDA model
# - For each class, data is generated using:
#   - DIFFERENT permutations.
#   - DIFFERENT covariance matrices.
#   - DIFFERENT means.
# This is the most flexible model, where each class has its unique mean,
# covariance matrix, and dependency structure (permutation), which is analogous
# to the assumptions in Quadratic Discriminant Analysis (QDA).

# Scenario 4: Classic LDA model
# - For each class, data is generated using:
#   - THE SAME covariance matrix.
#   - DIFFERENT means.
# This serves as a baseline for comparison, representing the classic LDA assumption

# Scenario 5: Classic QDA model
# - For each class, data is generated using:
#   - DIFFERENT covariance matrices.
#   - DIFFERENT means.
# This serves as a baseline for comparison, representing the classic QDA assumption.

# Set the main seed for the entire simulation process
set.seed(44)

# --- Execute Scenarios ---

cat("=========================================================\n")
cat("=                 STARTING SIMULATION                   =\n")
cat("=========================================================\n\n")

# --- Scenario 1: gipsLDA ---
# Same permutation, same covariance matrix, different means.
source("generate/scenario_1_gipsLDA.R")

# --- Scenario 2: gipsMultQDA ---
# Same permutation, DIFFERENT covariance matrices, different means.
source("generate/scenario_2_gipsMultQDA.R")

# --- Scenario 3: gipsQDA ---
# DIFFERENT permutations, DIFFERENT covariance matrices, different means.
source("generate/scenario_3_gipsQDA.R")

# --- Scenario 4: Classic LDA ---
# Same covariance matrix (no projection), different means.
source("generate/scenario_4_classicLDA.R")

# --- Scenario 5: Classic QDA ---
# DIFFERENT covariance matrices (no projection), different means.
source("generate/scenario_5_classicQDA.R")


cat("=========================================================\n")
cat("=                ALL SCENARIOS COMPLETE                 =\n")
cat("=========================================================\n\n")
