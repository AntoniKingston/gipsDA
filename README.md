# gipsDA

[![Build Status](https://img.shields.io/travis/com/your_username/gipsDA.svg?style=flat-square)](https://travis-ci.com/your_username/gipsDA)
[![PyPI version](https://img.shields.io/pypi/v/gipsDA.svg?style=flat-square)](https://pypi.org/project/gipsDA/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

An R/Python package for discriminant analysis classification using covariance matrices with permutation symmetries.

## About The Project

**gipsDA** extends classical Linear Discriminant Analysis (LDA) and Quadratic Discriminant Analysis (QDA) by incorporating permutation group structures into the estimation of covariance matrices. By leveraging the methodology of the `gips` library, this package aims to improve classification performance in scenarios where features (variables) exhibit underlying symmetries.

The core idea is to find and impose a permutation symmetry on the covariance matrix, which acts as a form of regularization and can lead to more stable and interpretable models, especially in high-dimensional settings.

### Key Features

*   Implementation of four novel `gips`-based discriminant analysis classifiers.
*   A flexible, user-friendly model API consistent with established machine learning libraries.
*   A specialized `gipsmult` module for modeling class-specific covariances under a shared symmetry.
*   Includes scripts for generating complex synthetic datasets to rigorously test model performance.

## Project Structure

The repository is organized into the following key directories:

-   #### `data/`
    This directory contains all scripts and documentation related to the datasets used for model evaluation.
    -   **Synthetic Data:** Includes R scripts for generating the five synthetic data scenarios described in our study. The generation process is designed to create datasets where classes have specific, controlled covariance structures (e.g., shared vs. unique permutations, shared vs. unique covariance matrices).
    -   **Real-World Data:** This sub-directory provides information and links to the real-world datasets used for benchmarking. These datasets cover domains such as medical diagnostics, finance, and industrial quality control, offering a diverse range of classification challenges.

-   #### `tests/`
    This directory houses all testing scripts to ensure the reliability and correctness of the package.
    -   **Unit Tests:** Automated tests that verify the functionality of individual functions and methods within the modules (e.g., matrix calculations, object creation, prediction outputs).
    -   **Manual Tests:** Higher-level scripts designed to be run manually to validate entire workflows, such as replicating the results from our simulation studies or ensuring that the models train and predict correctly on specific datasets.

-   #### `gipsmult/`
    This directory contains the source code for the `gipsmult` module. This module extends the core functionality of the `gips` library. Its primary purpose is to enable the estimation of different covariance matrices for each class under the constraint that they all share the same underlying permutation symmetry. This is designed for scenarios where the classes may have different variances but are assumed to possess a common dependency structure.

-   #### `gipsDA/`
    This is the main, user-facing module of the package. It provides the `gipsDA` class, which serves as the primary interface for creating, training, and using the classification models. The behavior of a `gipsDA` object is controlled by a hyperparameter that specifies which of the four implemented algorithms to use.

## Available Models in `gipsDA`

The `gipsDA` object can be configured to run one of four different classification algorithms, each with different assumptions about the data structure.

1.  **`gipsLDA_weighted_average`**
    In this approach, a separate covariance matrix is first estimated for each class using the standard `gips` library. A final, single covariance matrix is then computed as a **weighted average** of these individual matrices. This pooled matrix is subsequently used for classification within the LDA framework.

2.  **`gipsLDA_classic`**
    This model follows a more traditional approach by estimating a **single, pooled covariance matrix** directly from all classes at once using the `gips` library. This is analogous to the classic LDA assumption but with the added layer of permutation symmetry discovery.

3.  **`gipsMultQDA`**
    This model leverages the `gipsmult` module. The process involves first identifying the single most probable **permutation structure that is common across all classes**. Subsequently, a separate covariance matrix is estimated for each class, with each matrix being projected onto this shared permutation. This allows for class-specific covariance while maintaining a unified dependency model.

4.  **`gipsQDA`**
    This represents the most flexible model. The `gips` library is applied **independently to each class**. Consequently, each class can have its own uniquely estimated permutation structure and its own distinct covariance matrix. This is analogous to the classic QDA framework but with individualized symmetry discovery for each class.