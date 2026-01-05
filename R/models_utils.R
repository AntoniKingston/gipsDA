#' Project Empirical Covariance Matrices onto Permutation Symmetry
#'
#' @description
#' A core internal utility that optimizes the permutation symmetry for a list of
#' empirical covariance matrices and projects them onto the discovered structure.
#' It handles both the Maximum A Posteriori (MAP) estimation and Bayesian Model
#' Averaging (weighted average) approaches.
#'
#' @details
#' This function serves as a bridge between the user-facing model functions
#' (like \code{\link{gipsLDA}} or \code{\link{gipsMultQDA}}) and the \code{gipsmult}
#' backend.
#'
#' It operates in two modes controlled by the \code{MAP} parameter:
#' \itemize{
#'   \item \strong{MAP = TRUE}: Finds the single most probable permutation group
#'   (argmax of the posterior) using \code{\link{find_MAP}}. All input matrices are
#'   projected onto this single permutation.
#'   \item \strong{MAP = FALSE}: Estimates the posterior distribution of permutations.
#'   The resulting covariance matrices are calculated as a weighted average of
#'   matrices projected onto all visited permutations, weighted by their posterior
#'   probabilities.
#' }
#'
#' To improve computational efficiency when \code{MAP = FALSE}, permutations with
#' a posterior probability below \code{tol} are discarded. If no permutation exceeds
#' this threshold, the function falls back to the single best permutation (MAP)
#' and issues a warning.
#'
#' @param emp_covs A list of numeric matrices. Each matrix represents the empirical
#'   covariance matrix of a class or group. They must be square and symmetric.
#' @param ns_obs A numeric vector of integers indicating the number of observations
#'   used to estimate each matrix in \code{emp_covs}.
#' @param MAP Logical. If \code{TRUE}, projects onto the single best permutation.
#'   If \code{FALSE}, uses a weighted average over the posterior distribution.
#'   Defaults to \code{TRUE}.
#' @param optimizer Character string specifying the optimization algorithm:
#'   \code{"BF"} (Brute Force) or \code{"MH"} (Metropolis-Hastings).
#' @param max_iter Integer. The number of iterations for the Metropolis-Hastings
#'   algorithm. Ignored if \code{optimizer = "BF"}.
#' @param tol Numeric. The probability threshold for the weighted average estimator.
#'   Permutations with posterior probability less than \code{tol} are ignored
#'   during projection. Defaults to \code{1e-3}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{covs}}{A list of projected covariance matrices (same length as \code{emp_covs}).}
#'   \item{\code{opt_info}}{Depending on \code{MAP}:
#'     \itemize{
#'       \item If \code{MAP = TRUE}: The optimal \code{gips_perm} object found.
#'       \item If \code{MAP = FALSE}: A named numeric vector of posterior probabilities used for weighting.
#'     }
#'   }
#' }
#'
#' @seealso
#' \code{\link{gipsLDA}}, \code{\link{gipsMultQDA}}, \code{\link{find_MAP}}
#'
#' @keywords internal
project_covs <- function(emp_covs, ns_obs, MAP = TRUE, optimizer, max_iter, tol = 1e-3) {
    gg <- gipsmult(emp_covs, ns_obs, was_mean_estimated = TRUE)
    if (MAP) {
        gg <- find_MAP(gg, optimizer = optimizer, max_iter = max_iter)
        perm <- gg[[1]]
        return(list(covs = lapply(emp_covs, function(x) gips::project_matrix(x, perm)), opt_info = perm))
    }
    gg <- find_MAP(gg, optimizer = optimizer, max_iter = max_iter, return_probabilities = TRUE, save_all_perms = TRUE)
    probs <- get_probabilities_from_gipsmult(gg)
    if (all(probs <= tol)) {
        warning("There are no perms with estimated probability above threshold, projecting onto MAP")
        probs <- probs[1]
    }
    else {
        probs <- probs[probs > tol]
    }


    return(list(covs = lapply(emp_covs, function(x) project_matrix_multiperm(x, probs)), opt_info = probs))
}


#' Weighted Average Projection of Covariance Matrix
#'
#' @description
#' Calculates a weighted average of the empirical covariance matrix projected onto
#' multiple different permutation symmetries. This implements the Bayesian Model
#' Averaging (BMA) estimator.
#'
#' @details
#' Instead of selecting a single "best" permutation, this function iterates through
#' a set of permutations provided in the names of the \code{probs} vector.
#' For each permutation, it projects the \code{emp_cov} onto that symmetry using
#' \code{\link[gips]{project_matrix}} and scales it by the corresponding probability.
#'
#' The final result is normalized by dividing by the sum of the provided probabilities.
#' This ensures the output is a valid weighted average even if the input \code{probs}
#' have been filtered (e.g., by a tolerance threshold in \code{\link{project_covs}})
#' and do not sum exactly to 1.
#'
#' @param emp_cov A numeric matrix. The empirical covariance matrix to be projected.
#'   Must be square and symmetric.
#' @param probs A named numeric vector.
#'   \itemize{
#'     \item \strong{Values}: The posterior probabilities (weights) associated with each permutation.
#'     \item \strong{Names}: Character strings representing the permutations (e.g., \code{"(1,2)(3,4)"})
#'           that are passed to the \code{gips} projection function.
#'   }
#'
#' @return A numeric matrix of the same dimensions as \code{emp_cov}, representing
#'   the regularized covariance estimate.
#'
#' @seealso \code{\link{project_covs}}
#'
#' @keywords internal
project_matrix_multiperm <- function(emp_cov, probs) {
    perms <- names(probs)
    projected_matrix <- matrix(0, nrow = dim(emp_cov), dim(emp_cov))
    for (i in 1:length(probs)) {
        projected_matrix <- projected_matrix + probs[[i]] * gips::project_matrix(emp_cov, perms[i])
    }
    return(projected_matrix / sum(probs))
}


#' Prepare R Objects for Custom JSON Serialization
#'
#' @description
#' Recursively traverses an R object and converts complex types (such as formulas,
#' language calls, matrices, and custom \code{gips_perm} objects) into a simplified
#' list structure containing type metadata. This prepares the object for safe
#' conversion to JSON using \code{jsonlite}.
#'
#' @details
#' Standard JSON serializers often lose R-specific type information (e.g., distinguishing
#' a formula from a string, or a matrix from a flat array). This function acts as a
#' pre-processor that wraps specific types in a list with a \code{"__type"} tag:
#' \itemize{
#'   \item \strong{gips_perm}: Converted to its character representation with an
#'         added \code{size} field (calculated via \code{recursive_length}).
#'   \item \strong{formula/terms}: Converted to a string representation via \code{deparse}.
#'   \item \strong{language (calls)}: Converted to a string.
#'   \item \strong{matrix}: Converted to a list containing the flattened data vector
#'         and dimensions (\code{nrow}, \code{ncol}).
#' }
#' Lists are traversed recursively. Atomic vectors are returned as-is.
#'
#' @param x An arbitrary R object to be prepared for serialization. Typically
#'   a component of a trained \code{gipsDA} model object.
#'
#' @return A list structure (or atomic vector) ready to be passed to
#'   \code{jsonlite::toJSON}. Complex objects are replaced by lists containing
#'   reconstruction metadata.
#'
#' @seealso \code{\link{deserialize_from_json}}, \code{\link{gipsDA_to_json}}
#'
#' @keywords internal
serialize_for_json <- function(x) {
  # Handle gips_perm
  if (inherits(x, "gips_perm")) {
    return(list(
      `__type` = "gips_perm",
      value = as.character(x),
      size = recursive_length(x)
    ))
  }

  # Handle calls/language objects
  if (typeof(x) == "language") {
    return(list(
      `__type` = "call",
      value = paste(deparse(x), collapse = " ")
    ))
  }

  # Handle formulas
  if (inherits(x, "formula") || inherits(x, "terms")) {
    return(list(
      `__type` = "formula",
      value = paste(deparse(x), collapse = " ")
    ))
  }

  # Handle matrices explicitly
  if (is.matrix(x)) {
    return(list(
      `__type` = "matrix",
      data = as.vector(x),
      nrow = nrow(x),
      ncol = ncol(x)
    ))
  }

  # Recurse into lists
  if (is.list(x)) {
    return(lapply(x, serialize_for_json))
  }

  x
}


#' Reconstruct R Objects from Custom JSON Structure
#'
#' @description
#' Recursively traverses a list structure (typically imported from JSON) and
#' restores specific R objects based on the \code{"__type"} metadata tag. This is the
#' inverse operation of \code{\link{serialize_for_json}}.
#'
#' @details
#' This function inspects the input \code{x} for a specific list component named
#' \code{"__type"}. Based on its value, it reconstructs the original R object:
#' \itemize{
#'   \item \strong{"gips_perm"}: Recreates the permutation object using
#'         \code{\link[gips]{gips_perm}} with the stored character string and size.
#'   \item \strong{"formula"}: Converts the stored string back to a formula object
#'         using \code{\link{as.formula}}.
#'   \item \strong{"call"}: Converts the stored string back to a language call
#'         using \code{\link{str2lang}}.
#'   \item \strong{"matrix"}: Reconstructs the matrix from the flattened data vector
#'         and stored dimensions (\code{nrow}, \code{ncol}).
#' }
#' If the input is a standard list without the type tag, the function recurses
#' into its elements. Atomic vectors are returned unchanged.
#'
#' @param x A list (or nested list structure) resulting from parsing a JSON file,
#'   potentially containing the special \code{__type} metadata fields.
#'
#' @return An R object with restored classes and structures (e.g., a fully functional
#'   \code{gipsLDA} model object).
#'
#' @seealso \code{\link{serialize_for_json}}, \code{\link{gipsDA_from_json}}
#'
#' @keywords internal
deserialize_from_json <- function(x) {
  if (is.list(x) && !is.null(x$`__type`)) {
    type <- x$`__type`

    if (type == "gips_perm") {
      return(gips::gips_perm(x$value, x$size))
    }

    if (type == "call") {
      return(str2lang(x$value))
    }

    if (type == "formula") {
      return(as.formula(x$value))
    }

    if (type == "matrix") {
      return(matrix(x$data, nrow = x$nrow, ncol = x$ncol))
    }
  }

  if (is.list(x)) {
    return(lapply(x, deserialize_from_json))
  }

  x
}


#' Save gipsDA Model to a JSON File
#'
#' @description
#' Serializes a trained \code{gipsDA} model object (e.g., \code{gipsLDA}, \code{gipsQDA})
#' and writes it to a file in JSON format. This function ensures that complex R structures
#' like formulas, call objects, and custom \code{gips_perm} classes are correctly preserved.
#'
#' @details
#' This function acts as a wrapper around \code{\link[jsonlite]{write_json}}. Before writing,
#' it processes the object using an internal serialization routine to handle types that
#' are not natively supported by JSON (or lose information upon standard conversion).
#'
#' Specifically, it converts:
#' \itemize{
#'   \item \code{gips_perm} objects to their character representation and size.
#'   \item Matrices to a list containing data and dimensions.
#'   \item Formulas and calls to string representations.
#' }
#' The resulting JSON file is formatted to be human-readable (\code{pretty = TRUE}).
#'
#' @param obj A trained model object (of class \code{gipsLDA}, \code{gipsQDA}, or \code{gipsMultQDA}).
#' @param file A character string naming the file to write to.
#'
#' @return \code{NULL}, invisibly. The function is called for its side effect of writing a file.
#'
#' @seealso \code{\link{gipsDA_from_json}} for loading the model back into R.
#'
#' @examples
#' \dontrun{
#' # Train a model
#' data(iris)
#' model <- gipsLDA(Species ~ ., data = iris)
#'
#' # Save to a temporary file
#' tmp_file <- tempfile(fileext = ".json")
#' gipsDA_to_json(model, tmp_file)
#'
#' # Check if file exists
#' file.exists(tmp_file)
#'
#' # Clean up
#' unlink(tmp_file)
#' }
#'
#' @export
gipsDA_to_json <- function(obj, file) {
  jsonlite::write_json(
    serialize_for_json(obj),
    file,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}


#' Load a gipsDA Model from a JSON File
#'
#' @description
#' Reads a JSON file containing a serialized \code{gipsDA} model and reconstructs
#' the original R object. This function restores complex structures (formulas,
#' matrices, permutation objects) and reassigns the specified S3 class.
#'
#' @details
#' This function performs the inverse operation of \code{\link{gipsDA_to_json}}.
#' It uses \code{\link[jsonlite]{read_json}} to load the raw list structure and then
#' applies an internal deserializer to convert specific metadata-tagged lists back
#' into their original R types (e.g., converting strings back to formulas or
#' reconstructing \code{gips_perm} objects).
#'
#' Since the JSON format does not natively store R S3 class attributes, the
#' specific class of the model must be provided manually via the \code{classname}
#' argument to ensure method dispatch (e.g., for \code{predict}) works correctly.
#'
#' @param file A character string naming the JSON file to read.
#' @param classname A character vector specifying the S3 class to assign to the
#'   reconstructed object. Typically, this should match the class of the original
#'   model (e.g., \code{"gipsLDA"}, \code{"gipsQDA"}, or \code{c("gipsLDA", "lda")}).
#'
#' @return A trained model object with the restored structure and class, ready
#'   for use with standard methods like \code{predict}.
#'
#' @seealso \code{\link{gipsDA_to_json}}
#'
#' @examples
#' \dontrun{
#' # 1. Train and save a model
#' data(iris)
#' original_model <- gipsLDA(Species ~ ., data = iris)
#' tmp_file <- tempfile(fileext = ".json")
#' gipsDA_to_json(original_model, tmp_file)
#'
#' # 2. Load the model back
#' # Note: We explicitly provide the class name "gipsLDA"
#' loaded_model <- gipsDA_from_json(tmp_file, classname = "gipsLDA")
#'
#' # 3. Verify it works
#' predict(loaded_model, newdata = iris[1:5, ])
#'
#' # Clean up
#' unlink(tmp_file)
#' }
#'
#' @export
gipsDA_from_json <- function(file, classname) {
  raw <- jsonlite::read_json(file)
  obj <- deserialize_from_json(raw)
  class(obj) <- classname
  obj
}


#' Recursively Calculate Total Number of Atomic Elements
#'
#' @description
#' Traverses a potentially nested list structure and calculates the sum of lengths
#' of all atomic elements found at the leaves of the structure.
#'
#' @details
#' This function checks if the input `x` is atomic. If so, it returns its length.
#' If `x` is a list, it recursively applies itself to every element of the list
#' and sums the results. For other types (e.g., NULL), it returns 0.
#'
#' @param x An R object. Typically a list (which can be nested) or an atomic vector.
#'
#' @return An integer representing the total count of atomic elements within the structure.
#'
#' @examples
#' # Case 1: Atomic vector
#' recursive_length(c(1, 2, 3))
#' # Expected output: 3
#'
#' # Case 2: Simple list
#' recursive_length(list(a = 1, b = c(2, 3)))
#' # Expected output: 3
#'
#' # Case 3: Nested list
#' nested_list <- list(
#'   group1 = c(1, 2),
#'   group2 = list(
#'     subA = 3,
#'     subB = c(4, 5)
#'   )
#' )
#' recursive_length(nested_list)
#' # Expected output: 5
#'
#' @export
recursive_length <- function(x) {
  if (is.atomic(x)) {
    return(length(x))
  }
  if (is.list(x)) {
    return(sum(vapply(x, recursive_length, integer(1))))
  }
  return(0)
}


#' Regularize Matrix to Ensure Minimum Eigenvalue
#'
#' @description
#' Checks the smallest eigenvalue of a matrix and, if necessary, regularizes it
#' by mixing it with the Identity matrix (shrinkage) to ensure the smallest
#' eigenvalue is at least equal to a specified \code{target}.
#'
#' @details
#' The function computes the eigenvalues of \code{A}. If the smallest eigenvalue
#' \eqn{\lambda} is already greater than or equal to \code{target}, the matrix
#' is returned unchanged.
#'
#' If \eqn{\lambda < target}, the function calculates a scaling factor \eqn{s}
#' such that the modified matrix:
#' \deqn{A_{new} = \frac{A + sI}{1 + s}}
#' has a minimum eigenvalue exactly equal to \code{target}. This is effectively
#' a linear shrinkage towards the Identity matrix (which has all eigenvalues equal to 1).
#' This technique is used to ensure positive definiteness or numerical stability
#' when estimating covariance matrices from high-dimensional or singular data.
#'
#' @param A A numeric matrix (typically a covariance matrix). It is expected to be symmetric.
#' @param target A numeric scalar representing the minimum desired eigenvalue.
#'   Defaults to \code{0.05}. It must be strictly less than 1 to avoid division by zero.
#'
#' @return A numeric matrix of the same dimensions as \code{A}. If regularization
#'   was applied, the matrix is strictly positive definite with the smallest
#'   eigenvalue equal to \code{target}.
#'
#' @keywords internal
desingularize <- function(A, target = 0.05) {
  symmetric <- all.equal(A, t(A))
  eigvals <- eigen(A, symmetric = symmetric, only.values = TRUE)$values
  idx <- which.min(abs(eigvals))
  lambda <- eigvals[idx]

  if (abs(lambda) >= target) {
    return(A)
  }

  s <- (target - lambda) / (1 - target)

  if (1 + s <= 0) {
    stop("Invalid scaling: 1 + s <= 0")
  }

  return((A + diag(s, nrow(A))) / (1 + s))
}
