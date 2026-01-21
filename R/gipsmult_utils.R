#' List of matrices check
#'
#' Cheks whether provided argument is a list of matrices or not
#'
#' @noRd
list_of_matrices_check <- function(x) {
  if (!is.list(x)) {
    return(FALSE)
  }
  for (y in x) {
    if (!is.matrix(y)) {
      return(FALSE)
    }
  }
  return(TRUE)
}


#' SDN compatibility check
#'
#' Checks wheter 2 first arguments have the same length as the 3rd and all of their elements the same dimension or not.
#'
#' @noRd
SDN_compatibility_check <- function(Ss, D_matrices, numbers_of_observations) {
  n <- length(numbers_of_observations)
  if (length(Ss) != n | length(D_matrices) != n) {
    return(FALSE)
  }
  if (!all(dim(D_matrices[[1]]) == dim(Ss[[1]]))) {
    return(FALSE)
  }
  return(TRUE)
}

#' Noo check
#'
#' Checks wheter all elements of provided argument are whole numbers or not.
#'
#' @noRd
noo_check <- function(numbers_of_observations) {
  return(all(is.wholenumber(numbers_of_observations)))
}
#' Plot single gg
#'
#' Plots single matrix assuming one has ggplot2 installed.
#'
#' @noRd
plot_single_gg <- function(my_projected_matrix, perm) {
  p <- ncol(my_projected_matrix)

  if (is.null(colnames(my_projected_matrix))) {
    colnames(my_projected_matrix) <- paste0(seq(1, p))
  }
  if (is.null(rownames(my_projected_matrix))) {
    rownames(my_projected_matrix) <- paste0(seq(1, p))
  }

  my_rownames <- rownames(my_projected_matrix)
  my_colnames <- colnames(my_projected_matrix)
  rownames(my_projected_matrix) <- as.character(1:p)
  colnames(my_projected_matrix) <- as.character(1:p)

  # With this line, the R CMD check's "no visible binding for global variable" warning will not occur:
  col_id <- covariance <- row_id <- NULL

  # Life would be easier with pipes (%>%)
  my_transformed_matrix <- tibble::rownames_to_column(
    as.data.frame(my_projected_matrix),
    "row_id"
  )
  my_transformed_matrix <- tidyr::pivot_longer(my_transformed_matrix,
    -c(row_id),
    names_to = "col_id",
    values_to = "covariance"
  )
  my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
    col_id = as.numeric(col_id)
  )
  my_transformed_matrix <- dplyr::mutate(my_transformed_matrix,
    row_id = as.numeric(row_id)
  )
  g_plot <- ggplot2::ggplot(
    my_transformed_matrix,
    ggplot2::aes(x = col_id, y = row_id, fill = covariance)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(na.value = "white") +
    ggplot2::scale_x_continuous(breaks = 1:p, labels = my_rownames) +
    ggplot2::scale_y_reverse(breaks = 1:p, labels = my_colnames) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("Estimated covariance matrix\nprojected on permutation ", perm),
      x = "", y = ""
    )

  return(g_plot)
}

#' Plot single stats
#'
#' Plots single matrix assuming one does not have  ggplot2 installed.
#'
#' @noRd
plot_single_stats <- function(my_projected_matrix, color, ...) {
  if (is.null(color)) { # Setting col = NA or col = NULL turns off the whole plot.
    stats::heatmap(my_projected_matrix,
      symm = TRUE,
      Rowv = NA, Colv = NA, ...
    )
  } else {
    stats::heatmap(my_projected_matrix,
      symm = TRUE,
      Rowv = NA, Colv = NA, col = color, ...
    )
  }
}


#' Get diagonalized matrices for heatmap
#'
#' Applies get_diagonalized_matrix_for_heatmap() to a list od matrices of a gipsmult object.
#'
#' @noRd
get_diagonalized_matrices_for_heatmap <- function(x) {
  Ss <- attr(x, "Ss")
  lapply(Ss, get_diagonalized_matrix_for_heatmap)
}
