plot_matrix_columns <- function(mat,
                                main = "Matrix Column Plot",
                                xlab = "X",
                                ylab = "Value") {
  # Check inputs
  if (!is.matrix(mat)) stop("mat must be a matrix")
  if (ncol(mat) < 2) stop("Matrix must have at least 2 columns")

  # Extract x and y
  x <- mat[,1]
  Y <- mat[,-1, drop = FALSE]

  # Colors for each Y column
  cols <- grDevices::rainbow(ncol(Y))

  # First plot
  plot(x, Y[,1], type = "l", col = cols[1], lwd = 2,
       main = main, xlab = colnames(mat)[1], ylab = ylab)

  # Additional columns
  if (ncol(Y) > 1) {
    for (i in 2:ncol(Y)) {
      lines(x, Y[,i], col = cols[i], lwd = 2)
    }
  }

  # Legend based on column names (or indices if no names)
  legend_labels <- colnames(Y)
  if (is.null(legend_labels)) legend_labels <- paste0("Series ", 1:ncol(Y))

  legend("topright", legend = legend_labels,
         col = cols, lwd = 2, cex = 0.8, bg = "white")
}



colnames(results)[1] <- "number of observations"

plot_matrix_columns(results, ylab = "proportion of correctly classified", main = "methods comparison")