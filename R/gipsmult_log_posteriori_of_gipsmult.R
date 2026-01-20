#' A log of a posteriori that the covariance matrix is invariant under permutation
#'
#' More precisely, it is the logarithm of an unnormalized
#' posterior probability. It is the goal function for
#' optimization algorithms in the `find_MAP()` function.
#' The `perm_proposal` that maximizes this function is
#' the Maximum A Posteriori (MAP) Estimator.
#'
#' It is calculated using
#' [formulas (33) and (27) from references](https://arxiv.org/abs/2004.03503).
#'
#' If `Inf` or `NaN` is reached, it produces a warning.
#'
#' @export
#'
#' @param g An object of a `gipsmult` class.
#'
#' @references Piotr Graczyk, Hideyuki Ishi, Bartosz Kołodziejek, Hélène Massam.
#' "Model selection in the space of Gaussian models invariant by symmetry."
#' The Annals of Statistics, 50(3) 1747-1774 June 2022.
#' [arXiv link](https://arxiv.org/abs/2004.03503);
#' \doi{10.1214/22-AOS2174}
#'
#' @seealso
#' * [find_MAP()] - The function that optimizes
#'     the `log_posteriori_of_gips` function.
#' * [compare_posteriories_of_perms()] - Uses `log_posteriori_of_gips()`
#'     to compare a posteriori of two permutations.
#'
#' @returns Returns a value of
#'     the logarithm of an unnormalized A Posteriori.
#'
#' @examples
#' # In the space with p = 2, there is only 2 permutations:
#' perm1 <- permutations::as.cycle("(1)(2)")
#' perm2 <- permutations::as.cycle("(1,2)")
#' S1 <- matrix(c(1, 0.5, 0.5, 2), nrow = 2, byrow = TRUE)
#' S2 <- matrix(c(2, 1, 3, 7), nrow = 2, byrow = TRUE)
#' g1 <- gipsmult(list(S1,S2), c(100,100), perm = perm1)
#' g2 <- gipsmult(list(S1,S2), c(100,100), perm = perm2)
#' log_posteriori_of_gipsmult(g1) # -354.4394, this is the MAP Estimator
#' log_posteriori_of_gipsmult(g2) # -380.0079
#'
#' exp(log_posteriori_of_gipsmult(g1) - log_posteriori_of_gipsmult(g2)) # 127131902082
#' # g1 is 127131902082 times more likely than g2.
#' # This is the expected outcome because S1[1,1] and S2[1,1]
#' # differ significantly from S1[2,2] and S2[2,2] respectively.
#'
#' # ========================================================================
#'
#' S3 <- matrix(c(1, 0.5, 0.5, 1.1), nrow = 2, byrow = TRUE)
#' S4 <- matrix(c(2, 1, 3, 2.137), nrow = 2, byrow = TRUE)
#' g1 <- gipsmult(list(S3,S4), c(100,100), perm = perm1)
#' g2 <- gipsmult(list(S3,S4), c(100,100), perm = perm2)
#' log_posteriori_of_gipsmult(g1) # -148.6485
#' log_posteriori_of_gipsmult(g2) # -145.3019, this is the MAP Estimator
#'
#' exp(log_posteriori_of_gipsmult(g2) - log_posteriori_of_gipsmult(g1)) # 28.406
#' # g2 is 28.406 times more likely than g1.
#' # This is the expected outcome because S1[1,1] and S2[1,1]
#' # are very close to S1[2,2] and S2[2,2] respectively.
log_posteriori_of_gipsmult <- function(g) {

  numbers_of_observations <- attr(g, "numbers_of_observations")
  was_mean_estimated <- attr(g, "was_mean_estimated")

  if (was_mean_estimated) {
    edited_numbers_of_observations <- numbers_of_observations - 1
  } else {
    edited_numbers_of_observations <- numbers_of_observations
  }

  log_posteriori_of_perm(
    perm_proposal = g[[1]], Ss = attr(g, "Ss"),
    numbers_of_observations = edited_numbers_of_observations,
    delta = attr(g, "delta"), D_matrices = attr(g, "D_matrices")
  )
}

#' We recommend to use the `log_posteriori_of_gips()` function.
#'
#' If You really want to use `log_posteriori_of_perm()`, remember
#'  to edit `number_of_observations` if the mean was estimated!
#'
#' @noRd
log_posteriori_of_perm <- function(perm_proposal, Ss, numbers_of_observations,
                                         delta, D_matrices){
  log_values <- mapply(
    function(S, n_obs, D_matrix) {
      gips:::log_posteriori_of_perm(perm_proposal, S, n_obs, delta, D_matrix)
    },
    Ss,
    numbers_of_observations,
    D_matrices
  )

  sum(log_values)

}