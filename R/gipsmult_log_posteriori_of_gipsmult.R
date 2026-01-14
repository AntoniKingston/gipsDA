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
