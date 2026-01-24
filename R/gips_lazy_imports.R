## This file contains 1:1 pasted functions and expressions from gips library that are not exported, but needed
# To be changed when gipsmult gets merged into gips

OEIS_A000142 <- c(1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000, 51090942171709440000, 1124000727777607680000)
OEIS_A051625 <- c(1, 2, 5, 17, 67, 362, 2039, 14170, 109694, 976412, 8921002, 101134244, 1104940280, 13914013024, 191754490412, 2824047042632, 41304021782824, 708492417746000, 11629404776897384, 222093818836736752, 4351196253952132832, 88481681599705382144)

#' Calculate exact posterior probabilities
#'
#' We use the "second approach" from the paper.
#'
#' @param perms An output of `permutations::allperms()`.
#' @param log_posteriories A vector of all values of log posteriories of all `perms`.
#' @param show_progress_bar A boolean. Indicate whether or not to show two progress bars.
#'
#' @returns A named numeric vector. Names: character representations of permutations.
#' Elements: estimated posterior probabilities of permutations.
#' @noRd
calculate_probabilities <- function(perms, log_posteriories, show_progress_bar = FALSE) {
  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = length(perms), initial = 1)
  }

  for (i in 1:19) {
    perms_size <- i
    if (OEIS_A051625[i] == length(perms) || OEIS_A000142[i] == length(perms)) {
      break
    }
  }

  if (perms_size == 19) {
    rlang::abort("There is sth wrong with this sequence of permutations!",
      "i" = "The length of the permutation vector has to be an element of OEIS sequence A051625 or A000142",
      "x" = paste0("You have the length of permutation vector = ", length(perms))
    )
  }

  group_representatives <- character(0)
  for (i in 1:length(perms)) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }

    # We should use `g_perm <- gips_perm(perms[i], perms_size)`,
    # but this is faster, because it lacks safe checks:
    g_perm <- gips_perm_no_checks(perms[i], perms_size)

    group_representatives[i] <- as.character(get_group_representative(g_perm))
  }

  if (show_progress_bar) {
    close(progressBar)
  }


  # get rid of the repeated permutations:
  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = length(group_representatives), initial = 1)
  }
  groups_unrepeated_log_posteriories <- log_posteriories[1]
  names(groups_unrepeated_log_posteriories)[1] <- group_representatives[1]
  for (i in 2:length(group_representatives)) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }
    if (!(group_representatives[i] %in% names(groups_unrepeated_log_posteriories))) {
      groups_unrepeated_log_posteriories[length(groups_unrepeated_log_posteriories) + 1] <- log_posteriories[i]
      names(groups_unrepeated_log_posteriories)[length(groups_unrepeated_log_posteriories)] <- group_representatives[i]
    }
  }

  if (show_progress_bar) {
    close(progressBar)
  }

  sort(change_log_probabilities_unnorm_to_probabilities(groups_unrepeated_log_posteriories), decreasing = TRUE)
}

#' Compose permutation with transposition
#'
#' @param gips_perm Object of a `gips_perm` class.
#' @param transposition An integer vector of length 2. Transposition in a form of a
#' cycle.
#'
#' @returns An object of a `gips_perm` class. Composition of `gips_perm` parameter and `transposition`.
#'
#' @noRd
#' @examples
#' perm <- permutations::as.cycle("(1,2,3)(4,5)")
#' gperm <- gips_perm(perm, 6)
#' tr <- c(2, 3)
#' tr_perm <- permutations::as.cycle(tr)
#'
#' composed <- compose_with_transposition(gperm, tr)
#' composed2 <- perm * tr_perm
#'
#' # composed and composed 2 refer to the same permutation
compose_with_transposition <- function(gips_perm, transposition) {
  cycle_1_index <- which(sapply(gips_perm, function(cycle) {
    transposition[1] %in% cycle
  }))
  cycle_2_index <- which(sapply(gips_perm, function(cycle) {
    transposition[2] %in% cycle
  }))
  cycle_1 <- gips_perm[[cycle_1_index]]
  cycle_2 <- gips_perm[[cycle_2_index]]
  composed_gips_perm <- gips_perm[c(-cycle_1_index, -cycle_2_index)]
  if (cycle_1_index == cycle_2_index) {
    # We are breaking cycle into 2 cycles
    shifted_cycle <- shift_vector(cycle_1, which(cycle_1 == transposition[1]) - 1)
    new_cycle_1 <- shifted_cycle[1:(which(shifted_cycle == transposition[2]) - 1)]
    new_cycle_2 <- shifted_cycle[(which(shifted_cycle == transposition[2])):length(shifted_cycle)]
    composed_gips_perm <- add_cycle(composed_gips_perm, new_cycle_1)
    composed_gips_perm <- add_cycle(composed_gips_perm, new_cycle_2)
  } else {
    # We are merging 2 cycles
    ind <- which(cycle_1 == transposition[1])
    fragment_1 <- shift_vector(cycle_2, which(cycle_2 == transposition[2]) - 1)
    fragment_2 <- shift_vector(cycle_1, which(cycle_1 == transposition[1]) - 1)
    new_cycle <- c(fragment_1, fragment_2)
    composed_gips_perm <- add_cycle(composed_gips_perm, new_cycle)
  }
  gips::new_gips_perm(
    rearrange_cycles(composed_gips_perm),
    attr(gips_perm, "size")
  )
}

#' Convert the log difference to the appropriate string.
#' If bigger than 10 millions, use the scientific notation
#'
#' @param digits Number of digits after comma
#'
#' @examples
#' convert_log_diff_to_str(1009.5, 3) == "2.632e+438"
#' convert_log_diff_to_str(16.1, 3) == "9820670.922"
#' convert_log_diff_to_str(16.2, 3) == "1.085e+7"
#' convert_log_diff_to_str(-7.677, 3) == "4.634e-4"
#'
#' @noRd
convert_log_diff_to_str <- function(log_diff, digits) {
  if (is.infinite(log_diff)) {
    return(ifelse(log_diff > 0, "Inf", "-Inf"))
  }
  if (log_diff == 0) {
    return("1")
  }

  times_more_likely <- round(
    exp(log_diff),
    digits = digits
  )

  log10_diff <- log_diff * log10(exp(1))

  ifelse(0 < times_more_likely && times_more_likely < 10000000,
    as.character(times_more_likely),
    paste0(
      round(10^(log10_diff - floor(log10_diff)), digits = digits),
      "e",
      ifelse(log_diff > 0, "+", ""), # If log_diff < 0, then floor(log10_diff) will have "-" in front
      floor(log10_diff)
    )
  )
}

#' Estimate posterior probabilities from Metropolis-Hastings run
#'
#' We use the "second approach" from the paper.
#'
#' @param perms A list of `gips_perm` objects. Visited groups during the Metropolis-Hastings run.
#' @param show_progress_bar A boolean. Indicate whether or not to show the progress bar.
#'
#' @returns A named numeric vector. Names: character representations of permutations.
#' Elements: estimated posterior probabilities of permutations.
#' @noRd
estimate_probabilities <- function(perms, show_progress_bar = FALSE) {
  if (show_progress_bar) {
    progressBar <- utils::txtProgressBar(min = 0, max = length(perms), initial = 1)
  }
  group_representatives <- sapply(1:length(perms), function(i) {
    if (show_progress_bar) {
      utils::setTxtProgressBar(progressBar, i)
    }
    as.character(get_group_representative(perms[[i]]))
  })
  repr_counts <- table(group_representatives)
  if (show_progress_bar) {
    close(progressBar)
  }

  repr_weights <- sapply(names(repr_counts), function(p_str) {
    perm <- permutations::char2cycle(p_str)
    p_order <- permutations::permorder(perm)
    1 / numbers::eulersPhi(p_order)
  })
  unnormalized_probabilities <- repr_counts * repr_weights / length(perms)
  probabilities <- unnormalized_probabilities / sum(unnormalized_probabilities)
  probabilities <- as.numeric(probabilities)
  names(probabilities) <- names(unnormalized_probabilities)

  sort(probabilities, decreasing = TRUE)
}

#' Replace all non-block entries with NA
#'
#' Diagonalize matrix using found permutation and
#' replace all entries outside blocks (equal to 0) with NA.
#' This is done, because later these fields are plotted with background color.
#' It is more clear then.
#'
#' @param g `gips` object.
#' @noRd
get_diagonalized_matrix_for_heatmap <- function(g) {
  perm <- g[[1]]
  projected_matrix <- gips::project_matrix(attr(g, "S"), perm)
  diagonalising_matrix <- gips::prepare_orthogonal_matrix(perm)
  full_block_matrix <- t(diagonalising_matrix) %*% projected_matrix %*% diagonalising_matrix
  block_ends <- get_block_ends(gips::get_structure_constants(perm))
  block_starts <- c(1, block_ends[-length(block_ends)] + 1)
  block_matrix <- matrix(
    nrow = nrow(full_block_matrix),
    ncol = ncol(full_block_matrix)
  )
  for (i in 1:length(block_starts)) {
    slice <- block_starts[i]:block_ends[i]
    block_matrix[slice, slice] <- full_block_matrix[slice, slice, drop = FALSE]
  }
  block_matrix
}

#' Internal
#' @return (integer) n0
#' @noRd
get_n0_from_perm <- function(g_perm, was_mean_estimated) {
  structure_constants <- gips::get_structure_constants(g_perm)
  n0 <- max(structure_constants[["r"]] * structure_constants[["d"]] / structure_constants[["k"]])

  if (was_mean_estimated) { # correction for estimating the mean
    n0 <- n0 + 1
  }

  c(n0)
}

#' Is the numeric value representing the whole number
#'
#' This code is copied from [base::integer()] example.
#'
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (!is.numeric(x)) {
    return(rep(FALSE, length(x)))
  }
  abs(x - round(x)) < tol
}

#' We recommend to use the `log_posteriori_of_gips()` function.
#'
#' If You really want to use `log_posteriori_of_perm_gips()`, remember
#'  to edit `number_of_observations` if the mean was estimated!
#'
#' @noRd
log_posteriori_of_perm_gips <- function(perm_proposal, S, number_of_observations,
                                        delta, D_matrix) {
  U <- S * number_of_observations # in the paper there is U everywhere instead of S, so it is easier to use U matrix in the code
  perm_size <- dim(S)[1]

  if (!inherits(perm_proposal, "gips_perm")) {
    perm_proposal <- gips::gips_perm(perm_proposal, perm_size)
  }

  if (is.null(D_matrix)) {
    D_matrix <- diag(nrow = perm_size) # identity matrix
  }

  structure_constants <- gips::get_structure_constants(perm_proposal)

  # Ac_part
  Ac <- sum(structure_constants[["r"]] * structure_constants[["k"]] * log(structure_constants[["k"]])) # (20)
  Ac_part <- (-number_of_observations / 2 * Ac)

  # G_part and phi_part
  G_part <- G_function(structure_constants, delta + number_of_observations) -
    G_function(structure_constants, delta)

  # phi_part
  phi_part <- calculate_phi_part(
    perm_proposal, number_of_observations, U,
    delta, D_matrix, structure_constants
  )

  out <- Ac_part + G_part + phi_part

  if (is.infinite(out)) {
    rlang::warn("The infinite value of a posteriori was produced.")
  }
  if (is.nan(out)) {
    rlang::warn("The NaN value of a posteriori was produced.")
  }

  out
}

#' Uniformly random transposition of perm_size elements
#'
#' @param perm_size A size from which take transpositions.
#'
#' @noRd
runif_transposition <- function(perm_size) {
  sample(perm_size, 2, replace = FALSE)
}

#' The same as [gips_perm()], but lacks safety checks
#'
#' We highly advise against using this function.
#' Instead, use [gips_perm()]. Although, this one is slightly faster.
#' @noRd
gips_perm_no_checks <- function(x, size) {
  cycles <- unclass(x)[[1]]
  all_ints <- unlist(cycles)
  representatives <- permutations::get1(x)
  fixed_boolean <- permutations::fixed(x)
  if (length(fixed_boolean) < size) {
    fixed_boolean[(length(fixed_boolean) + 1):size] <- TRUE
  }
  fixed_elements <- which(fixed_boolean)
  subcycles <- c(cycles, as.list(fixed_elements))

  gips::new_gips_perm(rearrange_cycles(subcycles), size)
}

#' Get a representative of a cyclic permutation group
#'
#' Essentially a "nu" function from paper (beginning of section 4.1.1)
#' `get_representative(perm) == nu(< perm >)`, where "`nu`" is from the paper,
#' and "`< perm >`" is a cyclic group generated by the permutation.
#'
#' @param perm `gips_perm`
#'
#' @returns Object of a `gips_perm` class.
#' @noRd
get_group_representative <- function(perm) {
  size <- attr(perm, "size")
  if (size == 0) {
    return(perm)
  }
  perm <- permutations::as.cycle(perm)
  p_order <- permutations::permorder(perm)
  coprimes <- get_coprimes(p_order)
  all_perms <- lapply(coprimes, function(cp) {
    permutations::cycle_power(perm, cp)
  })
  if (length(all_perms) == 0) {
    return(gips::gips_perm(permutations::nullword, 0))
  }
  all_perms_as_vectors <- lapply(all_perms, function(p) {
    as.integer(permutations::cycle2word(p))
  })
  all_perms_as_strings <- sapply(all_perms_as_vectors, function(v) {
    paste(v, collapse = " ")
  })
  min_index <- stringi::stri_order(all_perms_as_strings)[1]
  gips::gips_perm(all_perms[[min_index]], size)
}

change_log_probabilities_unnorm_to_probabilities <- function(log_probabilities_unnorm) {
  log_probabilities_unnorm <- log_probabilities_unnorm - max(log_probabilities_unnorm)
  probabilities_unnorm <- exp(log_probabilities_unnorm)
  probabilities <- probabilities_unnorm / sum(probabilities_unnorm)

  probabilities
}

#' Shift vector
#'
#' Move k elements from the start of a vector to its end.
#'
#' @noRd
shift_vector <- function(v, k) {
  if (k == 0) {
    return(v)
  }
  c(v[-(1:k)], v[1:k])
}

#' Add a new cycle to permutation
#'
#' @param cycles A list of integer vectors. Each corresponds to cycles of a permutation.
#' @param new_cycle An integer vector. None of its elements are present in `cycles`.
#'
#' @noRd
add_cycle <- function(cycles, new_cycle) {
  # Assume, that cycles are sorted by their min element
  # new_cycle - not necessarily
  new_cycle <- rearrange_vector(new_cycle)
  min_representatives <- sapply(cycles, function(v) v[1])
  insert_index <- findInterval(new_cycle[1], min_representatives)
  append(cycles, list(new_cycle), after = insert_index)
}

#' Rearrange cycles
#'
#' `gips_perm` object stores permutations in cyclic form in following convention:
#' 1) cycles are ordered by their minimal element
#' 2) First element of a cycle is its minimal
#'
#' @param cycles A list of integer vectors.
#'
#' @examples
#' cycles <- list(c(2, 4, 3), c(5, 1))
#' rearranged <- rearrange_cycles(cycles)
#' # rearranged is list (c(1,5), c(2,4,3))
#' @noRd
rearrange_cycles <- function(cycles) {
  rearranged_cycles <- lapply(cycles, rearrange_vector)
  representatives <- sapply(rearranged_cycles, function(v) v[1])
  rearranged_cycles[order(representatives)]
}

#' Get a representative of a cyclic permutation group
#'
#' Essentially a "nu" function from paper (beginning of section 4.1.1)
#' `get_representative(perm) == nu(< perm >)`, where "`nu`" is from the paper,
#' and "`< perm >`" is a cyclic group generated by the permutation.
#'
#' @param perm `gips_perm`
#'
#' @returns Object of a `gips_perm` class.
#' @noRd
get_group_representative <- function(perm) {
  size <- attr(perm, "size")
  if (size == 0) {
    return(perm)
  }
  perm <- permutations::as.cycle(perm)
  p_order <- permutations::permorder(perm)
  coprimes <- get_coprimes(p_order)
  all_perms <- lapply(coprimes, function(cp) {
    permutations::cycle_power(perm, cp)
  })
  if (length(all_perms) == 0) {
    return(gips::gips_perm(permutations::nullword, 0))
  }
  all_perms_as_vectors <- lapply(all_perms, function(p) {
    as.integer(permutations::cycle2word(p))
  })
  all_perms_as_strings <- sapply(all_perms_as_vectors, function(v) {
    paste(v, collapse = " ")
  })
  min_index <- stringi::stri_order(all_perms_as_strings)[1]
  gips::gips_perm(all_perms[[min_index]], size)
}

get_block_ends <- function(structure_constants) {
  cumsum(structure_constants[["r"]] * structure_constants[["d"]])
}

#' G_function for `log_posteriori_of_gips()`
#'
#' @param delta Parameter of a method. The default is `3`.
#'     When `structure_constants` are from an id permutation, `delta <= 0` iff `G_function() = +Inf`.
#'     When `structure_constants` are from a permutation that is not id, `delta <= 1` iff `G_function() = +Inf`.
#' @param structure_constants Constants from `gips::get_structure_constants` function.
#'
#' @returns Sum of logarithms of elements of `calculate_gamma_omega` from i to L.
#' It is a log of a product part of the equation (27). For more information, see Issue #3 on `gips`' GitHub.
#'
#' @examples
#' perm_size <- 6
#' perm <- permutations::as.cycle(permutations::as.word(c(2, 3, 1, 5, 4, 6)))
#' my_gips_perm <- gips_perm(perm, perm_size)
#' structure_constants <- gips::get_structure_constants(my_gips_perm)
#' gips:::G_function(structure_constants, 3)
#'
#' @noRd
G_function <- function(structure_constants, delta = 3) {
  single_G_i <- sapply(1:structure_constants[["L"]], function(i) {
    lambda_i <- structure_constants[["k"]][i] * (delta - 2) / 2 + structure_constants[["dim_omega"]][i] / structure_constants[["r"]][i]

    calculate_gamma_omega(lambda_i, structure_constants[["dim_omega"]][i], structure_constants[["r"]][i], structure_constants[["d"]][i])
  })

  sum(single_G_i)
}

#' Calculate log phi_part of log_posteriori_of_gips
#'
#' @param structure_constants An output of
#' `gips::get_structure_constants(perm_proposal, perm_size)`.
#' Rest of params as in `log_posteriori_of_gips()`.
#'
#' @noRd
calculate_phi_part <- function(perm_proposal, number_of_observations, U,
                               delta, D_matrix, structure_constants) {
  # projection of matrices on perm_proposal
  equal_indices <- get_equal_indices_by_perm(perm_proposal)
  Dc <- gips::project_matrix(D_matrix, perm_proposal,
    precomputed_equal_indices = equal_indices
  )
  Uc <- gips::project_matrix(U, perm_proposal,
    precomputed_equal_indices = equal_indices
  )

  # divide by 2 - refer to newest version of the paper
  Dc <- Dc / 2
  Uc <- Uc / 2

  # diagonalization
  diagonalising_matrix <- gips::prepare_orthogonal_matrix(perm_proposal)
  Dc_diagonalised <- t(diagonalising_matrix) %*% Dc %*% diagonalising_matrix
  DcUc_diagonalised <- t(diagonalising_matrix) %*% (Uc + Dc) %*% diagonalising_matrix

  # block part
  block_ends <- get_block_ends(structure_constants)
  Dc_block_log_dets <- calculate_log_determinants_of_block_matrices(
    Dc_diagonalised,
    block_ends
  )
  DcUc_block_log_dets <- calculate_log_determinants_of_block_matrices(
    DcUc_diagonalised,
    block_ends
  )
  Dc_exponent <- (delta - 2) / 2 + structure_constants[["dim_omega"]] /
    (structure_constants[["r"]] * structure_constants[["k"]])
  DcUc_exponent <- -(number_of_observations + delta - 2) / 2 - structure_constants[["dim_omega"]] /
    (structure_constants[["r"]] * structure_constants[["k"]])

  out <- sum(Dc_block_log_dets * Dc_exponent + DcUc_block_log_dets * DcUc_exponent)

  out
}

#' Calculate the logarithm of a single Gamma omega function
#'
#' Using the formula (12) from the paper
#'
#' @inheritParams calculate_gamma_function
#' @param dim_omega_i Single element from `gips::get_structure_constants`.
#' @param r_i Single element from `gips::get_structure_constants`.
#' @param d_i Single element from `gips::get_structure_constants`.
#'
#' @returns Logarithm of the value of Gamma function.
#'
#' @noRd
calculate_gamma_omega <- function(lambda, dim_omega_i, r_i, d_i) {
  if (lambda <= dim_omega_i / r_i - 1) {
    rlang::warn(c("Gamma integral is divergent for the given lambda value and structure constants.",
      "i" = paste0(
        "Gamma(lambda = ", lambda,
        ", dim_omega_i = ", dim_omega_i,
        ", r_i = ", r_i,
        ", d_i = ", d_i,
        ") = Inf."
      )
    ))
    return(Inf) # the integral does not converge
  }

  sum(lgamma((0:(-(r_i - 1))) * d_i / 2 + lambda)) + (dim_omega_i - r_i) / 2 * log(2 * pi)
}

#' Calculate log of determinants of matrices from block decomposition
#'
#' Block decomposition 1 from paper
#'
#' @param diagonalized_matrix A middle matrix from decomposition 1.
#' @param block_ends The indices of last columns of block matrices.
#' Last element equals size of matrix.
#'
#' @returns A numeric vector.
#' @noRd
calculate_log_determinants_of_block_matrices <- function(diagonalised_matrix,
                                                         block_ends) {
  block_starts <- c(0, block_ends[-length(block_ends)] + 1)
  sapply(1:length(block_starts), function(i) {
    slice <- block_starts[i]:block_ends[i]
    block_matrix <- diagonalised_matrix[slice, slice, drop = FALSE]
    determinant(block_matrix, logarithm = TRUE)[["modulus"]]
  })
}

#' Get coprime numbers
#'
#' @param n A single integer.
#'
#' @returns All integers smaller than n, that are coprime to n. Exception: n = 1,
#' in which case 1 is returned.
#' @noRd
get_coprimes <- function(n) {
  if (n == 1) {
    return(1)
  }
  smaller_ints <- 1:(n - 1)
  are_coprime <- sapply(smaller_ints, function(m) numbers::coprime(n, m))
  smaller_ints[are_coprime]
}

#' Get indices of elements of perm_size x perm_size matrix, which should be equal
#'
#' @param perm An object of a `gips_perm` class.
#'
#' @returns A list of integer vectors. Each vector contains a SINGLE index
#' of elements, which should be equal in symmetrical matrix invariant
#' by permutation `perm`.
#'
#' @examples
#' perm <- gips_perm("(1,2,3)(4,5)", 6)
#' matrix_symvariant <- matrix(c(
#'   2, 1, 1, 3, 3, 4,
#'   1, 2, 1, 3, 3, 4,
#'   1, 1, 2, 3, 3, 4,
#'   3, 3, 3, 5, 6, 7,
#'   3, 3, 3, 6, 5, 7,
#'   4, 4, 4, 7, 7, 8
#' ), byrow = TRUE, ncol = 6)
#' out <- get_equal_indices_by_perm(perm, 6)
#' all(sapply(out, function(v) all.equal(matrix_symvariant[v]))) # TRUE
#' @noRd
get_equal_indices_by_perm <- function(perm) {
  perm_size <- attr(perm, "size")
  # We'll be iterating over pairs of subcycles
  subcycle_indice_pairs <- matrix(
    c(
      rep(1:length(perm), each = length(perm)),
      rep(1:length(perm), times = length(perm))
    ),
    ncol = 2
  )

  subcycle_indice_pairs <- subcycle_indice_pairs[subcycle_indice_pairs[, 1] <= subcycle_indice_pairs[, 2], ,
    drop = FALSE
  ]
  # subcycle_indice_pairs is a matrix of pairs like (5,6), where the second is not smaller than the first

  # Let's go
  nested_list <- lapply(1:nrow(subcycle_indice_pairs), function(pair_index) {
    i <- subcycle_indice_pairs[pair_index, 1]
    j <- subcycle_indice_pairs[pair_index, 2]
    subcycle_1 <- perm[[i]]
    subcycle_2 <- perm[[j]]

    # matrix_subcycle is a subcycle of permutation P defined as
    # P(k,l) = (perm(k), perm(l))
    matrix_subcycle_length <- numbers::LCM(
      length(subcycle_1),
      length(subcycle_2)
    )
    number_of_matrix_subcycles <- length(subcycle_1) * length(subcycle_2) /
      matrix_subcycle_length

    # Instead of operating on subcycle elements, we will be operating
    # on elements' indices
    # I.e. for subcycle s=[3,5,2] we operate on indices [1,2,3]
    # s[1] = 3 etc
    elements_indices <- 1:matrix_subcycle_length
    subcycle_1_indices <- elements_indices %% length(subcycle_1)
    subcycle_1_indices[subcycle_1_indices == 0] <- length(subcycle_1)
    subcycle_1_elements <- subcycle_1[subcycle_1_indices]

    subcycle_2_indices <- elements_indices %% length(subcycle_2)
    subcycle_2_indices[subcycle_2_indices == 0] <- length(subcycle_2)

    lapply(1:(number_of_matrix_subcycles), function(k) {
      subcycle_2_elements <- subcycle_2[shift_vector(
        subcycle_2_indices,
        k - 1
      )]

      double_indices <- matrix(c(
        subcycle_1_elements, subcycle_2_elements,
        subcycle_2_elements, subcycle_1_elements
      ), ncol = 2)

      single_indices <- get_single_from_double_indices(
        double_indices,
        perm_size
      )

      # We need to correct for matrix subcycles, that are symmetric
      if (i == j &&
        single_indices[matrix_subcycle_length + 1] %in%
          single_indices[1:matrix_subcycle_length]) {
        single_indices <- single_indices[1:matrix_subcycle_length]
      }
      shift_vector(single_indices, which.min(single_indices) - 1)
    })
  })
  unlist(nested_list, recursive = FALSE)
}

#' Rearrange vector
#'
#' Move elements from the start of a vector to its end, so that the minimal
#' element will be first.
#'
#' @examples
#' v <- c(5, 3, 2, 1, 4)
#' rearranged <- rearrange_vector(v)
#' all(rearranged == c(1, 4, 5, 3, 2)) # TRUE
#'
#' @noRd
rearrange_vector <- function(v) {
  shift_vector(v, which.min(v) - 1)
}

get_single_from_double_indices <- function(indices, matrix_size) {
  (indices[, 2] - 1) * matrix_size + indices[, 1]
}
