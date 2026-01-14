project_covs <- function(emp_covs, ns_obs, MAP = TRUE, optimizer, max_iter, tol = 1e-3) {
    gg <- gipsmult(emp_covs, ns_obs, was_mean_estimated = TRUE)
    if (MAP) {
        gg <- find_MAP(gg, optimizer = optimizer, max_iter = max_iter, show_progress_bar = FALSE)
        perm <- gg[[1]]
        return(list(covs = lapply(emp_covs, function(x) gips::project_matrix(x, perm)), opt_info = perm))
    }
    gg <- find_MAP(gg, optimizer = optimizer, max_iter = max_iter, return_probabilities = TRUE, save_all_perms = TRUE, show_progress_bar = FALSE)
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



project_matrix_multiperm <- function(emp_cov, probs) {
    perms <- names(probs)
    projected_matrix <- matrix(0, nrow = dim(emp_cov), dim(emp_cov))
    for (i in 1:length(probs)) {
        projected_matrix <- projected_matrix + probs[[i]] * gips::project_matrix(emp_cov, perms[i])
    }
    return(projected_matrix / sum(probs))
}



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

gipsDA_to_json <- function(obj, file) {
  jsonlite::write_json(
    serialize_for_json(obj),
    file,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

gipsDA_from_json <- function(file, classname) {
  raw <- jsonlite::read_json(file)
  obj <- deserialize_from_json(raw)
  class(obj) <- classname
  obj
}

recursive_length <- function(x) {
  if (is.atomic(x)) {
    return(length(x))
  }
  if (is.list(x)) {
    return(sum(vapply(x, recursive_length, integer(1))))
  }
  return(0)
}

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


