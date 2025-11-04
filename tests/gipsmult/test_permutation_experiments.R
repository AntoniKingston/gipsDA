correct_perm_pro <- function(counts, correct_perm, tot) {
  if (!(correct_perm %in% names(counts))){
    return(0)
  }
  return(counts[[correct_perm]]/tot)
}

mat_gen <- function(m, nu, V) {
  Ss <- rWishart(m, nu, V) / nu
  Ss <- lapply(seq_len(dim(Ss)[3]), function(i) Ss[,,i])
  Ss
}

avg_list_of_mat <- function(Ss) {
  suma <- Ss[[1]]
  for (S in Ss[-1]) {
    suma <- suma + S
  }
  suma / length(Ss)
}

single_experiment <- function(n, m, nu, V, optimizer, max_iter, perm) {
  # Generowanie danych jest teraz losowe dla każdego wywołania, co jest poprawne
  Sss <- replicate(n = n, expr = mat_gen(m, nu, V), simplify = FALSE)
  perms <- list()
  perms_gips <- list()
  perms_ggips <- list()

  for (i in 1:n) {
    Ss <- Sss[[i]]
    g <- gipsmult:::gipsmult(Ss, nu, was_mean_estimated = FALSE)
    gg <- gips::gips(Ss[[1]], nu, was_mean_estimated = FALSE)
    ggg <- gips::gips(avg_list_of_mat(Ss), nu*m, was_mean_estimated = FALSE)

    g_map <- gipsmult:::find_MAP(g, optimizer = optimizer, max_iter = max_iter, show_progress_bar = FALSE)
    gg_map <- gips::find_MAP(gg, optimizer = optimizer,  max_iter = max_iter, show_progress_bar = FALSE)
    ggg_map <- gips::find_MAP(ggg, optimizer = optimizer, max_iter = max_iter, show_progress_bar = FALSE)

    perms <- c(perms, as.character(g_map[[1]]))
    perms_gips <- c(perms_gips, as.character(gg_map[[1]]))
    perms_ggips <- c(perms_ggips, as.character(ggg_map[[1]]))
  }

  counts <- sort(table(unlist(perms)), decreasing = TRUE)
  counts_gips <- sort(table(unlist(perms_gips)), decreasing = TRUE)
  counts_ggips <- sort(table(unlist(perms_ggips)), decreasing = TRUE)

  return(c(gipsmult = correct_perm_pro(counts, perm, n),
           gips_single = correct_perm_pro(counts_gips, perm, n),
           gips_avg = correct_perm_pro(counts_ggips, perm, n)))
}


# ==============================================================================
# 2. Główna funkcja eksperymentu (MAŁA, ALE WAŻNA ZMIANA)
# ==============================================================================

run_permutation_experiment <- function(n, m, n_nus, nu_lb, nu_ub,  dim, optimizer, max_iter, perm) {

  lwblg <- log(nu_lb)
  upblg <- log(nu_ub)
  nus <- round(exp(seq(lwblg, upblg, length.out = n_nus)))

  # Generowanie V jest teraz unikalne dla każdego eksperymentu, co jest bardziej poprawne metodologicznie.
  # Ponieważ ziarno jest zarządzane na zewnątrz, wyniki wciąż będą powtarzalne.
  X <- matrix(rnorm(dim * dim), nrow = dim)
  Q <- qr.Q(qr(X))
  lambdas <- rgamma(dim, shape = 2, rate = 1)
  V <- Q %*% diag(lambdas) %*% t(Q)
  V <- gips::project_matrix(V, perm)

  # Ta pętla po 'nu' jest wystarczająco szybka, zrównoleglamy pętlę zewnętrzną
  message(paste("Uruchamiam eksperyment dla nu =", nus[[1]]))
  props <- single_experiment(n, m, nus[[1]], V, optimizer, max_iter, perm)

  for (nu in nus[-1]) {
    message(paste("Uruchamiam eksperyment dla nu =", nu))
    props <- rbind(props, single_experiment(n, m, nu, V, optimizer, max_iter, perm))
  }

  rownames(props) <- nus

  # Dostosuj ścieżkę do swoich potrzeb
  folder_docelowy <- "repo/gipsDA/tests/gipsmult/test_csv/"
  if (!dir.exists(folder_docelowy)) {
    dir.create(folder_docelowy, recursive = TRUE)
  }

  perm_filename <- gsub("[^a-zA-Z0-9_]", "", perm)
  perm_filename <- substr(perm_filename, 1, 15)

  nazwa_pliku <- paste0(
    "m-", m,
    "_dim-", dim,
    "_opt-", optimizer,
    "_iter-", max_iter,
    "_perm-", perm_filename,
    ".csv"
  )

  sciezka_pliku <- file.path(folder_docelowy, nazwa_pliku)
  write.csv(props, file = sciezka_pliku, row.names = TRUE)
  message(paste("Wyniki zapisano do pliku:", sciezka_pliku))

  return(props)
}

# ==============================================================================
# 3. Definicja parametrów (bez zmian)
# ==============================================================================
n <- 10
n_nus <- 10
nu_lb <- 50
nu_ub <- 1000

perm_dim5_1 <- "(2,3,4)"
perm_dim5_2 <- "(1,3)(2,5)"
perm_dim50_1 <- "(1,14)(2,8,9,37)(3,43,34,21,28)(4,19,39,36)(5,27)(7,41,10,15)(11,23,20,50,40)(13,33,38,45)(16,48)(17,18,35,46,44)(22,42)(24,31)(25,49,29)(26,32)"
perm_dim50_2 <- "(1,49,40,42,11,19)(2,39,22)(3,38,45,48)(4,7,16,47,37,15)(5,27,32,46,28)(6,14,24,20,10,17)(8,12,43,31,41,36)(9,21,34)(18,33)(23,26,30)(25,29,44)"

params <- expand.grid(
  m = c(5, 50),
  dim = c(5, 50),
  optimizer = c("BF", "MH", "HC"),
  max_iter = c(100, 500, 1000),
  perm = c("perm1", "perm2"),
  stringsAsFactors = FALSE
)

params <- params[!(params$dim == 5 & params$optimizer %in% c("MH", "HC")),]
params <- params[!(params$dim == 50 & params$optimizer == "BF"),]
params <- params[!(params$optimizer != "MH" & params$max_iter != 100),]
params$perm[params$dim == 5 & params$perm == "perm1"] <- perm_dim5_1
params$perm[params$dim == 5 & params$perm == "perm2"] <- perm_dim5_2
params$perm[params$dim == 50 & params$perm == "perm1"] <- perm_dim50_1
params$perm[params$dim == 50 & params$perm == "perm2"] <- perm_dim50_2

params <- unique(params)
rownames(params) <- NULL

# ==============================================================================
# 4. Uruchomienie wszystkich kombinacji W SPOSÓB RÓWNOLEGŁY
# ==============================================================================

# Krok 1: Wczytaj biblioteki do przetwarzania równoległego
library(future)
library(future.apply)

# Krok 2: Ustaw strategię równoległą.
# 'multisession' uruchomi tyle procesów R w tle, ile masz rdzeni.
# Możesz też ręcznie ustawić liczbę, np. plan(multisession, workers = 4)
plan(multisession)
message(paste("Uruchamiam obliczenia równoległe na", future::availableCores(), "rdzeniach."))

# Krok 3: Zastąp pętlę 'for' przez 'future_lapply'
# Ta funkcja wykona operacje dla każdego wiersza 'params' na osobnym rdzeniu.
# 'future.seed = TRUE' gwarantuje powtarzalność wyników losowych.
future_lapply(1:nrow(params), function(i) {

  current_params <- params[i, ]

  # Komunikaty będą się teraz mieszać, bo wiele procesów pisze na raz.
  # To normalne zachowanie.
  message(paste("\n--- Rozpoczynam przebieg", i, "z", nrow(params), "na jednym z rdzeni ---"))
  message(paste(names(current_params), current_params, collapse = " | "))

  run_permutation_experiment(
    n = n,
    m = current_params$m,
    n_nus = n_nus,
    nu_lb = nu_lb,
    nu_ub = nu_ub,
    dim = current_params$dim,
    optimizer = current_params$optimizer,
    max_iter = current_params$max_iter,
    perm = current_params$perm
  )
}, future.seed = TRUE) # <-- KLUCZOWY ELEMENT DLA POWTARZALNOŚCI!

message("\n--- Wszystkie eksperymenty zostały zlecone i zakończone! ---")

# Zamknij procesy w tle po zakończeniu pracy (dobra praktyka)
plan(sequential)