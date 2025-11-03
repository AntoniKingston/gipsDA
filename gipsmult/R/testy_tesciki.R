
n <- 10 #number of repetitions for each nu
m <- 7 #number of classes for gipsmult scenario
dim <- 5 #diemnsionality of data
n_nus <- 10 #number of different numbers of observations
nu_lb <- 50
nu_ub <- 5000 # lower and upper bounds for numbers of observations
lwblg <- log(nu_lb)
upblg <- log(nu_up)
nus <- round(exp(seq(lwblg, upblg, length.out = n_nus))) # final vector of numbers of observations (spread out logarithmicaly)
perm <- "(2,3,4)"
optimizer <- "BF"
max_iter <- 100
set.seed(12345)

#Generating theoretical covariance matrix
X <- matrix(rnorm(dim * dim), nrow = dim)
Q <- qr.Q(qr(X))
lambdas <- rgamma(dim, shape = 2, rate = 1)
V <- Q %*% diag(lambdas) %*% t(Q)
V <- gips::project_matrix(V, perm)

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
  Sss <- replicate(n = n, expr = mat_gen(m, nu, V), simplify = FALSE)
  perms <- list()
  perms_gips <- list()
  perms_ggips <- list()
  for (i in 1:n) {
    Ss <- Sss[[i]]
    g <- gipsmult:::gipsmult(Ss, nu, was_mean_estimated = FALSE)
    gg <- gips::gips(Ss[[1]], nu, was_mean_estimated = FALSE)
    ggg <- gips::gips(avg_list_of_mat(Ss), nu*m, was_mean_estimated = FALSE)
    g_map <- gipsmult:::find_MAP(g, optimizer = optimizer, max_iter = max_iter)
    gg_map <- gips::find_MAP(gg, optimizer = optimizer,  max_iter = max_iter)
    gg_map <- gips::find_MAP(ggg, optimizer = optimizer, max_iter = max_iter)
    perms <- c(perms, as.character(g_map[[1]]))
    perms_gips <- c(perms_gips, as.character(gg_map[[1]]))
    perms_ggips <- c(perms_ggips, as.character(gg_map[[1]]))
  }
  counts <- sort(table(unlist(perms)), decreasing = TRUE)
  counts_gips <- sort(table(unlist(perms_gips)), decreasing = TRUE)
  counts_ggips <- sort(table(unlist(perms_ggips)), decreasing = TRUE)
  return(c(correct_perm_pro(counts, perm, n), correct_perm_pro(counts_gips, perm, n), correct_perm_pro(counts_ggips, perm, n)))
}



props <- single_experiment(n,m,nus[[1]],V, optimizer, max_iter, perm)

for (nu in nus[-1]) {
  props <- rbind(props, single_experiment(n, m, nu, V, optimizer, max_iter, perm))
}

#na wejściu to co na wejściu
#na wyjściu macierz len(nus)x3
#1 kolumna - gipsmult
#2 kolumna - gips (na 1 macierzy)
#3 kolumna - gips na macierzy uśrednionej z wszystkich macierzy
#5 kolumna - ...
#wyjście do wyplotowania w zaimplementowany przez ciebie lub przez dowolny model sposób
pojedynczy_bardzo_kurwa_indywidualny_zajebiscie_wyekstrahowany_ekperyment_zwany_rowniez_testem_jak_bedziesz_chcial_zmienic_nazwe_tej_funkcji_to_przestaje_pisac_inzynierke_i_zostajesz_sam_z_chojem_ghaagaghaghah <- function (n, m, n_nus, nu_lb, nu_ub,  dim, optimizer, max_iter, perm) {
  #liczby obserwacji po logarytmie (prawie jak po roku w rosji)
  lwblg <- log(nu_lb)
  upblg <- log(nu_ub)
  nus <- round(exp(seq(lwblg, upblg, length.out = n_nus))) # final vector of numbers of observations (spread out logarithmicaly)
  #twoja mać
  X <- matrix(rnorm(dim * dim), nrow = dim)
  Q <- qr.Q(qr(X))
  lambdas <- rgamma(dim, shape = 2, rate = 1)
  V <- Q %*% diag(lambdas) %*% t(Q)
  V <- gips::project_matrix(V, perm)
  #jegzberymendy
  props <- single_experiment(n,m,nus[[1]],V, optimizer, max_iter, perm)
  for (nu in nus[-1]) {
    props <- rbind(props, single_experiment(n, m, nu, V, optimizer, max_iter, perm))
  }
  props
}
#rób co chcesz z maciorą
plotuj_macierz_3naen <- function(maciora) {
  #nie uważałem na TWD to powiem szczerze, że akurat w R to zmówię niebagatelnie długie pacierze i niewyplotowane zostaną macierze

}
