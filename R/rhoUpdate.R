rhoUpdate <- function(t, y, c, phi, rho, tt_basis, pi, gamma1, gamma2, knots_shape,
                      degree = 3, intercept = F, var_e) {

  n <- length(t)
  N <- length(c)
  Q <- ncol(phi)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  P1 <- 0
  P0 <- 0
  rho_old <- rho
  # rho_new <- stats::rbeta(n = 1, shape1 = 1 / (1 - rho_old), shape2 = 2)
  rho_new <- extraDistr::rtnorm(n = 1, mean = rho_old, sd = 0.01, a = 0, b = Inf)

  for (i in 1:N) {
    modelMean_old <- meanWarp(t, c[i], phi[i, ], rho_old, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp(t, c[i], phi[i, ], rho_new, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)

    y_i <- y[i, ]
    P1 <- P1 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_new, sigma = var_e * diag(n), log = T)
    P0 <- P0 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_old, sigma = var_e * diag(n), log = T)
  }

  # Q1 <- stats::dbeta(x = rho_new, shape1 = 1/ (1 - rho_old), shape2 = 2, log = T)
  # Q0 <- stats::dbeta(x = rho_old, shape1 = 1/ (1 - rho_new), shape2 = 2, log = T)

  # P1 <- P1 + stats::dbeta(rho_new, shape1 = 1 / 97 , shape2 = 1 / 97 , log = T)
  # P0 <- P0 + stats::dbeta(rho_old, shape1 = 1 / 97, shape2 = 1 / 97, log = T)

  P1 <- P1 + stats::dgamma(rho_new, shape = 0.5 , rate = 0.5 , log = T)
  P0 <- P0 + stats::dgamma(rho_old, shape = 0.5, rate = 0.5, log = T)

  Q1 <- extraDistr::dtnorm(x = rho_new, mean = rho_old, sd = 0.01, a = 0, b = Inf, log = T)
  Q0 <- extraDistr::dtnorm(x = rho_old, mean = rho_new, sd = 0.01, a = 0, b = Inf, log = T)

  ratio <- (P1 - Q1) - (P0 - Q0)
  if (log(stats::runif(1)) < ratio) {
    rho <- rho_new
  }
  return(rho)
}

kfeature_rhoUpdate <- function(t, y, c, a, phi, rho, tt_basis, pi, gamma, knots_shape,
                      degree = 3, intercept = F, var_e) {

  n <- length(t)
  N <- length(c)
  Q <- ncol(phi)
  warp_num <- length(rho)
  K <- ncol(gamma)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  P1 <- 0
  P0 <- 0
  rho_old <- rho

  rho_new <- rho_old
  for (w in 1:warp_num){
    rho_new[w] <- extraDistr::rtnorm(n = 1, mean = rho_old[w], sd = 0.01, a = 0, b = 1)
  }


  for (i in 1:N) {
    if (K > 1){
      modelMean_old <- kfeature_meanWarp(t, c[i], a, phi[i, ], rho_old, tt_basis, gamma, pi[i, ],
                                         knots_shape, degree, intercept)
      modelMean_new <- kfeature_meanWarp(t, c[i], a,phi[i, ], rho_new, tt_basis, gamma, pi[i, ],
                                         knots_shape, degree, intercept)
    } else {
      modelMean_old <- kfeature_meanWarp(t, c[i], a, phi[i, ], rho_old, tt_basis, gamma, pi[i],
                                         knots_shape, degree, intercept)
      modelMean_new <- kfeature_meanWarp(t, c[i], a, phi[i, ], rho_new, tt_basis, gamma, pi[i],
                                         knots_shape, degree, intercept)
    }

    y_i <- y[i, ]
    P1 <- P1 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_new, sigma = var_e * diag(n), log = T)
    P0 <- P0 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_old, sigma = var_e * diag(n), log = T)
  }

  # Q1 <- stats::dbeta(x = rho_new, shape1 = 1/ (1 - rho_old), shape2 = 2, log = T)
  # Q0 <- stats::dbeta(x = rho_old, shape1 = 1/ (1 - rho_new), shape2 = 2, log = T)

  # P1 <- P1 + stats::dbeta(rho_new, shape1 = 1 / 97 , shape2 = 1 / 97 , log = T)
  # P0 <- P0 + stats::dbeta(rho_old, shape1 = 1 / 97, shape2 = 1 / 97, log = T)

  Q1 <- 0
  Q0 <- 0

  for (w in 1:warp_num) {
    P1 <- P1 + stats::dbeta(rho_new[w], shape1 = .5 , shape2 = .5 , log = T)
    P0 <- P0 + stats::dbeta(rho_old[w], shape1 = .5, shape2 =  .5, log = T)

    Q1 <- Q1 + extraDistr::dtnorm(x = rho_new[w], mean = rho_old[w], sd = 0.01, a = 0, b = 1, log = T)
    Q0 <- Q0 + extraDistr::dtnorm(x = rho_old[w], mean = rho_new[w], sd = 0.01, a = 0, b = 1, log = T)
  }

  ratio <- (P1 - Q1) - (P0 - Q0)
  if (log(stats::runif(1)) < ratio) {
    rho <- rho_new
  }
  return(rho)
}
