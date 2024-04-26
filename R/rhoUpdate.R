rhoUpdate <- function(t, y, c, phi, rho, tt_basis, gamma1, gamma2, knots_shape,
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
  rho_new <- stats::rbeta(n = 1, shape1 = 1 / (1 - rho_old), shape2 = 2)

  for (i in 1:N) {
    modelMean_old <- meanWarp(t, c[i], phi[i, ], rho_old, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp(t, c[i], phi[i, ], rho_new, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)

    y_i <- y[i, ]
    P1 <- P1 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_new, sigma = var_e * diag(n), logged = T)
    P0 <- P0 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_old, sigma = var_e * diag(n), logged = T)
  }

  Q1 <- stats::dbeta(x = rho_new, shape1 = 1/ (1 - rho_old), shape2 = 2, log = T)
  Q0 <- stats::dbeta(x = rho_old, shape1 = 1/ (1 - rho_new), shape2 = 2, log = T)

  ratio <- (P1 - Q1) - (P0 - Q0)
  if (log(stats::runif(1)) < ratio) {
    rho <- rho_new
  }
  return(rho)
}
