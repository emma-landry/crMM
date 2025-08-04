kfeature_aUpdate <- function(t, y, c, a, phi, rho, tt_basis, gamma, pi, knots_shape, degree = 3,
                             intercept = F, var_e, common_a = F) {
  n <- length(t)
  N <- length(c)
  K <- ncol(pi)

  a_old <- a

  if (common_a == F) {
    a_new <- a_old
    for (k in 1:K){
      a_new[k] <- extraDistr::rtnorm(n = 1, mean = a_old[k], sd = 0.1, a = 0)
    }


    P1 <- 0
    P0 <- 0

    for (i in 1:N) {
      if (K > 1){
        modelMean_old <- kfeature_meanWarp(t, c[i], a_old, phi[i, ], rho, tt_basis, gamma, pi[i, ],
                                           knots_shape, degree, intercept)
        modelMean_new <- kfeature_meanWarp(t, c[i], a_new, phi[i, ], rho, tt_basis, gamma, pi[i, ],
                                           knots_shape, degree, intercept)
      } else {
        modelMean_old <- kfeature_meanWarp(t, c[i], a_old, phi[i, ], rho, tt_basis, gamma, pi[i],
                                           knots_shape, degree, intercept)
        modelMean_new <- kfeature_meanWarp(t, c[i], a_new, phi[i, ], rho, tt_basis, gamma, pi[i],
                                           knots_shape, degree, intercept)
      }

      y_i <- y[i, ]
      P1 <- P1 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_new, sigma = var_e * diag(n), log = T)
      P0 <- P0 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_old, sigma = var_e * diag(n), log = T)
    }

    Q1 <- 0
    Q0 <- 0

    for (k in 1:K){
      P1 <- P1 + stats::dlnorm(a_new[k], meanlog = 0, sdlog = 0.5, log = T)
      P0 <- P0 + stats::dlnorm(a_old[k], meanlog = 0, sdlog = 0.5, log = T)

      Q1 <- Q1 + extraDistr::dtnorm(x = a_new[k], mean = a_old[k], sd = 0.1, a = 0, log = T)
      Q0 <- Q0 + extraDistr::dtnorm(x = a_old[k], mean = a_new[k], sd = 0.1, a = 0, log = T)
    }

    ratio <- (P1 - Q1) - (P0 - Q0)
    if (log(stats::runif(1)) < ratio) {
      a <- a_new
    }
  } else {
    a_new <- extraDistr::rtnorm(n = 1, mean = a_old[1], sd = 0.1, a = 0)
    a_new <- rep(a_new, K)

    P1 <- 0
    P0 <- 0

    for (i in 1:N) {
      if (K > 1){
        modelMean_old <- kfeature_meanWarp(t, c[i], a_old, phi[i, ], rho, tt_basis, gamma, pi[i, ],
                                           knots_shape, degree, intercept)
        modelMean_new <- kfeature_meanWarp(t, c[i], a_new, phi[i, ], rho, tt_basis, gamma, pi[i, ],
                                           knots_shape, degree, intercept)
      } else {
        modelMean_old <- kfeature_meanWarp(t, c[i], a_old, phi[i, ], rho, tt_basis, gamma, pi[i],
                                           knots_shape, degree, intercept)
        modelMean_new <- kfeature_meanWarp(t, c[i], a_new, phi[i, ], rho, tt_basis, gamma, pi[i],
                                           knots_shape, degree, intercept)
      }

      y_i <- y[i, ]
      P1 <- P1 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_new, sigma = var_e * diag(n), log = T)
      P0 <- P0 + mvtnorm::dmvnorm(x = y_i, mean = modelMean_old, sigma = var_e * diag(n), log = T)
    }

    P1 <- P1 + stats::dlnorm(a_new[1], meanlog = 0, sdlog = 0.5, log = T)
    P0 <- P0 + stats::dlnorm(a_old[1], meanlog = 0, sdlog = 0.5, log = T)

    Q1 <- extraDistr::dtnorm(x = a_new[1], mean = a_old[1], sd = 0.1, a = 0, log = T)
    Q0 <- extraDistr::dtnorm(x = a_old[1], mean = a_new[1], sd = 0.1, a = 0, log = T)
    ratio <- (P1 - Q1) - (P0 - Q0)
    if (log(stats::runif(1)) < ratio) {
      a <- a_new
    }
    }

  return(a)
}
