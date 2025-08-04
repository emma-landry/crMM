piUpdate_NoWarp <- function(t, y, c, gamma1, gamma2, pi, shape_basis,
                            var_e, alpha, rescale = T, label1 = NULL, label2 = NULL,
                            tuning_param = 1000) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (!is.null(label1)) {
    if (!(all(label1 >= 1 & label1 <= N))) {
      stop("The indexes in 'label1' fall out of the range of the total number of observations.")
    }

    if (!is.null(label2)) {
      if (length(intersect(label1, label2)) > 0) {
        stop("Observations cannot be labelled as belonging fully to both feature 1 and feature 2.")
      }
    }
  }

  if (!is.null(label2)) {
    if (!(all(label2 >= 1 & label2 <= N))) {
      stop("The indexes in 'label2' fall out of the range of the total number of observations.")
    }
  }

  for (i in 1:N) {
    if (!is.null(label1) & (i %in% label1)) {
      next
    }

    if (!is.null(label2) & (i %in% label2)) {
      next
    }

    pi_old <- pi[i, ]

    if(pi_old[1] < 0.0001) {
      pi_old[1] <- 0.01
      pi_old[2] <- 1 - pi_old[1]
    }

    if (pi_old[2] < 0.0001) {
      pi_old[2] <- 0.01
      pi_old[1] <- 1 - pi_old[2]
    }

    pi_new <- LaplacesDemon::rdirichlet(n = 1, alpha = tuning_param * pi_old)

    if(pi_new[1] < 0.0001) {
      pi_new[1] <- 0.01
      pi_new[2] <- 1 - pi_new[1]
    }

    if (pi_new[2] < 0.0001) {
      pi_new[2] <- 0.01
      pi_new[1] <- 1 - pi_new[2]
    }

    modelMean_old <- meanNoWarp(t, c[i], gamma1, gamma2, pi_old, shape_basis)
    modelMean_new <- meanNoWarp(t, c[i], gamma1, gamma2, pi_new, shape_basis)
    y_i <- y[i, ]

    P1 <- -1 / (2 * var_e) * sum((y_i - modelMean_new) ^ 2) +
          (alpha[1] - 1) * log(pi_new[1]) +
          (alpha[2] - 1) * log(pi_new[2])
    Q1 <- LaplacesDemon::ddirichlet(pi_new, alpha = tuning_param * pi_old, log = T)

    P0 <- -1 / (2 * var_e) * sum((y_i - modelMean_old) ^ 2) +
          (alpha[1] - 1) * log(pi_old[1]) +
          (alpha[2] - 1) * log(pi_old[2])
    Q0 <- LaplacesDemon::ddirichlet(pi_old, alpha = tuning_param * pi_new, log = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      pi_old <- pi_new
    }
    pi[i, ] <- pi_old
  }

  if (rescale == T) {
    pi[, 1] <- (pi[, 1] - min(pi[, 1])) / (max(pi[, 1]) - min(pi[, 1]))
    pi[, 2] <- 1 - pi[, 1]
  }
  return(pi)
}

piUpdate_Warp <- function(t, y, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree = 3,
                          intercept = F, var_e, alpha, rescale = T, label1 = NULL, label2 = NULL,
                          tuning_param = 1000) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (!is.null(label1)) {
    if (!(all(label1 >= 1 & label1 <= N))) {
      stop("The indexes in 'label1' fall out of the range of the total number of observations.")
    }

    if (!is.null(label2)) {
      if (length(intersect(label1, label2)) > 0) {
        stop("Observations cannot be labelled as belonging fully to both feature 1 and feature 2.")
      }
    }
  }

  if (!is.null(label2)) {
    if (!(all(label2 >= 1 & label2 <= N))) {
      stop("The indexes in 'label2' fall out of the range of the total number of observations.")
    }
  }

  for (i in 1:N) {
    if (!is.null(label1) & (i %in% label1)) {
      next
    }

    if (!is.null(label2) & (i %in% label2)) {
      next
    }

    pi_old <- pi[i, ]

    if(pi_old[1] < 0.0001) {
      pi_old[1] <- 0.01
      pi_old[2] <- 1 - pi_old[1]
    }

    if (pi_old[2] < 0.0001) {
      pi_old[2] <- 0.01
      pi_old[1] <- 1 - pi_old[2]
    }

    pi_new <- LaplacesDemon::rdirichlet(n = 1, alpha = tuning_param * pi_old)

    if(pi_new[1] < 0.0001) {
      pi_new[1] <- 0.01
      pi_new[2] <- 1 - pi_new[1]
    }

    if (pi_new[2] < 0.0001) {
      pi_new[2] <- 0.01
      pi_new[1] <- 1 - pi_new[2]
    }

    modelMean_old <- meanWarp(t, c[i], phi[i, ], rho, tt_basis, gamma1, gamma2, pi_old,
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp(t, c[i], phi[i, ], rho, tt_basis, gamma1, gamma2, pi_new,
                              knots_shape, degree, intercept)
    y_i <- y[i, ]

    P1 <- -1 / (2 * var_e) * sum((y_i - modelMean_new) ^ 2) +
      (alpha[1] - 1) * log(pi_new[1]) +
      (alpha[2] - 1) * log(pi_new[2])
    Q1 <- LaplacesDemon::ddirichlet(pi_new, alpha = tuning_param * pi_old, log = T)

    P0 <- -1 / (2 * var_e) * sum((y_i - modelMean_old) ^ 2) +
      (alpha[1] - 1) * log(pi_old[1]) +
      (alpha[2] - 1) * log(pi_old[2])
    Q0 <- LaplacesDemon::ddirichlet(pi_old, alpha = tuning_param * pi_new, log = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      pi_old <- pi_new
    }
    pi[i, ] <- pi_old
  }

  if (rescale == T) {
    pi[, 1] <- (pi[, 1] - min(pi[, 1])) / (max(pi[, 1]) - min(pi[, 1]))
    pi[, 2] <- 1 - pi[, 1]
  }
  return(pi)
}

piUpdate_Warp_alt <- function(t, y, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape, degree = 3,
                              intercept = F, var_e, alpha, rescale = T, label1 = NULL, label2 = NULL,
                              tuning_param = 1000) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (!is.null(label1)) {
    if (!(all(label1 >= 1 & label1 <= N))) {
      stop("The indexes in 'label1' fall out of the range of the total number of observations.")
    }

    if (!is.null(label2)) {
      if (length(intersect(label1, label2)) > 0) {
        stop("Observations cannot be labelled as belonging fully to both feature 1 and feature 2.")
      }
    }
  }

  if (!is.null(label2)) {
    if (!(all(label2 >= 1 & label2 <= N))) {
      stop("The indexes in 'label2' fall out of the range of the total number of observations.")
    }
  }

  for (i in 1:N) {
    if (!is.null(label1) & (i %in% label1)) {
      next
    }

    if (!is.null(label2) & (i %in% label2)) {
      next
    }

    pi_old <- pi[i, ]

    if(pi_old[1] < 0.0001) {
      pi_old[1] <- 0.01
      pi_old[2] <- 1 - pi_old[1]
    }

    if (pi_old[2] < 0.0001) {
      pi_old[2] <- 0.01
      pi_old[1] <- 1 - pi_old[2]
    }

    pi_new <- LaplacesDemon::rdirichlet(n = 1, alpha = tuning_param * pi_old)

    if(pi_new[1] < 0.0001) {
      pi_new[1] <- 0.01
      pi_new[2] <- 1 - pi_new[1]
    }

    if (pi_new[2] < 0.0001) {
      pi_new[2] <- 0.01
      pi_new[1] <- 1 - pi_new[2]
    }

    modelMean_old <- meanWarp_alt(t, c[i], phi[i, ], tt_basis, gamma1, gamma2, pi_old,
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp_alt(t, c[i], phi[i, ], tt_basis, gamma1, gamma2, pi_new,
                              knots_shape, degree, intercept)
    y_i <- y[i, ]

    P1 <- -1 / (2 * var_e) * sum((y_i - modelMean_new) ^ 2) +
      (alpha[1] - 1) * log(pi_new[1]) +
      (alpha[2] - 1) * log(pi_new[2])
    Q1 <- LaplacesDemon::ddirichlet(pi_new, alpha = tuning_param * pi_old, log = T)

    P0 <- -1 / (2 * var_e) * sum((y_i - modelMean_old) ^ 2) +
      (alpha[1] - 1) * log(pi_old[1]) +
      (alpha[2] - 1) * log(pi_old[2])
    Q0 <- LaplacesDemon::ddirichlet(pi_old, alpha = tuning_param * pi_new, log = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      pi_old <- pi_new
    }
    pi[i, ] <- pi_old
  }

  if (rescale == T) {
    pi[, 1] <- (pi[, 1] - min(pi[, 1])) / (max(pi[, 1]) - min(pi[, 1]))
    pi[, 2] <- 1 - pi[, 1]
  }
  return(pi)
}

kfeature_piUpdate_Warp <- function(t, y, c, a, phi, rho, tt_basis, gamma, pi, knots_shape, degree = 3,
                          intercept = F, var_e, alpha, tuning_param = 1000, reg = 1, repulsive = T,
                          temperature = 1) {
  n <- length(t)
  N <- length(c)
  K <- ncol(pi)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (nrow(pi) != N) {
    stop("The number of rows in 'pi' must match the length of 'c'.")
  }


  for (i in 1:N) {
    pi_old <- pi[i, ]
    pi_new <- LaplacesDemon::rdirichlet(n = 1, alpha = tuning_param * pi_old)

    pi_new[pi_new < 1e-4] <- 1e-4
    pi_new <- pi_new / sum(pi_new)

    modelMean_old <- kfeature_meanWarp(t, c[i], a, phi[i, ], rho, tt_basis, gamma, pi_old,
                              knots_shape, degree, intercept)
    modelMean_new <- kfeature_meanWarp(t, c[i], a, phi[i, ], rho, tt_basis, gamma, pi_new,
                              knots_shape, degree, intercept)
    y_i <- y[i, ]

    P1 <- -1 / (2 * var_e) * sum((y_i - modelMean_new) ^ 2)
    P1 <- P1 / temperature

    for (k in 1:K) {
      P1 <- P1 + (alpha[k] - 1) * log(pi_new[k])
    }

    Q1 <- LaplacesDemon::ddirichlet(pi_new, alpha = tuning_param * pi_old, log = T)

    P0 <- -1 / (2 * var_e) * sum((y_i - modelMean_old) ^ 2)
    P0 <- P0 / temperature

    for (k in 1:K) {
      P0 <- P0 + (alpha[k] - 1) * log(pi_old[k])
    }
    if (repulsive == T){
      for (j in 1:N) {
        if(j == i ){
          next
        } else {
          P1 <- P1 - (reg / N) / sum((pi_new - pi[j, ])^2)
          P0 <- P0 - (reg / N) / sum((pi_old - pi[j, ])^2)
        }
      }
    }

    Q0 <- LaplacesDemon::ddirichlet(pi_old, alpha = tuning_param * pi_new, log = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      pi_old <- pi_new
    }
    pi[i, ] <- pi_old
  }

  return(pi)
}

