phiUpdate_NoReg <- function(t, y, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape,
                            degree = 3, intercept = F, var_e, var_phi, Upsilon,
                            tau, it_num, acceptance_sums) {
  n <- length(t)
  N <- length(c)
  Q <- ncol(phi)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (length(tau) != N) {
    stop("The length of 'tau' must match the length of 'c'.")
  }

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  if (is.unsorted(Upsilon)) {
    stop("'Upsilon' needs to have ordered values.")
  }

  if (Upsilon[1] != t[1] | Upsilon[Q] != t[n]) {
    stop("'Upsilon' needs to satisfy the image constraint.")
  }

  it_num <- it_num + 1

  for (i in 1:N) {
    phi_old <- phi[i, ]
    eta_old <- jupp(phi_old)

    eta_new <- mvtnorm::rmvnorm(n = 1, mean = eta_old, sigma = tau[i] * diag(Q))
    eta_new[1] <- 0
    eta_old[Q] <- 1

    phi_new <- juppinv(eta_new)

    modelMean_old <- meanWarp(t, c[i], phi_old, rho, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp(t, c[i], phi_new, rho, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)

    y_i <- y[i, ]

    P1 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_new, sigma = var_e * diag(n), logged = T) +
          mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)],
                           sigma = var_phi * diag(Q - 2), logged = T)
    Q1 <- mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = eta_old[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    P0 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_old, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q0 <- mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = eta_new[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      phi[i, ] <- phi_new
      acceptance_sums[i] <- acceptance_sums[i] + 1
    }
    tau[i] <- tau[i] * (1 + acceptance_sums[i] / it_num - 0.3) / sqrt(it_num)
  }
  phiT <- t(phi)
  jupp_phiT <- apply(phiT, 2, jupp)
  jupp_means <- apply(jupp_phiT, 1, mean)
  jupp_phiT <- jupp_phiT - jupp_means + jupp(Upsilon)
  phiT <- apply(jupp_phiT, 2, juppinv)
  phi <- t(phiT)

  metropolis <- list("phi" = phi, "acceptance" = acceptance_sums, "tau" = tau)
  return(metropolis)
}

phiUpdate_NoReg_alt <- function(t, y, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape,
                            degree = 3, intercept = F, var_e, var_phi, Upsilon,
                            tau, it_num, acceptance_sums) {
  n <- length(t)
  N <- length(c)
  Q <- ncol(phi)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (length(tau) != N) {
    stop("The length of 'tau' must match the length of 'c'.")
  }

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  if (is.unsorted(Upsilon)) {
    stop("'Upsilon' needs to have ordered values.")
  }

  if (Upsilon[1] != t[1] | Upsilon[Q] != t[n]) {
    stop("'Upsilon' needs to satisfy the image constraint.")
  }

  it_num <- it_num + 1

  for (i in 1:N) {
    phi_old <- phi[i, ]
    eta_old <- jupp(phi_old)

    eta_new <- mvtnorm::rmvnorm(n = 1, mean = eta_old, sigma = tau[i] * diag(Q))
    eta_new[1] <- 0
    eta_old[Q] <- 1

    phi_new <- juppinv(eta_new)

    modelMean_old <- meanWarp_alt(t, c[i], phi_old, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp_alt(t, c[i], phi_new, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)

    y_i <- y[i, ]

    P1 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_new, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q1 <- mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = eta_old[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    P0 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_old, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q0 <- mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = eta_new[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      phi[i, ] <- phi_new
      acceptance_sums[i] <- acceptance_sums[i] + 1
    }
    tau[i] <- tau[i] * (1 + acceptance_sums[i] / it_num - 0.3) / sqrt(it_num)
  }
  phiT <- t(phi)
  jupp_phiT <- apply(phiT, 2, jupp)
  jupp_means <- apply(jupp_phiT, 1, mean)
  jupp_phiT <- jupp_phiT - jupp_means + jupp(Upsilon)
  phiT <- apply(jupp_phiT, 2, juppinv)
  phi <- t(phiT)

  metropolis <- list("phi" = phi, "acceptance" = acceptance_sums, "tau" = tau)
  return(metropolis)
}

phiUpdate_Reg <- function(t, y, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape,
                            degree = 3, intercept = F, var_e, var_phi, Upsilon, B, X,
                            tau, it_num, acceptance_sums) {
  n <- length(t)
  N <- length(c)
  Q <- ncol(phi)
  l <- ncol(X)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (length(tau) != N) {
    stop("The length of 'tau' must match the length of 'c'.")
  }

  if (nrow(X) != N) {
    stop("The number of rows in 'X' must match the length of 'c'.")
  }

  if (nrow(B) != l) {
    stop("The number of rows in 'B' must match the number of columns in 'X'.")
  }

  if (ncol(B) != Q - 2) {
    stop("The number of columns in 'B' must be 'ncol(phi) - 2'.")
  }

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  if (is.unsorted(Upsilon)) {
    stop("'Upsilon' needs to have ordered values.")
  }

  if (Upsilon[1] != t[1] | Upsilon[Q] != t[n]) {
    stop("'Upsilon' needs to satisfy the image constraint.")
  }

  it_num <- it_num + 1

  for (i in 1:N) {
    phi_old <- phi[i, ]
    eta_old <- jupp(phi_old)

    eta_new <- mvtnorm::rmvnorm(n = 1, mean = eta_old, sigma = tau[i] * diag(Q))
    eta_new[1] <- 0
    eta_old[Q] <- 1

    phi_new <- juppinv(eta_new)

    modelMean_old <- meanWarp(t, c[i], phi_old, rho, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)
    modelMean_new <- meanWarp(t, c[i], phi_new, rho, tt_basis, gamma1, gamma2, pi[i, ],
                              knots_shape, degree, intercept)

    y_i <- y[i, ]

    P1 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_new, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)] + t(B) %*% X[i, ],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q1 <- mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = eta_old[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    P0 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_old, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)] + t(B) %*% X[i, ],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q0 <- mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = eta_new[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      phi[i, ] <- phi_new
      acceptance_sums[i] <- acceptance_sums[i] + 1
    }
    tau[i] <- tau[i] * (1 + acceptance_sums[i] / it_num - 0.3) / sqrt(it_num)
  }
  phiT <- t(phi)
  jupp_phiT <- apply(phiT, 2, jupp)
  jupp_means <- apply(jupp_phiT, 1, mean)
  jupp_phiT <- jupp_phiT - jupp_means + jupp(Upsilon)
  phiT <- apply(jupp_phiT, 2, juppinv)
  phi <- t(phiT)

  metropolis <- list("phi" = phi, "acceptance" = acceptance_sums, "tau" = tau)
  return(metropolis)
}

phiUpdate_NoReg_alt <- function(t, y, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape,
                                degree = 3, intercept = F, var_e, var_phi, Upsilon, B, X,
                                tau, it_num, acceptance_sums) {
  n <- length(t)
  N <- length(c)
  Q <- ncol(phi)
  l <- ncol(X)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (length(tau) != N) {
    stop("The length of 'tau' must match the length of 'c'.")
  }

  if (nrow(X) != N) {
    stop("The number of rows in 'X' must match the length of 'c'.")
  }

  if (nrow(B) != l) {
    stop("The number of rows in 'B' must match the number of columns in 'X'.")
  }

  if (ncol(B) != Q - 2) {
    stop("The number of columns in 'B' must be 'ncol(phi) - 2'.")
  }

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  if (is.unsorted(Upsilon)) {
    stop("'Upsilon' needs to have ordered values.")
  }

  if (Upsilon[1] != t[1] | Upsilon[Q] != t[n]) {
    stop("'Upsilon' needs to satisfy the image constraint.")
  }

  it_num <- it_num + 1

  for (i in 1:N) {
    phi_old <- phi[i, ]
    eta_old <- jupp(phi_old)

    eta_new <- mvtnorm::rmvnorm(n = 1, mean = eta_old, sigma = tau[i] * diag(Q))
    eta_new[1] <- 0
    eta_old[Q] <- 1

    phi_new <- juppinv(eta_new)

    modelMean_old <- meanWarp_alt(t, c[i], phi_old, tt_basis, gamma1, gamma2, pi[i, ],
                                  knots_shape, degree, intercept)
    modelMean_new <- meanWarp_alt(t, c[i], phi_new, tt_basis, gamma1, gamma2, pi[i, ],
                                  knots_shape, degree, intercept)

    y_i <- y[i, ]

    P1 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_new, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)] + t(B) %*% X[i, ],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q1 <- mvtnorm::dmvnorm(x = eta_new[-c(1, Q)], mu = eta_old[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    P0 <- mvtnorm::dmvnorm(x = y_i, mu = modelMean_old, sigma = var_e * diag(n), logged = T) +
      mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = jupp(Upsilon)[-c(1, Q)] + t(B) %*% X[i, ],
                       sigma = var_phi * diag(Q - 2), logged = T)
    Q0 <- mvtnorm::dmvnorm(x = eta_old[-c(1, Q)], mu = eta_new[-c(1, Q)],
                           sigma = tau[i] * diag(Q - 2), logged = T)

    ratio <- (P1 - Q1) - (P0 - Q0)

    if (log(stats::runif(1)) < ratio) {
      phi[i, ] <- phi_new
      acceptance_sums[i] <- acceptance_sums[i] + 1
    }
    tau[i] <- tau[i] * (1 + acceptance_sums[i] / it_num - 0.3) / sqrt(it_num)
  }
  phiT <- t(phi)
  jupp_phiT <- apply(phiT, 2, jupp)
  jupp_means <- apply(jupp_phiT, 1, mean)
  jupp_phiT <- jupp_phiT - jupp_means + jupp(Upsilon)
  phiT <- apply(jupp_phiT, 2, juppinv)
  phi <- t(phiT)

  metropolis <- list("phi" = phi, "acceptance" = acceptance_sums, "tau" = tau)
  return(metropolis)
}
