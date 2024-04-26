var_cUpdate <- function(c, a_c, b_c) {
  N <- length(c)
  var_c <- 1 / stats::rgamma(n = 1, shape = a_c + N / 2, rate = b_c + sum(c ^ 2) / 2)
  return(var_c)
}

lambdaUpdate <- function(gamma1, gamma2, a_l, b_l, Omega){
  df <- nrow(Omega)
  gamma1 <- matrix(gamma1, nrow = df)
  gamma2 <- matrix(gamma2, nrow = df)
  lambda1 <- 1 / stats::rgamma(n = 1, shape = a_l + df / 2,
                               rate = b_l + t(gamma1) %*% Omega %*% gamma1 / 2)
  lambda2 <- 1 / stats::rgamma(n = 1, shape = a_l + df / 2,
                               rate = b_l + t(gamma2) %*% Omega %*% gamma2 / 2)
  lambdas <- list(lambda1, lambda2)
  return(lambdas)
}

var_eUpdate_NoWarp <- function(t, y, c, gamma1, gamma2, pi, knots_shape, degree = 3, intercept = F,
                               a_e, b_e) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  modelMean <- meanNoWarp(t, c, gamma1, gamma2, pi, knots_shape, degree, intercept)
  squaredSums <- sum(rowSums((y - modelMean) ^ 2))
  var_e <- 1 / stats::rgamma(n = 1, shape = a_e + n * N / 2, rate = b_e + squaredSums / 2)
  return(var_e)
}

var_eUpdate_Warp <- function(t, y, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape,
                             degree = 3, intercept = F, a_e, b_e) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  modelMean <- meanWarp(t, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept)
  squaredSums <- sum(rowSums((y - modelMean) ^ 2))
  var_e <- 1 / stats::rgamma(n = 1, shape = a_e + n * N / 2, rate = b_e + squaredSums / 2)
  return(var_e)
}

var_eUpdate_Warp_alt <- function(t, y, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape,
                             degree = 3, intercept = F, a_e, b_e) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  modelMean <- meanWarp_alt(t, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept)
  squaredSums <- sum(rowSums((y - modelMean) ^ 2))
  var_e <- 1 / stats::rgamma(n = 1, shape = a_e + n * N / 2, rate = b_e + squaredSums / 2)
  return(var_e)
}

var_phiUpdate_NoReg <- function(phi, Upsilon, a_phi, b_phi) {
  Q <- ncol(phi)
  N <- nrow(phi)

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  jupp_mean <- jupp(Upsilon)[-c(1, Q)]
  var_sum <- 0

  for (i in 1:N) {
    phi_i <- phi[i, ]
    eta_i <- jupp(phi_i)[-c(1, Q)]
    var_sum <- var_sum + t(eta_i - jupp_mean) %*% (eta_i - jupp_mean)
  }
  var_phi <- 1 / stats::rgamma(n = 1, shape = a_phi + N * (Q - 2) / 2, rate = b_phi + var_sum / 2)
  return(var_phi)
}

var_phiUpdate_Reg <- function(phi, Upsilon, a_phi, b_phi, X, B, B0, V0) {
  Q <- ncol(phi)
  N <- nrow(phi)
  l <- ncol(X)

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  if (nrow(X) != N) {
    stop("The number of rows in 'X' must match the number of rows in 'phi'.")
  }

  if (nrow(B) != l) {
    stop("The number of rows in 'B' must match the number of columns in 'X'.")
  }

  if (ncol(B) != Q - 2) {
    stop("The number of columns in 'B' must be 'ncol(phi) - 2'.")
  }

  if (ncol(B) != ncol(B0) | nrow(B) != nrow(B0)) {
    stop("The dimensions of 'B' and 'B0' must match.")
  }

  if (ncol(V0) != l | nrow(V0) != l) {
    stop("'V0' must be a square matrix of dimension matching the number of columns in 'X'.")
  }

  jupp_mean <- jupp(Upsilon)[-c(1, Q)]
  var_sum <- 0

  for (i in 1:N) {
    phi_i <- phi[i, ]
    eta_i <- jupp(phi_i)[-c(1, Q)]
    var_sum <- var_sum + t(eta_i - jupp_mean - t(B) %*% X[i, ]) %*% (eta_i - jupp_mean - t(B) %*% X[i, ])
  }
  posterior_a <- a_phi + (Q - 2) * (N + l) / 2
  B_vec <- matrix(B, nrow = (Q - 2) * l)
  B0_vec <- matrix(B0, nrow = (Q - 2) * l)
  posterior_b <- b_phi + var_sum / 2 +(B_vec - B0_vec) %*% kronecker(diag(Q - 2), V0) %*% (B_vec - B0_vec) / 2
  var_phi <- 1 / stats::rgamma(n = 1, shape = posterior_a, rate = posterior_b)
  return(var_phi)
}
