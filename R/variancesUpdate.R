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
