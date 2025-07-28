cUpdate_NoWarp <- function(t, y, gamma1, gamma2, pi, shape_basis, var_c, var_e) {
  n <- length(t)
  N <- nrow(y)

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  c <- rep(0, N)
  modelMean <- meanNoWarp(t, c, gamma1, gamma2, pi, shape_basis)

  mean_c <- var_c / (var_e / n + var_c) * rowSums(y - modelMean) / n
  var_c <- 1 / (1 / var_c + n / var_e)

  c <- stats::rnorm(n = N, mean = mean_c, sd = sqrt(var_c))
  c <- c - mean(c)

  return(c)
}

cUpdate_Warp <- function(t, y, phi, rho, tt_basis, gamma1, gamma2, pi,
                         knots_shape, var_c, var_e, degree = 3, intercept = F) {
  n <- length(t)
  N <- nrow(y)

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  c <- rep(0, N)
  modelMean <- meanWarp(t, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept)

  mean_c <- var_c / (var_e / n + var_c) * rowSums(y - modelMean) / n
  var_c <- 1 / (1 / var_c + n / var_e)

  c <- stats::rnorm(n = N, mean = mean_c, sd = sqrt(var_c))
  c <- c - mean(c)

  return(c)
}

cUpdate_Warp_alt <- function(t, y, phi, tt_basis, gamma1, gamma2, pi,
                         knots_shape, var_c, var_e, degree = 3, intercept = F) {
  n <- length(t)
  N <- nrow(y)

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  c <- rep(0, N)
  modelMean <- meanWarp_alt(t, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept)

  mean_c <- var_c / (var_e / n + var_c) * rowSums(y - modelMean) / n
  var_c <- 1 / (1 / var_c + n / var_e)

  c <- stats::rnorm(n = N, mean = mean_c, sd = sqrt(var_c))
  c <- c - mean(c)

  return(c)
}

kfeature_cUpdate_Warp <- function(t, y, phi, rho, tt_basis, gamma, pi,
                         knots_shape, var_c, var_e, degree = 3, intercept = F) {
  n <- length(t)
  N <- nrow(y)

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  c <- rep(0, N)
  modelMean <- meanWarp(t, c, phi, rho, tt_basis, gamma, pi, knots_shape, degree, intercept)

  mean_c <- var_c / (var_e / n + var_c) * rowSums(y - modelMean) / n
  var_c <- 1 / (1 / var_c + n / var_e)

  c <- stats::rnorm(n = N, mean = mean_c, sd = sqrt(var_c))
  c <- c - mean(c)

  return(c)
}
