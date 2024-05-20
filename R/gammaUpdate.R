gammaUpdate_NoWarp <- function(t, y, c, pi, shape_basis, Omega, lambda1, lambda2, var_e) {
  n <- length(t)
  N <- length(c)
  df <- ncol(shape_basis)

  if (nrow(Omega) != df) {
    stop("The dimension of 'Omega' needs to match the number of columns of 'shape_basis'.")
  }

  if (nrow(pi) != N) {
    stop("The number of rows of 'pi' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  preMean <- matrix(rep(0, 2 * df), nrow = 2 * df)
  preCov <- matrix(rep(0, 4 * df ^ 2), nrow = 2 * df)

  for (i in 1:N) {
    S_i <- cbind(pi[i, 1] * shape_basis, pi[i, 2] * shape_basis)
    y_i <- matrix(y[i, ] - c[i], nrow = n)
    preCov <- preCov + t(S_i) %*% S_i
    preMean <- preMean + t(S_i) %*% y_i
  }

  priorPrecision <- matrix(rep(0, 4 * df * df), nrow = 2 * df)
  priorPrecision[1:df, 1:df] <- Omega / lambda1
  priorPrecision[(df + 1):(2 * df), (df + 1):(2 * df)] <- Omega / lambda2

  cov_gamma <- MASS::ginv(priorPrecision + preCov / var_e)
  mean_gamma <- cov_gamma %*% preMean / var_e

  gamma <- mvtnorm::rmvnorm(n = 1, mean = mean_gamma, sigma = cov_gamma)
  gamma[1:df] <- gamma[1:df] - mean(gamma[1:df])

  return(list(gamma[1:df], gamma[(df + 1):(2 * df)]))
}

gammaUpdate_Warp <- function(t, y, c, phi, rho, tt_basis, pi, knots_shape, Omega,
                             lambda1, lambda2, var_e, degree = 3, intercept = F) {
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)
  df <- nrow(Omega)
  Q <- ncol(phi)

  if (intercept == F) {
    if (df != p + degree) {
      stop("The dimension of 'Omega' needs to match the B-spline dimensions implied by the length
           of 'knots_shape', and values of'degree' and 'intercept'.")
    }
  } else {
    if (df != p + degree + 1) {
      stop("The dimension of 'Omega' needs to match the B-spline dimensions implied by the length
           of 'knots_shape', and values of'degree' and 'intercept'.")
    }
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (nrow(pi) != N) {
    stop("The number of rows of 'pi' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (nrow(phi) != N) {
    stop("The number of rows of 'phi' must match the length of 'c'.")
  }

  if (ncol(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same number of columns as 'phi'.")
  }

  if (rho < 0 | rho > 1) {
    stop("'rho' must be between 0 and 1.")
  }

  preMean <- matrix(rep(0, 2 * df), nrow = 2 * df)
  preCov <- matrix(rep(0, 4 * df ^ 2), nrow = 2 * df)

  shape_basis <- splines::bs(x = t, knots = knots_shape, degree = degree, intercept = intercept)

  for (i in 1:N) {
    phi_i <- phi[i, ]
    tWarp1 <- tt_basis %*% phi_i
    tWarp2 <- rho * (tWarp1 - t) + t

    shape_basis1 <- splines::bs(x = tWarp1, knots = knots_shape, degree = degree, intercept = intercept)
    shape_basis2 <- splines::bs(x = tWarp2, knots = knots_shape, degree = degree, intercept = intercept)

    S_i <- cbind(pi[i, 1] * shape_basis1, pi[i, 2] * shape_basis2)
    y_i <- matrix(y[i, ] - c[i], nrow = n)
    preCov <- preCov + t(S_i) %*% S_i
    preMean <- preMean + t(S_i) %*% y_i
  }

  priorPrecision <- matrix(rep(0, 4 * df * df), nrow = 2 * df)
  priorPrecision[1:df, 1:df] <- Omega / lambda1
  priorPrecision[(df + 1):(2 * df), (df + 1):(2 * df)] <- Omega / lambda2

  cov_gamma <- MASS::ginv(priorPrecision + preCov / var_e)
  mean_gamma <- cov_gamma %*% preMean / var_e

  gamma <- mvtnorm::rmvnorm(n = 1, mean = mean_gamma, sigma = cov_gamma)
  gamma[1:df] <- gamma[1:df] - mean(gamma[1:df])

  return(list(gamma[1:df], gamma[(df + 1):(2 * df)]))
}

gammaUpdate_Warp_alt <- function(t, y, c, phi, tt_basis, pi, knots_shape, Omega,
                             lambda1, lambda2, var_e, degree = 3, intercept = F) {
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)
  df <- nrow(Omega)
  Q <- ncol(phi)

  if (intercept == F) {
    if (df != p + degree) {
      stop("The dimension of 'Omega' needs to match the B-spline dimensions implied by the length
           of 'knots_shape', and values of'degree' and 'intercept'.")
    }
  } else {
    if (df != p + degree + 1) {
      stop("The dimension of 'Omega' needs to match the B-spline dimensions implied by the length
           of 'knots_shape', and values of'degree' and 'intercept'.")
    }
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (nrow(pi) != N) {
    stop("The number of rows of 'pi' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (nrow(phi) != N) {
    stop("The number of rows of 'phi' must match the length of 'c'.")
  }

  if (ncol(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same number of columns as 'phi'.")
  }

  preMean <- matrix(rep(0, 2 * df), nrow = 2 * df)
  preCov <- matrix(rep(0, 4 * df ^ 2), nrow = 2 * df)

  shape_basis <- splines::bs(x = t, knots = knots_shape, degree = degree, intercept = intercept)

  for (i in 1:N) {
    phi_i <- phi[i, ]
    tWarp <- tt_basis %*% phi_i

    shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)

    S_i <- cbind(pi[i, 1] * shape_basis, pi[i, 2] * shape_basis)
    y_i <- matrix(y[i, ] - c[i], nrow = n)
    preCov <- preCov + t(S_i) %*% S_i
    preMean <- preMean + t(S_i) %*% y_i
  }

  priorPrecision <- matrix(rep(0, 4 * df * df), nrow = 2 * df)
  priorPrecision[1:df, 1:df] <- Omega / lambda1
  priorPrecision[(df + 1):(2 * df), (df + 1):(2 * df)] <- Omega / lambda2

  cov_gamma <- MASS::ginv(priorPrecision + preCov / var_e)
  mean_gamma <- cov_gamma %*% preMean / var_e

  gamma <- mvtnorm::rmvnorm(n = 1, mean = mean_gamma, sigma = cov_gamma)
  gamma[1:df] <- gamma[1:df] - mean(gamma[1:df])

  return(list(gamma[1:df], gamma[(df + 1):(2 * df)]))
}

