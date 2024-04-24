meanNoWarp <- function(t, c, gamma1, gamma2, pi, knots_shape, degree, intercept = F) {
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)

  if (intercept == F) {
    df <- p + degree
  } else {
    df <- p + degree + 1
  }

  if (length(gamma1) != df | length(gamma2) != df){
    stop("The length of 'gamma1' and 'gamma2' needs to match the B-spline dimensions.")
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (nrow(pi) != N) {
    stop("The number of rows of 'pi' must match the length of 'c'.")
  }

  spline_basis <- splines::bs(x = t, knots = knots_shape, degree = degree, intercept = intercept)
  f1 <- spline_basis %*% gamma1
  f2 <- spline_basis %*% gamma2

  modelMean <- c + outer(pi[, 1], gamma1) + outer(pi[, 2], gamma2)
  return(modelMean)
}

meanWarp <- function(t, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept = F){
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)
  Q <- ncol(phis)

  if (intercept == F) {
    df <- p + degree
  } else {
    df <- p + degree + 1
  }

  if (length(gamma1) != df | length(gamma2) != df) {
    stop("The length of 'gamma1' and 'gamma2' needs to match the B-spline dimensions.")
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (nrow(pi) != N) {
    stop("The number of rows of 'pi' must match the length of 'c'.")
  }

  if (nrow(phis) != N) {
    stop("The number of rows of 'phi' must match the length of 'c'.")
  }

  if (length(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same length as the
         number of columns in 'phi'.")
  }

  if (rho < 0 | rho > 1) {
    stop("'rho' must be between 0 and 1.")
  }

  modelMean <- matrix(data = NA, nrow = N, ncol = n)
  for (i in 1:N) {
      tWarp2 <- tt_basis %*% phi[j, ]
      tWarp1 <- rho * (tWarp1 - t) + t

      splines_basis1 <- splines::bs(x = tWarp1, knots = knots_shape, degree = degree, intercept = intercept)
      splines_basis2 <- splines::bs(x = tWarp2, knots = knots_shape, degree = degree, intercept = intercept)

      f1 <- spline_basis1 %*% gamma1
      f2 <- spline_basis2 %*% gamma2

      modelMean[i, ] <- c[i] + pi[i, 1] * f1 + pi[i, 2] * f2
  }
  return(modelMean)
}
