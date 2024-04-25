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

  if (N > 1) {
    modelMean <- c + outer(pi[, 1], gamma1) + outer(pi[, 2], gamma2)
  } else {
    modelMean <- c + pi[1] * gamma1 + pi[2] * gamma2
  }

  return(modelMean)
}

meanWarp <- function(t, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept = F){
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)
  Q <- ncol(phi)

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

  if (nrow(phi) != N) {
    stop("The number of rows of 'phi' must match the length of 'c'.")
  }

  if (length(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same length as the
         number of columns in 'phi'.")
  }

  if (rho < 0 | rho > 1) {
    stop("'rho' must be between 0 and 1.")
  }

  if (N == 1) {
    tWarp2 <- tt_basis %*% phi
    tWarp1 <- rho * (tWarp1 - t) + t

    spline_basis1 <- splines::bs(x = tWarp1, knots = knots_shape, degree = degree, intercept = intercept)
    spline_basis2 <- splines::bs(x = tWarp2, knots = knots_shape, degree = degree, intercept = intercept)

    f1 <- spline_basis1 %*% gamma1
    f2 <- spline_basis2 %*% gamma2

    modelMean <- c + pi[1] * f1 + pi[2] * f2

  } else {
    modelMean <- matrix(data = NA, nrow = N, ncol = n)
    for (i in 1:N) {
      tWarp2 <- tt_basis %*% phi[i, ]
      tWarp1 <- rho * (tWarp1 - t) + t

      spline_basis1 <- splines::bs(x = tWarp1, knots = knots_shape, degree = degree, intercept = intercept)
      spline_basis2 <- splines::bs(x = tWarp2, knots = knots_shape, degree = degree, intercept = intercept)

      f1 <- spline_basis1 %*% gamma1
      f2 <- spline_basis2 %*% gamma2

      modelMean[i, ] <- c[i] + pi[i, 1] * f1 + pi[i, 2] * f2
    }
  }
  return(modelMean)
}

meanWarp_alt <- function(t, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept = F){
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)
  Q <- ncol(phi)

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

  if (nrow(phi) != N) {
    stop("The number of rows of 'phi' must match the length of 'c'.")
  }

  if (length(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same length as the
         number of columns in 'phi'.")
  }

  if (N == 1){
    tWarp <- tt_basis %*% phi
    spline_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
    f1 <- spline_basis %*% gamma1
    f2 <- spline_basis %*% gamma2

    modelMean <- c + pi[1] * f1 + pi[2] * f2
  } else {
    modelMean <- matrix(data = NA, nrow = N, ncol = n)
    for (i in 1:N) {
      tWarp <- tt_basis %*% phi[i, ]
      splines_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
      f1 <- spline_basis %*% gamma1
      f2 <- spline_basis %*% gamma2

      modelMean[i, ] <- c[i] + pi[i, 1] * f1 + pi[i, 2] * f2
    }
  }

  return(modelMean)
}

Likelihood <- function(t, y, c, gamma1, gamma2, pi, knots_shape, degree, var_e, phi = NULL,
                       tt_basis = NULL, rho = NULL, intercept = F, log = F){
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (is.null(phi) & is.null(tt_basis) & is.null(rho)) {
    modelMean <- meanNoWarp(t, c, gamma1, gamma2, pi, knots_shape, degree, intercept)
  } else if (!is.null(phi) & !is.null(tt_basis) & !is.null(rho)){
    modelMean <- meanWarp(t, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept)
  } else if (!is.null(phi) & !is.null(tt_basis) & is.null(rho)){
    modelMean <- meanWarp_alt(t, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape, degree, intercept)
  } else if (!is.null(rho) & (is.null(phi) | is.null(tt_basis))) {
    stop("When providing 'rho' you must also provide 'phi' and 'tt_basis'.")
  } else if ((!is.null(phi) & is.null(tt_basis)) | (is.null(phi) & !is.null(tt_basis))) {
    stop ("'phi' and 'tt_basis' need to be provided together.")
  }

  if (N == 1){
    likelihood <- mvtnorm::dmvnorm(x = y, mean = modelMean, sigma = var_e * diag(n), log = log)
  } else {
    if (log == F){
      likelihood <- 1
      for (i in 1:N) {
        likelihood_i <-mvtnorm::dmvnorm(x = y[i, ], mean = modelMean, sigma = var_e * diag(n), log = log)
        likelihood <- likelihood  * likelihood_i
      }
    } else {
      likelihood <- 0
      for (i in 1:N) {
        likelihood_i <-mvtnorm::dmvnorm(x = y[i, ], mean = modelMean, sigma = var_e * diag(n), log = log)
        likelihood <- likelihood + likelihood_i
      }
    }
  }
  return(likelihood)
}






