meanNoWarp <- function(t, c, gamma1, gamma2, pi, shape_basis) {
  n <- length(t)
  N <- length(c)
  df <- ncol(shape_basis)

  if (length(gamma1) != df | length(gamma2) != df){
    stop("The length of 'gamma1' and 'gamma2' needs to match the B-spline dimensions
         implied by the length of 'knots_shape', and values of 'degree' and 'intercept'.")
  }

  if (is.matrix(pi)) {
    if (nrow(pi) != N) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }
  } else {
    if (N != 1) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }
  }

  f1 <- as.numeric(shape_basis %*% gamma1)
  f2 <- as.numeric(shape_basis %*% gamma2)

  if (N > 1) {
    modelMean <- c + outer(pi[, 1], f1) + outer(pi[, 2], f2)
  } else {
    modelMean <- c + pi[1] * f1 + pi[2] * f2
  }

  return(modelMean)
}

meanWarp <- function(t, c, phi, rho, tt_basis, gamma1, gamma2, pi, knots_shape, degree = 3, intercept = F) {
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)

  if (is.matrix(phi)){
    Q <- ncol(phi)
  } else {
    Q <- length(phi)
  }

  if (intercept == F) {
    df <- p + degree
  } else {
    df <- p + degree + 1
  }

  if (length(gamma1) != df | length(gamma2) != df) {
    stop("The length of 'gamma1' and 'gamma2' needs to match the B-spline dimensions
         implied by the length of 'knots_shape', and values of 'degree' and 'intercept'.")
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (is.matrix(pi)) {
    if (nrow(pi) != N) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }
  } else {
    if (N != 1) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }
  }

  if (is.matrix(phi)) {
    if (nrow(phi) != N) {
      stop("The number of rows of 'phi' must match the length of 'c'.")
    }
  } else {
    if (N != 1) {
      stop("The number of rows of 'phi' must match the length of 'c'.")
    }
  }

  if (ncol(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same number of columns as 'phi'.")
  }

  # if (rho < 0 | rho > 1) {
  #   stop("'rho' must be between 0 and 1.")
  # }

  if (N == 1) {
    tWarp1 <- tt_basis %*% phi
    tWarp2 <- rho * (tWarp1 - t) + t

    shape_basis1 <- splines::bs(x = tWarp1, knots = knots_shape, degree = degree, intercept = intercept)
    shape_basis2 <- splines::bs(x = tWarp2, knots = knots_shape, degree = degree, intercept = intercept)

    f1 <- shape_basis1 %*% gamma1
    f2 <- shape_basis2 %*% gamma2

    modelMean <- c + pi[1] * f1 + pi[2] * f2

  } else {
    modelMean <- matrix(data = NA, nrow = N, ncol = n)
    for (i in 1:N) {
      tWarp1 <- tt_basis %*% phi[i, ]
      tWarp2 <- rho * (tWarp1 - t) + t

      shape_basis1 <- splines::bs(x = tWarp1, knots = knots_shape, degree = degree, intercept = intercept)
      shape_basis2 <- splines::bs(x = tWarp2, knots = knots_shape, degree = degree, intercept = intercept)

      f1 <- shape_basis1 %*% gamma1
      f2 <- shape_basis2 %*% gamma2

      modelMean[i, ] <- c[i] + pi[i, 1] * f1 + pi[i, 2] * f2
    }
  }
  return(modelMean)
}

meanWarp_alt <- function(t, c, phi, tt_basis, gamma1, gamma2, pi, knots_shape, degree = 3, intercept = F) {
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)

  if (is.matrix(phi)){
    Q <- ncol(phi)
  } else {
    Q <- length(phi)
  }

  if (intercept == F) {
    df <- p + degree
  } else {
    df <- p + degree + 1
  }

  if (length(gamma1) != df | length(gamma2) != df) {
    stop("The length of 'gamma1' and 'gamma2' needs to match the B-spline dimensions
         implied by the length of 'knots_shape', and values of 'degree' and 'intercept'.")
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (is.matrix(pi)) {
    if (nrow(pi) != N) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }
  } else {
    if (N != 1) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }
  }

  if (is.matrix(phi)) {
    if (nrow(phi) != N) {
      stop("The number of rows of 'phi' must match the length of 'c'.")
    }
  } else {
    if (N != 1) {
      stop("The number of rows of 'phi' must match the length of 'c'.")
    }
  }

  if (ncol(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same number of columns as 'phi'.")
  }

  if (N == 1){
    tWarp <- tt_basis %*% phi
    shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
    f1 <- shape_basis %*% gamma1
    f2 <- shape_basis %*% gamma2

    modelMean <- c + pi[1] * f1 + pi[2] * f2
  } else {
    modelMean <- matrix(data = NA, nrow = N, ncol = n)
    for (i in 1:N) {
      tWarp <- tt_basis %*% phi[i, ]
      shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
      f1 <- shape_basis %*% gamma1
      f2 <- shape_basis %*% gamma2

      modelMean[i, ] <- c[i] + pi[i, 1] * f1 + pi[i, 2] * f2
    }
  }

  return(modelMean)
}

kfeature_meanWarp <- function(t, c, a, phi, rho, tt_basis, gamma, pi, knots_shape, degree = 3, intercept = F) {
  n <- length(t)
  N <- length(c)
  p <- length(knots_shape)

  if (is.matrix(gamma)) {
    K <- ncol(gamma)
  } else {
    K <- 1
  }

  if (!is.null(rho)) {
    warp_num <- length(rho)
  }

  if (is.matrix(phi)){
    Q <- ncol(phi)
  } else {
    Q <- length(phi)
  }

  if (intercept == F) {
    df <- p + degree
  } else {
    df <- p + degree + 1
  }

  if (is.matrix(gamma)) {
    if (nrow(gamma) != df) {
      stop("The number of rows of 'gamma' must match the B-spline dimensions
           implied by the length of 'knots_shape', and values of 'degree' and 'intercept'.")
    }
  } else {
    if (K > 1) {
      stop("The number of columns of 'gamma' must match the number of latent features.")
    }
    if (length(gamma) != df) {
      stop("The length of 'gamma' must match the B-spline dimensions
           implied by the length of 'knots_shape', and values of 'degree' and 'intercept'.")
    }
  }

  if (min(knots_shape) <= min(t) | max(knots_shape) >= max(t)) {
    stop("The inner knots in 'knots_shape' should all fall strictly between the smallest
         and largest value of 't'.")
  }

  if (is.matrix(pi)) {
    if (nrow(pi) != N) {
      stop("The number of rows of 'pi' must match the length of 'c'.")
    }

    if (ncol(pi) != K) {
      stop("The number of columns of 'pi' must match the number of latent features.")
    }
  } else {
    if (K > 1) {
      if (N != 1) {
        stop("The number of rows of 'pi' must match the length of 'c'.")
      }

      if (length(pi) != K) {
        stop("The number of columns of 'pi' must match the number of latent features.")
      }
    }
  }

  if (is.matrix(phi)) {
    if (nrow(phi) != N) {
      stop("The number of rows of 'phi' must match the length of 'c'.")
    }
  } else {
    if (N != 1) {
      stop("The number of rows of 'phi' must match the length of 'c'.")
    }
  }

  if (ncol(tt_basis) != Q) {
    stop("The time-transformation spline basis 'tt_basis' must have the same number of columns as 'phi'.")
  }

  if (N == 1) {
    if (K == 1){
      tWarp <- tt_basis %*% phi
      shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
      f <- shape_basis %*% gamma
      modelMean <- c + a * pi * f
    } else {
      modelMean <- c
      for (k in 1:K) {
        if (k == 1){
          tWarp <- tt_basis %*% phi
        } else {
          if (is.null(rho)){
            tWarp <- t
          } else {
            if (k > 1 & warp_num >= (k - 1)) {
              tWarp <- rho[k - 1] * (tt_basis %*% phi - t) + t
            } else {
              tWarp <- t
            }
          }
        }
        shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
        f <- shape_basis %*% gamma[, k]
        modelMean <- modelMean + a[k] * pi[k] * f
      }
    }

  } else {
    modelMean <- matrix(data = NA, nrow = N, ncol = n)
    for (i in 1:N) {
      if (K == 1){
        tWarp <- tt_basis %*% phi[i, ]
        shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
        f <- shape_basis %*% gamma
        modelMean_temp <- c[i] + a * pi[i] * f
      } else {
        modelMean_temp <- c[i]
        for (k in 1:K) {
          if (k == 1){
            tWarp <- tt_basis %*% phi[i, ]
          } else {
            if (is.null(rho)){
              tWarp <- t
            } else {
              if (k > 1 & warp_num >= (k - 1)) {
                tWarp <- rho[k - 1] * (tt_basis %*% phi[i, ] - t) + t
              } else {
                tWarp <- t
              }
            }
          }
          shape_basis <- splines::bs(x = tWarp, knots = knots_shape, degree = degree, intercept = intercept)
          f <- shape_basis %*% gamma[, k]
          modelMean_temp <- modelMean_temp + a[k] * pi[i, k] * f
        }
      }
      modelMean[i, ] <- modelMean_temp
    }
  }
  return(modelMean)
}

Likelihood <- function(t, y, c, gamma1, gamma2, pi, shape_basis = NULL, knots_shape = NULL, degree = 3,
                       var_e, phi = NULL, tt_basis = NULL, rho = NULL, intercept = F, log = F) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (is.null(shape_basis)) {
    if (is.null(phi) & is.null(tt_basis) & is.null(rho)) {
      stop("'shape_basis' needs to be provided in the no warping case, i.e. when 'phi',
           'tt_basis' and 'rho' are not provided")
    } else if (is.null(knots_shape)) {
      stop("'knots_shape' knots_shape needs to be provided for the warping case.")
    }
  }

  if (is.null(phi) & is.null(tt_basis) & is.null(rho)) {
    modelMean <- meanNoWarp(t, c, gamma1, gamma2, pi, shape_basis)
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
        likelihood_i <-mvtnorm::dmvnorm(x = y[i, ], mean = modelMean[i, ], sigma = var_e * diag(n), log = log)
        likelihood <- likelihood  * likelihood_i
      }
    } else {
      likelihood <- 0
      for (i in 1:N) {
        likelihood_i <-mvtnorm::dmvnorm(x = y[i, ], mean = modelMean[i, ], sigma = var_e * diag(n), log = log)
        likelihood <- likelihood + likelihood_i
      }
    }
  }
  return(likelihood)
}

kfeature_Likelihood <- function(t, y, c, a, gamma, pi, knots_shape, degree = 3,
                       var_e, phi, tt_basis, rho, intercept = F, log = F) {
  n <- length(t)
  N <- length(c)

  if (nrow(y) != N){
    stop("The number of rows in 'y' must match the length of 'c'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  modelMean <- kfeature_meanWarp(t, c, a, phi, rho, tt_basis, gamma, pi, knots_shape, degree = degree, intercept = intercept)

  if (N == 1){
    likelihood <- mvtnorm::dmvnorm(x = y, mean = modelMean, sigma = var_e * diag(n), log = log)
  } else {
    if (log == F){
      likelihood <- 1
      for (i in 1:N) {
        likelihood_i <-mvtnorm::dmvnorm(x = y[i, ], mean = modelMean[i, ], sigma = var_e * diag(n), log = log)
        likelihood <- likelihood  * likelihood_i
      }
    } else {
      likelihood <- 0
      for (i in 1:N) {
        likelihood_i <-mvtnorm::dmvnorm(x = y[i, ], mean = modelMean[i, ], sigma = var_e * diag(n), log = log)
        likelihood <- likelihood + likelihood_i
      }
    }
  }
  return(likelihood)
}













