#' MCMC sampler for curve registration with regression for functional mixed membership models
#'
#' @description
#' `crMM_WarpReg()` runs a Metropolis-within-Gibbs sampler for the analysis of functional data under the two-feature
#' mixed membership model assumption. This function incorporates curve registration to account for phase variation.
#' It differs from [crMM_Warp()] in that it includes covariate information in the estimation of the
#' time-transformation functions through linear regression.
#'
#'
#' @inheritParams crMM_Warp
#' @param X Matrix of covariates for time-transformation regression, where each row corresponds to covariates for
#' one observation.
#' @param B0 Prior mean matrix for the regression coefficient.
#' @param V0 Prior row covariance matrix for the regression coefficient.
#' @param B_init Initial value for the regression coefficient. Default is `B0`.
#'
#' @return
#' A list with the following items:
#'
#' * `gamma1`: a matrix with `num_it` rows of the posterior samples for the B-spline coefficient
#' for the shape associated with the first feature.
#' * `gamma2`: a matrix with `num_it` rows of the posterior samples for the B-spline coefficient
#' for the shape associated with the second feature.
#' * `c`: a matrix with `num_it` rows of the posterior samples for the individual intercepts. The number
#' of columns corresponds to the number of subjects, that is, `nrow(y)`.
#' * `variance`: a matrix with `num_it` rows of the posterior samples for the variance parameters.
#' The five columns, in order, correspond to the variances for the first and second feature B-spline
#' coefficients, the intercepts, the error term, and the time-transformation spline coefficients.
#' * `phi`: a matrix with `num_it` rows of the posterior samples for the B-spline coefficients for the
#' individual time-transformation functions. The coefficients across individuals are stacked in order into
#' one vector corresponding to a row at each iteration.
#' * `pi1`: a matrix with `num_it` rows of the posterior samples of the mixed membership component for the
#' first feature. The number of columns corresponds to the number of subjects, that is,
#' the number of rows in `y`.
#' * `B`: a matrix with `num_it` rows of the posterior samples of the time-transformation regression
#' coefficient matrix. The entries of the matrix are column-stacked into a vector, the columns correspond
#' to the entries.
#' * `fit_sample`: matrix with `num_it` rows of the posterior samples of individual fits. The first
#'  `length(t)` columns give the fit for the first observation, the second `length(t)` columns give
#'   the fit for the second observation, and so on for all the data.
#' * `registered_fit`: matrix with `num_it` rows of the posterior samples of individual aligned fits. The first
#'  `length(t)` columns give the fit for the first observation, the second `length(t)` columns give
#'   the fit for the second observation, and so on for all the data.
#' * `stochastic_time`: matrix with `num_it` rows of the posterior samples of the individual stochastic
#' (aligned) time. The first `length(t)` columns give the times for the first observation, the second
#' `length(t)` columns give the times for the second observation, and so on for all the data.
#' * `fit1`, `fit2`: matrices with `nrow(y)` rows corresponding to the ergodic first and second
#' sample moment of the fit at points `t` for each of the observations.
#'
#' Additionally, if `wantPAF` is `TRUE` then the list also includes a matrix with `num_it` rows of the
#' posterior samples of the Peak Alpha Frequency.
#'
#' Also, if `inc_rho` is `TRUE`, the list includes a matrix with `num_it` rows of the posterior samples
#' for the time-transformation rescaling parameter.
#'
#' @export
#'
crMM_WarpReg <- function(num_it, burnin = 0.2, t, y, X, p, degree_shape = 3, intercept_shape = FALSE,
                      inc_rho = TRUE, rho_init = 0.5, h, degree_tt = 3, intercept_tt = FALSE,
                      a_e, b_e, a_c, b_c, a_l, b_l, a_phi, b_phi, rescale_pi = TRUE, B0, V0, B_init = B0,
                      tuning_pi = 1000, alpha, label1 = NULL, label2 = NULL, wantPAF = FALSE,
                      gamma1_init = NULL, gamma2_init = NULL, lambda1_init = 0.1,
                      lambda2_init = 0.1, var_c_init = 0.1, var_e_init = 1, var_phi_init = 1) {

  # Lengths ---------------------------------------------------------
  n <- length(t)
  N <- nrow(y)
  l <- ncol(X)

  # Stop messages ---------------------------------------------------
  if (is.unsorted(t)) {
    stop("Provide ordered observation times 't'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  # Rescale time points ---------------------------------------------
  t1 <- t[1]
  tn <- t[n]
  t <- (t - t1) / (tn - t1)

  # Pre-calculated quantities ---------------------------------------
  knots_shape <- seq(0, 1, length.out = p + 2)[-c(1, p + 2)]
  shape_basis <- splines::bs(x = t, knots = knots_shape, degree = degree_shape, intercept = intercept_shape)
  Omega <- BsplinePrecision(knots = knots_shape, degree = degree_shape, intercept = intercept_shape)

  if (intercept_tt == T) {
    Q <- h + degree_tt + 1
  } else{
    Q <- h + degree_tt
  }

  knots_tt <- seq(0, 1, length.out = h + 2)[-c(1, h + 2)]
  tt_basis <- splines::bs(x = t, knots = knots_tt, degree = degree_tt, intercept = intercept_tt)

  U <- MASS::ginv(t(X) %*% X + V0)

  # Initialize parameters -------------------------------------------
  lambda1 <- lambda1_init
  lambda2 <- lambda2_init
  var_c <- var_c_init
  var_e <- var_e_init
  var_phi <- var_phi_init

  y_max <- apply(y, 2, max)
  y_min <- apply(y, 2, min)

  if (is.null(gamma1_init)) {
    gamma1 <- MASS::ginv(t(shape_basis) %*% shape_basis + Omega / lambda1) %*% t(shape_basis) %*% y_max
  }

  if (is.null(gamma2_init)) {
    gamma2 <- MASS::ginv(t(shape_basis) %*% shape_basis + Omega / lambda2) %*% t(shape_basis) %*% y_max
  }

  c <- stats::rnorm(n = N, mean = 0, sd = sqrt(var_c))
  c <- c - mean(c)

  pi <- matrix(rep(1/2, 2 * N), nrow = N)
  if (!is.null(label1)) {
    pi[label1, 1] <- 1
    pi[label1, 2] <- 0
  }
  if (!is.null(label2)) {
    pi[label2, 1] <- 0
    pi[label2, 2] <- 1
  }
  pi[, 1] <- (pi[, 1] - min(pi[, 1])) / (max(pi[, 1]) - min(pi[, 1]))
  pi[, 2] <- 1 - pi[, 1]

  Upsilon <- identityTT(Boundary.knots = c(0, 1), knots = knots_tt,
                        degree = degree_tt, intercept = intercept_tt)

  Upsilon[1] <- t[1]
  Upsilon[Q] <- t[n]

  phi <- matrix(rep(Upsilon, N), nrow = N, byrow = T)
  tau <- rep(0.05, N)
  acceptance_sums <- rep(0, N)

  rho <- rho_init

  B <- B_init


  # Storage matrices ------------------------------------------------
  burn_it      <- round(num_it * burnin)
  total_it     <- burn_it + num_it
  gamma1_mat   <- matrix(nrow = num_it, ncol = p + 4, data = NA)
  gamma2_mat   <- matrix(nrow = num_it, ncol = p + 4, data = NA)
  c_mat        <- matrix(nrow = num_it, ncol = N, data = NA)
  var_mat      <- matrix(nrow = num_it, ncol = 5, data = NA)
  pi1_mat      <- matrix(nrow = num_it, ncol = N, data = NA)
  phi_mat      <- matrix(nrow = num_it, ncol = N * Q, data = NA)
  B_mat        <- matrix(nrow = num_it, ncol = l * (Q - 2), data = NA)
  fit_mat      <- matrix(nrow = num_it, ncol = N * n, data = NA)
  register_mat <- matrix(nrow = num_it, ncol = N * n, data = NA)
  tt_mat       <- matrix(nrow = num_it, ncol = N * n, data = NA)
  current_fit  <- matrix(0, nrow = N, ncol = n)
  register_fit <- matrix(0, nrow = N, ncol = n)
  fit          <- matrix(0, nrow = N, ncol = n)
  fit2         <- matrix(0, nrow = N, ncol = n)

  if (wantPAF == T) {
    paf_mat <- matrix(nrow = num_it, ncol = 1, data = NA)
  }

  if (inc_rho == T) {
    rho_mat <- matrix(nrow = num_it, ncol = 1, data = NA)
  }

  # Sampling --------------------------------------------------------
  for (i in 1:total_it) {
    if (inc_rho == T) {
      c <- cUpdate_Warp(t = t, y = y, phi = phi, rho = rho, tt_basis = tt_basis,
                        gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                        var_c = var_c, var_e = var_e, degree = degree_shape, intercept = intercept_shape)
    } else {
      c <- cUpdate_Warp_alt(t = t, y = y, phi = phi, tt_basis = tt_basis,
                            gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                            var_c = var_c, var_e = var_e, degree = degree_shape, intercept = intercept_shape)
    }

    var_c <- var_cUpdate(c = c, a_c = a_c, b_c = b_c)

    if (inc_rho == T) {
      gamma <- gammaUpdate_Warp(t = t, y = y, c = c, phi = phi, rho = rho, tt_basis = tt_basis, pi = pi,
                                knots_shape = knots_shape, Omega = Omega, lambda1 = lambda1, lambda2 = lambda2,
                                var_e = var_e, degree = degree_shape, intercept = intercept_shape)
    } else {
      gamma <- gammaUpdate_Warp_alt(t = t, y = y, c = c, phi = phi, tt_basis = tt_basis, pi = pi,
                                    knots_shape = knots_shape, Omega = Omega,
                                    lambda1 = lambda1, lambda2 = lambda2, var_e = var_e,
                                    degree = degree_shape, intercept = intercept_shape)
    }

    gamma1 <- gamma[[1]]
    gamma2 <- gamma[[2]]

    lambda <- lambdaUpdate(gamma1 = gamma1, gamma2 = gamma2, a_l = a_l, b_l = b_l, Omega = Omega)
    lambda1 <- lambda[[1]]
    lambda2 <- lambda[[2]]

    if (inc_rho == T) {
      pi <- piUpdate_Warp(t = t, y = y, c = c, phi = phi, rho = rho, tt_basis = tt_basis,
                          gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                          degree = degree_shape, intercept = intercept_shape, var_e = var_e,
                          alpha = alpha, rescale = rescale_pi, label1 = label1, label2 = label2,
                          tuning_param = tuning_pi)
    } else {
      pi <- piUpdate_Warp_alt(t = t, y = y, c = c, phi = phi, tt_basis = tt_basis,
                              gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                              degree = degree_shape, intercept = intercept_shape, var_e = var_e,
                              alpha = alpha, rescale = rescale_pi, label1 = label1, label2 = label2,
                              tuning_param = tuning_pi)
    }

    if (inc_rho == T) {
      var_e <- var_eUpdate_Warp(t = t, y = y, c = c, phi = phi, rho = rho, tt_basis = tt_basis,
                                gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                                degree = degree_shape, intercept = intercept_shape, a_e = a_e, b_e = b_e)
    } else {
      var_e <- var_eUpdate_Warp_alt(t = t, y = y, c = c, phi = phi, tt_basis = tt_basis,
                                    gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                                    degree = degree_shape, intercept = intercept_shape, a_e = a_e, b_e = b_e)
    }
    if (inc_rho == T) {
      rho <- rhoUpdate(t = t, y = y, c = c, phi = phi, rho = rho, tt_basis = tt_basis, pi = pi,
                       gamma1 = gamma1, gamma2 = gamma2, knots_shape = knots_shape,
                       degree = degree_shape, intercept = intercept_shape)
    }

    if (inc_rho == T) {
      phi_out <- phiUpdate_Reg(t = t, y = y, c = c, phi = phi, rho = rho, tt_basis = tt_basis,
                                 gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                                 degree = degree_shape, intercept = intercept_shape, var_e = var_e,
                                 var_phi = var_phi, Upsilon = Upsilon, B = B, X = X, tau = tau, it_num = i,
                                 acceptance_sums = acceptance_sums)
    } else {
      phi_out <-  phiUpdate_Reg_alt(t = t, y = y, c = c, phi = phi, tt_basis = tt_basis,
                                      gamma1 = gamma1, gamma2 = gamma2, pi = pi, knots_shape = knots_shape,
                                      degree = degree_shape, intercept = intercept_shape, var_e = var_e,
                                      var_phi = var_phi, Upsilon = Upsilon, tau = tau, B = B, X = X,
                                      it_num = i, acceptance_sums = acceptance_sums)
    }
    phi <- phi_out$phi
    acceptance_sums <- phi_out$acceptance
    tau <- phi_out$tau

    B <- BetaUpdate(phi = phi, X = X, Upsilon = Upsilon, B0 = B0, V0 = V0, var_phi = var_phi, U = U)

    var_phi <- var_phiUpdate_Reg(phi = phi, Upsilon = Upsilon, a_phi = a_phi, b_phi = b_phi,
                                 X = X, B = B, B0 = B0, V0 = V0)

    if (wantPAF == T) {
      eval_grid <- seq(t1, tn, length = 5000)
      spline_eval_max <- splines::bs(x = eval_grid, knots = knots_shape,
                                     degree = degree_shape, intercept = intercept_shape)
      shape2_eval <- gamma2 %*% t(spline_eval_max)
      peak_location <- eval_grid[(order(shape2_eval, decreasing = T)[1])]
    }

    # Storing the samples
    if (i > burn_it) {
      indexing <- i - burn_it
      gamma1_mat[indexing, ] <- gamma1
      gamma2_mat[indexing, ] <- gamma2
      c_mat[indexing, ]      <- c
      pi1_mat[indexing, ]    <- pi[, 1]
      var_mat[indexing, ]    <- c(lambda1, lambda2, var_c, var_e, var_phi)
      phi_mat[indexing, ]    <- matrix(phi, ncol = N * Q)
      B_mat[indexing, ]      <- matrix(B, nrow = l * (Q - 2))

      if (inc_rho == T) {
        rho_mat[indexing, ] <- rho
      }

      if (wantPAF == T) {
        paf_mat[indexing, ] <- peak_location
      }

      # Ergodic mean and current fit
      for (j in 1:N) {
        phi_j <- phi[j, ]
        tWarp1 <- tt_basis %*% phi_j
        shape_basis1 <- splines::bs(x = tWarp1, knots = knots_shape,
                                    degree = degree_shape, intercept = intercept_shape)

        if (inc_rho == T) {
          tWarp2 <- rho * (tWarp1 - t) + t
          shape_basis2 <- splines::bs(x = tWarp2, knots = knots_shape,
                                      degree = degree_shape, intercept = intercept_shape)
        }

        if (inc_rho == T) {
          current_fit[j, ] <- c[j] + pi[j, 1] * shape_basis1 %*% gamma1 + pi[j, 2] * shape_basis2 %*% gamma2
        } else {
          current_fit[j, ] <- c[j] + pi[j, 1] * shape_basis2 %*% gamma1 + pi[j, 2] * shape_basis2 %*% gamma2
        }

        register_fit[j, ] <- c[j] + pi[j, 1] * shape_basis %*% gamma1 + pi[j, 2] * shape_basis %*% gamma2

        fit_mat[indexing, ((j - 1) * n + 1): (j * n)] <- current_fit[j, ]
        register_mat[indexing, ((j - 1) * n + 1): (j * n)] <- register_fit[j, ]
        tt_mat[indexing, ((j - 1) * n + 1): (j * n)] <- tWarp1
      }
      r1   <- (indexing - 1.0) / indexing
      r2   <- 1 / indexing
      fit  <- r1 * fit + r2 * current_fit
      fit2 <- r1 * fit2 + r2* current_fit ^ 2
    }
  }
  final <- list("gamma1" = gamma1_mat,
                "gamma2" = gamma2_mat,
                "c" = c_mat,
                "variance" = var_mat,
                "phi" = phi_mat,
                "pi1" = pi1_mat,
                "B" = B_mat,
                "fit_sample" = fit_mat,
                "registered_fit" = register_mat,
                "stochastic_time" = tt_mat,
                "fit" = fit,
                "fit2" = fit2)

  if(inc_rho == T) {
    final$rho <- rho_mat
  }

  if (wantPAF == T) {
    final$PAF <- paf_mat
  }

  return(final)
}




