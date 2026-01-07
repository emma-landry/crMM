#' MCMC sampler for curve registration for functional mixed membership models
#'
#' @description
#' `crMM_kfeature()` runs a Metropolis-within-Gibbs sampler for the analysis of functional data under the K-
#' feature mixed membership model assumption. This function incorporates curve registration to account for
#' phase variation.
#'
#' @inheritParams crMM_NoWarp
#' @param K The number of latent features assumed by the model.
#' @param warp_num The model assumes that the first feature is warped. `warp_num` specifies the number
#' of additional features that are warped, hence it can take a value between `0` and `K - 1`. `0` means
#' only the first feature is warped, and `K - 1` implies all are warped. Differential warping is assumed,
#' that is the level of warping of each features relative to the reference one is scaled. In practice, this
#' may be overly restrictive, and setting `warp_num` to `NULL` indicates that all `K` features are warped
#' by the same function (no `rho` parameters are estimated).
#' @param rho_init Initial value for the time-transformation rescaling parameters allowing for differential
#' warping of the `warp_num` extra features. Default is `rep(0.5, warp_num)`.
#' @param h The number of inner knots used to define the B-splines for the time-transformation functions.
#' The sampler assumes equally spaced knots between the first and last value of `t`.
#' @param degree_tt Degree of the piecewise polynomial of B-splines for the time-transformation functions.
#' The default is `3` for cubic splines.
#' @param intercept_tt If `TRUE` an intercept is included in the B-spline basis for the time-transformation
#'functions. The default is `FALSE`.
#' @param a_e,a_c,a_l,a_phi Shape parameters of the Inverse Gamma priors.
#' @param b_e,b_c,b_l,b_phi Rate parameters of the Inverse Gamma priors.
#' @param common_a The `K` shape functions are learned as having normalized amplitude. If `TRUE`, the amplitude
#' parameter is assumed equal across all shapes. If `FALSE`, a different amplitude parameter can be obtained for each feature.
#' Default is `TRUE`.
#' @param repulsive_pi If `TRUE`, a repulsive prior is assumed for the mixed membership coefficients.
#' Default is `TRUE`.
#' @param reg Hyperparameter for the repulsive prior. Default is `1`.
#' @param gamma_init Initial value of the B-spline coefficients, where each column corresponds to the coefficients for one
#' feature. Default is `NULL`, in which case the initial values are determined automatically based on provided data.
#' @param lambda_init Initial values for the prior variance of each of the spline coefficients. Default is
#' `rep(0.1, K)`.
#' @param var_c_init,var_e_init,var_phi_init Initial values for prior variances.
#' Default is `0.1` for `var_c_init`, and `1` for `var_e_init` and `var_phi_init`.
#'
#'  @returns An object of class `crMM_Obj`, consisting of a list with the following items:
#'
#' * `gamma`: a list containing `K` matrices of with `num_it` rows of the posterior samples for the
#' B-spline coefficients associated with each feature.
#' * `c`: a matrix with `num_it` rows of the posterior samples for the individual intercepts. The number
#' of columns corresponds to the number of subjects, that is, `nrow(y)`.
#' * `variance`: a matrix with `num_it` rows of the posterior samples for the variance parameters.
#' The `K + 3` columns, in order, correspond to the variances for the B-spline coefficients,
#' the intercepts, the error term, and the time-transformation spline coefficients.
#' * `phi`: a matrix with `num_it` rows of the posterior samples for the B-spline coefficients for the
#' individual time-transformation functions. The coefficients across individuals are stacked in order into
#' one vector corresponding to a row at each iteration.
#' * `pi`: a list containing `K` matrices with `num_it` rows of the mixed membership component for each
#' of the features. The number of columns corresponds to the number of subjects, that is,
#' the number of rows in `y`.
#' * `fit_sample`: matrix with `num_it` rows of the posterior samples of individual fits. The first
#' `length(t)` columns give the fit for the first observation, the second `length(t)` columns give
#' the fit for the second observation, and so on for all the data.
#' * `registered_fit`: matrix with `num_it` rows of the posterior samples of individual aligned fits. The first
#' `length(t)` columns give the fit for the first observation, the second `length(t)` columns give
#' the fit for the second observation, and so on for all the data.
#' * `stochastic_time`: matrix with `num_it` rows of the posterior samples of the individual stochastic
#' (aligned) time. The first `length(t)` columns give the times for the first observation, the second
#' `length(t)` columns give the times for the second observation, and so on for all the data.
#' * `fit1`, `fit2`: matrices with `nrow(y)` rows corresponding to the ergodic first and second
#' sample moment of the fit at points `t` for each of the observations.
#' * `loglik` : vector with `num_it` rows of the full model loglikelihood at each iteration.
#' * `loglik_i` : matrix with `num_it` rows and `nrow(y)` columns of the individual loglikelihood of each observation.
#'
#'
#' Additionally, if `warp_num` is >0, the list includes a matrix `rho` with `num_it` rows and `warp_num` columns
#' of the posterior samples for the time-transformation rescaling parameter.
#' `stochastic_time_rho` is then a list of `warp_num` matrices with `num_it` rows of the posterior samples
#' of the individual stochastic (aligned) time for the scaled additional features.
#' The first `length(t)` columns give the times for the first observation, the second
#' `length(t)` columns give the times for the second observation, and so on for all the data.
#'
#' @export
#'
crMM_kfeature <- function(num_it, burnin = 0.2, t, y, p, degree_shape = 3, intercept_shape = FALSE,
                      K, warp_num, rho_init = rep(0.5, warp_num), h, degree_tt = 3, intercept_tt = FALSE,
                      a_e, b_e, a_c, b_c, a_l, b_l, a_phi, b_phi, common_a = T, repulsive_pi = T, reg = 1,
                      tuning_pi = 1000, alpha, gamma_init = NULL, lambda_init = rep(0.1, K),
                      var_c_init = 0.1, var_e_init = 1, var_phi_init = 1,
                      process_id = NULL) {

  # Lengths ---------------------------------------------------------
  n <- length(t)
  N <- nrow(y)

  # Stop messages ---------------------------------------------------
  if (is.unsorted(t)) {
    stop("Provide ordered observation times 't'.")
  }

  if (ncol(y) != n) {
    stop("The number columns in 'y' must match the length of 't'.")
  }

  if (!is.null(warp_num) && warp_num > (K - 1)) {
    stop("`warp_num` must be between 0 and K - 1.")
  }

  if (is.null(warp_num)) {
    rho_init <- rep(1, K - 1)
    warp_num <- -1
  } else if (warp_num == 0) {
    rho_init <- NULL
  } else {
    if (length(rho_init) != warp_num) {
      stop(
        sprintf(
          "`rho_init` must have length %d (equal to `warp_num`).",
          warp_num
        )
      )
    }
    if (any(rho_init < 0 | rho_init > 1)) {
      rho_init[rho_init < 0 | rho_init > 1] <- 0.5
      warning("All entries of `rho_init` must be between 0 and 1. Invalid entries have been set to the default, 0.5.")
    }
  }

  # Rescale time points ---------------------------------------------
  t1 <- t[1]
  tn <- t[n]
  t <- (t - t1) / (tn - t1)

  # Pre-calculated quantities ---------------------------------------
  knots_shape <- seq(0, 1, length.out = p + 2)[-c(1, p + 2)]
  shape_basis <- splines::bs(x = t, knots = knots_shape, degree = degree_shape, intercept = intercept_shape)
  Omega <- BsplinePrecision(knots = knots_shape, degree = degree_shape, intercept = intercept_shape)

  if (intercept_shape == T) {
    df<- p + degree_shape + 1
  } else{
    df <- p + degree_shape
  }

  if (intercept_tt == T) {
    Q <- h + degree_tt + 1
  } else{
    Q <- h + degree_tt
  }

  knots_tt <- seq(0, 1, length.out = h + 2)[-c(1, h + 2)]
  tt_basis <- splines::bs(x = t, knots = knots_tt, degree = degree_tt, intercept = intercept_tt)

  # Initialize parameters -------------------------------------------
  lambda <- lambda_init
  var_c <- var_c_init
  var_e <- var_e_init
  var_phi <- var_phi_init

  y_max <- apply(y, 2, max)
  y_min <- apply(y, 2, min)

  if (is.null(gamma_init)){
    if (K == 1){
      gamma <-  MASS::ginv(t(shape_basis) %*% shape_basis + Omega / lambda) %*% t(shape_basis) %*% y_max
    } else {
      gamma <- matrix(nrow = df, ncol = K)
      for (k in 1:K) {
        gamma[, k] <-  MASS::ginv(t(shape_basis) %*% shape_basis + Omega / lambda[k]) %*% t(shape_basis) %*% (colMeans(y) + stats::rnorm(n, 0, 0.1))
      }
    }

  } else {
    gamma <- gamma_init
  }

  c <- stats::rnorm(n = N, mean = 0, sd = sqrt(var_c))
  c <- c - mean(c)

  if (K == 1){
    pi <- rep(1, N)
  } else {
    pi <- matrix(rep(1 / K, K * N), nrow = N)
  }

  Upsilon <- identityTT(Boundary.knots = c(0, 1), knots = knots_tt,
                        degree = degree_tt, intercept = intercept_tt)

  Upsilon[1] <- 0
  Upsilon[Q] <- 1

  phi <- matrix(rep(Upsilon, N), nrow = N, byrow = T)
  tau <- rep(0.05, N)
  acceptance_sums <- rep(0, N)

  rho <- rho_init

  a <- rep(1, K)

  if (0 <= burnin & burnin <= 1 ){
    burn_it <- round(num_it * burnin)
  } else {
    burn_it <- burnin
  }

  total_it     <- burn_it + num_it


  # Storage matrices ------------------------------------------------
  gamma_mats <- vector("list", K)
  pi_mats <- vector("list", K)
  for (k in 1:K) {
    gamma_mats[[k]] <- matrix(nrow = num_it, ncol = df, data = NA)
    pi_mats[[k]] <- matrix(nrow = num_it, ncol = N, data = NA)
  }

  c_mat        <- matrix(nrow = num_it, ncol = N, data = NA)
  a_mat        <- matrix(nrow = num_it, ncol = K, data = NA)
  var_mat      <- matrix(nrow = num_it, ncol = (3 + K), data = NA)
  phi_mat      <- matrix(nrow = num_it, ncol = N * Q, data = NA)
  fit_mat      <- matrix(nrow = num_it, ncol = N * n, data = NA)
  register_mat <- matrix(nrow = num_it, ncol = N * n, data = NA)
  tt_mat       <- matrix(nrow = num_it, ncol = N * n, data = NA)
  loglik       <- matrix(nrow = num_it, ncol = 1, data = NA)
  loglik_i     <- matrix(nrow = num_it, ncol = N, data = NA)
  current_fit  <- matrix(0, nrow = N, ncol = n)
  register_fit <- matrix(0, nrow = N, ncol = n)
  fit          <- matrix(0, nrow = N, ncol = n)
  fit2         <- matrix(0, nrow = N, ncol = n)

  if (warp_num > 0) {
    rho_mat <- matrix(nrow = num_it, ncol = warp_num, data = NA)

    tt_mats <- vector("list", warp_num)
    for (w in 1:warp_num) {
      tt_mats[[w]] <- matrix(nrow = num_it, ncol = N * n, data = NA)
    }
  }

  # Sampling --------------------------------------------------------
  for (i in 1:total_it) {
    if (i %% 100 == 0) print(i)

    if (i %% 2000 == 0 & !is.null(process_id)) {
      system(sprintf('echo "\n%s - Process %s - Completed %d iterations\n"', Sys.time(), process_id, i))
    }

    if (i <= burn_it) {
      temperature_gamma <- max(1, 1 + 5 * (1 - i / burn_it) ^ 0.5)
    } else {
      temperature_gamma <- 1
    }

    gamma <- kfeature_gammaUpdate_Warp(
      t = t, y = y, c = c, a = a, K = K, phi = phi, rho = rho, tt_basis = tt_basis, pi = pi,
      knots_shape = knots_shape, Omega = Omega, lambda = lambda,
      var_e = var_e, degree = degree_shape, intercept = intercept_shape, temperature = temperature_gamma
    )

    lambda <- kfeature_lambdaUpdate(gamma = gamma, a_l = a_l, b_l = b_l, Omega = Omega)

    c <- kfeature_cUpdate_Warp(t = t, y = y, a = a, phi = phi, rho = rho, tt_basis = tt_basis,
                               gamma = gamma, pi = pi, knots_shape = knots_shape,
                               var_c = var_c, var_e = var_e, degree = degree_shape, intercept = intercept_shape)


    var_c <- var_cUpdate(c = c, a_c = a_c, b_c = b_c)

    a <- kfeature_aUpdate(t = t, y = y, c = c, a = a, phi = phi, rho = rho, tt_basis = tt_basis,
                 gamma = gamma, pi = pi, knots_shape = knots_shape, degree = degree_shape,
                  intercept = intercept_shape, var_e = var_e, common_a = common_a)

    if (K > 1) {
      if (i <= burn_it) {
        temperature_pi <- max(1, 1 + 5 * (1 - i / burn_it) ^ 0.5)
      } else {
        temperature_pi <- 1
      }

      pi <- kfeature_piUpdate_Warp(t = t, y = y, c = c, a = a, phi = phi, rho = rho, tt_basis = tt_basis,
                                   gamma = gamma, pi = pi, knots_shape = knots_shape,
                                   degree = degree_shape, intercept = intercept_shape, var_e = var_e,
                                   alpha = alpha, tuning_param = tuning_pi, repulsive = repulsive_pi, reg = reg,
                                   temperature = temperature_pi)
    }


    var_e <- kfeature_var_eUpdate_Warp(t = t, y = y, c = c, a = a, phi = phi, rho = rho, tt_basis = tt_basis,
                                      gamma = gamma, pi = pi, knots_shape = knots_shape,
                                      degree = degree_shape, intercept = intercept_shape, a_e = a_e, b_e = b_e)

    if (warp_num > 0) {
      if (i <= burn_it) {
        temperature_rho <- max(1, 1 + 5 * (1 - i / burn_it) ^ 0.5)
      } else {
        temperature_rho <- 1
      }

      rho <- kfeature_rhoUpdate(t = t, y = y, c = c, a = a, phi = phi, rho = rho, tt_basis = tt_basis, pi = pi,
                       gamma = gamma, knots_shape = knots_shape,
                       degree = degree_shape, intercept = intercept_shape, var_e = var_e, temperature = temperature_rho)

    }

    phi_out <- kfeature_phiUpdate_NoReg(t = t, y = y, c = c, a = a, phi = phi, rho = rho, tt_basis = tt_basis,
                                 gamma = gamma, pi = pi, knots_shape = knots_shape,
                                 degree = degree_shape, intercept = intercept_shape, var_e = var_e,
                                 var_phi = var_phi, Upsilon = Upsilon, tau = tau, it_num = i,
                                 acceptance_sums = acceptance_sums)

    phi <- phi_out$phi
    acceptance_sums <- phi_out$acceptance
    tau <- phi_out$tau

    var_phi <- var_phiUpdate_NoReg(phi = phi, Upsilon = Upsilon, a_phi = a_phi, b_phi = b_phi)

    # Storing the samples
    if (i > burn_it) {
      indexing <- i - burn_it
      if (K > 1) {
        for (k in 1:K) {
          gamma_mats[[k]][indexing, ] <- gamma[, k]
          pi_mats[[k]][indexing, ] <- pi[, k]
        }
      } else {
        gamma_mats[[k]][indexing, ] <- gamma
        pi_mats[[k]][indexing, ] <- pi
      }

      c_mat[indexing, ]     <- c
      a_mat[indexing, ]     <- a
      var_mat[indexing, ]   <- c(lambda, var_c, var_e, var_phi)
      phi_mat[indexing, ]   <- matrix(t(phi), ncol = N * Q, byrow = T)

      if (warp_num > 0) {
        rho_mat[indexing, ] <- rho
      }

      # Ergodic mean and current fit
      for (j in 1:N) {
        phi_j <- phi[j, ]
        if (warp_num == 0) {
          tWarp <- tt_basis %*% phi_j
          tt_mat[indexing, ] <- tWarp
          if (K == 1) {
            shape_basis_k <- splines::bs(x = tWarp, knots = knots_shape,
                                       degree = degree_shape, intercept = intercept_shape)
            current_fit[j, ] <- c[j] + a * pi[j] * shape_basis_k %*% gamma
          } else {
            current_fit[j, ] <- c[j]
            for (k in 1: K) {
              if (k == 1){
                shape_basis_k <- splines::bs(x = tWarp, knots = knots_shape,
                                           degree = degree_shape, intercept = intercept_shape)
                current_fit[j, ] <- current_fit[j, ] + a[k] * pi[j, k] * shape_basis_k %*% gamma[, k]
              } else {
                shape_basis_k <- splines::bs(x = t, knots = knots_shape,
                                           degree = degree_shape, intercept = intercept_shape)
                current_fit[j, ] <- current_fit[j, ] + a[k] * pi[j, k] * shape_basis_k %*% gamma[, k]
              }
            }

          }
        } else {
          for (k in 1:K) {
            if (k == 1) {
              tWarp <- tt_basis %*% phi_j
              tt_mat[indexing, ] <- tWarp
              shape_basis_k <- splines::bs(x = tWarp, knots = knots_shape,
                                         degree = degree_shape, intercept = intercept_shape)
              current_fit[j, ] <- current_fit[j, ] + a * pi[j, k] * shape_basis_k %*% gamma[, k]
            } else {
              if (k > 1 & warp_num >= (k - 1)) {
                tWarp <- rho[k - 1] * (tt_basis %*% phi_j - t) + t
                tt_mats[[k - 1]][indexing, ] <- tWarp
                shape_basis_k <- splines::bs(x = tWarp, knots = knots_shape,
                                           degree = degree_shape, intercept = intercept_shape)
                current_fit[j, ] <- current_fit[j, ] + a[k] * pi[j, k] * shape_basis_k %*% gamma[, k]
              } else {
                shape_basis_k <- splines::bs(x = t, knots = knots_shape,
                                           degree = degree_shape, intercept = intercept_shape)
                current_fit[j, ] <- current_fit[j, ] + a[k] * pi[j, k] * shape_basis_k %*% gamma[, k]
              }
            }
          }
        }

        if (K == 1){
          register_fit[j, ] <- c[j] + a * pi[j] * shape_basis %*% gamma
        } else {
          register_fit[j, ] <- c[j]
          for (k in 1:K){
            register_fit[j, ] <- register_fit[j, ] + a[k] * pi[j, k] * shape_basis %*% gamma[, k]
          }
        }


        fit_mat[indexing, ((j - 1) * n + 1): (j * n)] <- current_fit[j, ]
        register_mat[indexing, ((j - 1) * n + 1): (j * n)] <- register_fit[j, ]
      }
      r1   <- (indexing - 1.0) / indexing
      r2   <- 1 / indexing
      fit  <- r1 * fit + r2 * current_fit
      fit2 <- r1 * fit2 + r2* current_fit ^ 2

      loglikelihood <- kfeature_Likelihood(t = t, y = y, c = c, a = a, gamma= gamma,
                                  pi = pi,knots_shape = knots_shape,
                                  degree = degree_shape, var_e = var_e, phi = phi, tt_basis = tt_basis,
                                  rho = rho, intercept = T, log = T)

      loglik[indexing, ] <- loglikelihood$likelihood
      loglik_i[indexing, ] <- loglikelihood$likelihoods
  }
  }

  if (warp_num > 0) {
    final <- construct_crMM(gamma = gamma_mats,
                            c = c_mat,
                            a = a_mat,
                            variance = var_mat,
                            phi = phi_mat,
                            rho = rho_mat,
                            pi = pi_mats,
                            fit_sample = fit_mat,
                            registered_fit = register_mat,
                            stochastic_time = tt_mat,
                            stochastic_time_rho = tt_mats,
                            fit = fit,
                            fit2 = fit2,
                            loglik = loglik,
                            loglik_i = loglik_i)
  } else {
    final <- construct_crMM(gamma = gamma_mats,
                            c = c_mat,
                            a = a_mat,
                            variance = var_mat,
                            phi = phi_mat,
                            pi = pi_mats,
                            fit_sample = fit_mat,
                            registered_fit = register_mat,
                            stochastic_time = tt_mat,
                            fit = fit,
                            fit2 = fit2,
                            loglik = loglik,
                            loglik_i = loglik_i)
  }
  return(final)
}




