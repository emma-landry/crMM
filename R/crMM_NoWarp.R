#' MCMC sampler for functional mixed membership functional without warping
#'
#' @description
#' `crMM_NoWarp()` runs a Metropolis-within-Gibbs sampler for the analysis of functional data under the two
#' feature mixed membership model assumption. This function does not account for phase variation, see
#' [crMM_Warp()] and [crMM_WarpReg()] for samplers that include estimation of time-transformation functions.
#'
#' @param num_it The number of MCMC iterations (after burn-in).
#' @param burnin If a number between 0 and 1, the proportion of the iterations eliminated for burn-in.
#' The first `burnin * num_it` iterations are discarded, and the sampler then runs for another `num_it`
#' iterations from which samples are recorded. If an integer, the number of iterations to discard, before
#' running for another `num_it` iterations from which samples are recorded. The default is `0.2`.
#' @param t The time points at which the data is observed.
#' @param y A matrix containing the functional observations. Each row corresponds to one observation.
#' @param p The number of inner knots used to define the B-splines for the common shape functions.
#' The sampler assumes equally spaced knots between the first and last value of `t`.
#' @param degree_shape Degree of the piecewise polynomial of B-splines for the common shape functions.
#' The default is `3` for cubic splines.
#' @param intercept_shape If `TRUE` an intercept is included in the B-spline basis for the common
#' shape functions. The default is `FALSE`.
#' @param a_e,a_c,a_l Shape parameters of the Inverse Gamma priors.
#' @param b_e,b_c,b_l Rate parameters of the Inverse Gamma priors.
#' @param rescale_pi If `TRUE`, the values of the mixed membership vector pi are rescaled such that
#' at least one observation belongs entirely to each cluster. Default is `TRUE`.
#' @param tuning_pi Scaling for the shape parameter of the Metropolis proposal Dirichlet distribution
#' for sampling of pi. Default is `1000`.
#' @param alpha Shape parameter for the Dirichlet prior for pi.
#' @param label1,label2 Vector containing the indices of observations labelled as belonging to
#' feature 1 and feature 2 respectively. Default is `NULL`, in which case observations are not labelled.
#' @param wantPAF If `TRUE`, the sampler records the PAF value at each iteration (only relevant for
#' the EEG case study). Default is `FALSE`.
#' @param gamma1_init,gamma2_init Initial value the B-spline coefficients for the common shape functions.
#' Default is `NULL`, in which case the sampler obtains coefficients based on the extreme observations.
#' @param lambda1_init,lambda2_init,var_c_init,var_e_init Initial values for prior variances.
#' Default is `0.1`, except for `var_e_init` which is `1`.
#'
#' @returns An object of class `crMM_Obj`, consisting of a list with the following items:
#'
#' * `gamma1`: a matrix with `num_it` rows of the posterior samples for the B-spline coefficient
#' for the shape associated with the first feature.
#' * `gamma2`: a matrix with `num_it` rows of the posterior samples for the B-spline coefficient
#' for the shape associated with the second feature.
#' * `c`: a matrix with `num_it` rows of the posterior samples for the individual intercepts. The number
#' of columns corresponds to the number of subjects, that is, `nrow(y)`.
#' * `variance`: a matrix with `num_it` rows of the posterior samples for the variance parameters.
#' The four columns, in order, correspond to the variances for the first and second feature B-spline
#' coefficients, the intercepts, and the error term.
#' * `pi1`: a matrix with `num_it`rows of the posterior samples of the mixed membership component for the
#' first feature. The number of columns corresponds to the number of subjects, that is,
#' the number of rows in `y`.
#' * `fit_sample`: matrix with `num_it` rows of the posterior samples of individual fits. The first
#'  `lenght(t)` columns give the fit for the first observation, the second `lenght(t)` columns give
#'   the fit for the second observation, and so on for all the data.
#' * `fit1`, `fit2`: matrices with `nrow(y)` rows corresponding to the ergodic first and second
#' sample moment of the fit at points `t` for each of the observations.
#'
#' Additionally, if `wantPAF` is `TRUE` then the list also includes a matrix with `num_it` rows of the
#' posterior samples of the Peak Alpha Frequency.
#'
#' @export
#'
crMM_NoWarp <- function(num_it, burnin = 0.2, t, y, p, degree_shape = 3, intercept_shape = FALSE,
                        a_e, b_e, a_c, b_c, a_l, b_l, rescale_pi = TRUE,
                        tuning_pi = 1000, alpha, label1 = NULL, label2 = NULL, wantPAF = FALSE,
                        gamma1_init = NULL, gamma2_init = NULL, lambda1_init = 0.1,
                        lambda2_init = 0.1, var_c_init = 0.1, var_e_init = 1) {

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

  # Rescale time points ---------------------------------------------
  t1 <- t[1]
  tn <- t[n]
  t <- (t - t1) / (tn - t1)

  # Pre-calculated quantities ---------------------------------------
  knots_shape <- seq(0, 1, length.out = p + 2)[-c(1, p + 2)]
  shape_basis <- splines::bs(x = t, knots = knots_shape, degree = degree_shape, intercept = intercept_shape)
  Omega <- BsplinePrecision(knots = knots_shape, degree = degree_shape, intercept = intercept_shape)

  # Initialize parameters -------------------------------------------
  lambda1 <- lambda1_init
  lambda2 <- lambda2_init
  var_c <- var_c_init
  var_e <- var_e_init

  y_max <- apply(y, 2, max)
  y_min <- apply(y, 2, min)

  if (is.null(gamma1_init)) {
    gamma1 <- MASS::ginv(t(shape_basis) %*% shape_basis + Omega / lambda1) %*% t(shape_basis) %*% y_max
  } else {
    gamma1 <- gamma1_init
  }

  if (is.null(gamma2_init)) {
    gamma2 <- MASS::ginv(t(shape_basis) %*% shape_basis + Omega / lambda2) %*% t(shape_basis) %*% y_max
  } else {
    gamma2 <- gamma2_init
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
  #pi[, 1] <- (pi[, 1] - min(pi[, 1])) / (max(pi[, 1]) - min(pi[, 1]))
  #pi[, 2] <- 1 - pi[, 1]

  if (0 <= burnin & burnin <= 1 ){
    burn_it <- round(num_it * burnin)
  } else {
    burn_it <- burnin
  }

  total_it     <- burn_it + num_it

  # Storage matrices ------------------------------------------------
  gamma1_mat  <- matrix(nrow = num_it, ncol = p + 4, data = NA)
  gamma2_mat  <- matrix(nrow = num_it, ncol = p + 4, data = NA)
  c_mat       <- matrix(nrow = num_it, ncol = N, data = NA)
  var_mat     <- matrix(nrow = num_it, ncol = 4, data = NA)
  pi1_mat     <- matrix(nrow = num_it, ncol = N, data = NA)
  fit_mat     <- matrix(nrow = num_it, ncol = N * n, data = NA)
  current_fit <- matrix(0, nrow = N, ncol = n)
  fit         <- matrix(0, nrow = N, ncol = n)
  fit2        <- matrix(0, nrow = N, ncol = n)

  if (wantPAF == T) {
    paf_mat <- matrix(nrow = num_it, ncol = 1, data = NA)
    eval_grid <- seq(t1, tn, length = 5000)
    eval_knots <- seq(t1, tn, length = p + 2)[-c(1, p + 2)]
    spline_eval_max <- splines::bs(x = eval_grid, knots = eval_knots,
                                   degree = degree_shape, intercept = intercept_shape)
  }

  # Sampling --------------------------------------------------------
  for (i in 1:total_it) {
    c <- cUpdate_NoWarp(t = t, y = y, gamma1 = gamma1, gamma2 = gamma2, pi = pi,
                        shape_basis = shape_basis, var_c = var_c, var_e = var_e)

    var_c <- var_cUpdate(c = c, a_c = a_c, b_c = b_c)

    gamma <- gammaUpdate_NoWarp(t = t, y = y, c = c, pi = pi, shape_basis = shape_basis,
                                Omega = Omega, lambda1 = lambda1, lambda2 = lambda2, var_e = var_e)
    gamma1 <- gamma[[1]]
    gamma2 <- gamma[[2]]

    lambda <- lambdaUpdate(gamma1 = gamma1, gamma2 = gamma2, a_l = a_l, b_l = b_l, Omega = Omega)
    lambda1 <- lambda[[1]]
    lambda2 <- lambda[[2]]

    pi <- piUpdate_NoWarp(t = t, y = y, c = c, gamma1 = gamma1, gamma2 = gamma2, pi = pi,
                          shape_basis = shape_basis, var_e = var_e, alpha = alpha, rescale = rescale_pi,
                          label1 = label1, label2 = label2, tuning_param = tuning_pi)

    var_e <- var_eUpdate_NoWarp(t = t, y = y, c = c, gamma1 = gamma1, gamma2 = gamma2, pi = pi,
                                shape_basis = shape_basis, a_e = a_e, b_e = b_e)

    if (wantPAF == T) {
      shape1_eval <- gamma1 %*% t(spline_eval_max)
      peak_location <- eval_grid[(order(shape1_eval, decreasing = T)[1])]
    }

    # Storing the samples
    if (i > burn_it) {
      indexing <- i - burn_it
      gamma1_mat[indexing, ] <- gamma1
      gamma2_mat[indexing, ] <- gamma2
      c_mat[indexing, ]      <- c
      pi1_mat[indexing, ]    <- pi[, 1]
      var_mat[indexing, ]    <- c(lambda1, lambda2, var_c, var_e)

      if (wantPAF == T) {
        paf_mat[indexing, ] <- peak_location
      }

      # Ergodic mean and current fit
      for (j in 1:N) {
        current_fit[j, ] <- c[j] + pi[j, 1] * shape_basis %*% gamma1 + pi[j, 2] * shape_basis %*% gamma2
        fit_mat[indexing, ((j - 1) * n + 1): (j * n)] <- current_fit[j, ]
      }
      r1   <- (indexing - 1.0) / indexing
      r2   <- 1 / indexing
      fit  <- r1 * fit + r2 * current_fit
      fit2 <- r1 * fit2 + r2 * current_fit ^ 2
    }
  }

  if (wantPAF == T) {
    final <- construct_crMM(gamma1 = gamma1_mat,
                            gamma2 = gamma2_mat,
                            c = c_mat,
                            variance = var_mat,
                            pi1 = pi1_mat,
                            fit_sample = fit_mat,
                            fit = fit,
                            fit2 = fit2,
                            PAF = paf_mat)
  } else {
    final <- construct_crMM(gamma1 = gamma1_mat,
                            gamma2 = gamma2_mat,
                            c = c_mat,
                            variance = var_mat,
                            pi1 = pi1_mat,
                            fit_sample = fit_mat,
                            fit = fit,
                            fit2 = fit2)
  }
  return(final)
}




