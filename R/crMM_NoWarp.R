crMM_NoWarp <- function(num_it, burnin = 0.2, t, y, p, degree_shape = 3, intercept_shape = F,
                        a_e, b_e, a_c, b_c, a_l, b_l, rescale_pi = T,
                        tuning_pi, alpha, label1 = NULL, label2 = NULL, wantPAF = T,
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
    pi[label2, 1] <- 1
    pi[label2, 2] <- 0
  }
  pi[, 1] <- (pi[, 1] - min(pi[, 1])) / (max(pi[, 1]) - min(pi[, 1]))
  pi[, 2] <- 1 - pi[, 1]

  # Storage matrices ------------------------------------------------
  burn_it     <- round(num_it * burnin)
  total_it    <- burn_it + num_it
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
      fit2 <- r1 * fit2 + r2* current_fit ^ 2
    }
  }
  final <- list("gamma1" = gamma1_mat,
                "gamma2" = gamma2_mat,
                "c" = c_mat,
                "variance" = var_mat,
                "pi1" = pi1_mat,
                "fit_sample" = fit_mat,
                "fit" = fit,
                "fit2" = fit2)

  if (wantPAF == T) {
    final$PAF <- paf_mat
  }

  return(final)
}




