#' Posterior estimates of the shape functions
#'
#' @description
#' `posterior_shape()` obtains posterior summaries for the shape functions. It calculates the desired summary
#' statistic/ moments of the posterior MCMC samples for the B-spline coefficients of the K features,
#' and fits the associated posterior function.
#'
#' @param crMM_samples An object of class `crMM_Obj`.
#' @param t The time points at which the data is observed.
#' @param p The number of inner knots used to define the B-splines for the common shape functions.
#' @param degree Degree of the piecewise polynomial of B-splines for the common shape functions.
#' The default is `3` for cubic splines.
#' @param intercept  If `TRUE` an intercept is included in the B-spline basis for the common
#' shape functions. The default is `FALSE`.
#' @param moments Vector with the posterior moments of interest. Options that can be included are:
#' `"mean"`, `"median"`, `"sd"`, or values between 0 and 1 for quantiles. The default is `c("mean")`.
#' @param eval_points Number of points at which to evaluated the function. Default is `1000`.
#'
#' @return A list containing a list for each of the inputted moments, that consists of:
#'
#' * `gamma` : a list of `K` elements containing the posterior summaries of the samples for
#' the B-spline coefficient associated with each of the `K` features.
#' * `f` : a list containing the evaluated functions associated with each features using the
#' estimate from `gamma`.
#'
#' @export
#'
posterior_shape <- function(crMM_samples, t, p, degree = 3, intercept = FALSE,
                            moments = c('mean'), eval_points = 1000) {

  if (!inherits(crMM_samples, "crMM_Obj")) {
   stop("'crMM_samples' must be an object of class 'crMM_Obj'.")
  }

  if (!is.null(crMM_samples$gamma)) {
    gamma_list <- crMM_samples$gamma
    K <- length(gamma_list)
  } else {
    gamma_list <- list(crMM_samples$gamma1, crMM_samples$gamma2)
    K <- 2
  }

  df <- ncol(gamma_list[[1]])

  if (intercept == F){
    if ((p + degree) != df) {
      stop("The provided B-spline parameters do not match the dimesion implied by the sample of
           shape coefficients.")
    }
  } else {
      if ((p + degree + 1) != df) {
        stop("The provided B-spline parameters do not match the dimesion implied by the sample of
           shape coefficients.")
      }
  }

  n <- length(t)
  t1 <- t[1]
  tn <- t[n]

  eval_t <- seq(t1, tn, length.out = eval_points)
  knots <- seq(t1, tn, length.out = p + 2)[-c(1, p + 2)]
  basis <- splines::bs(x = eval_t, knots = knots, degree = degree, intercept = intercept)

  number_moments <- length(moments)
  outputs <- list()

  warning_count <- 0
  for (i in 1:number_moments) {
    moment <- moments[i]
    moment_label <- moment
    valid_moment <- TRUE

    gamma_hat <- vector("list", K)
    f_hat     <- vector("list", K)

    for (k in seq_len(K)) {
      gamma_k <- gamma_list[[k]]

      if (is.character(moment)) {
        if (moment == "mean") {
          gamma_hat[[k]] <- colMeans(gamma_k)
        } else if (moment == "median") {
          gamma_hat[[k]] <- apply(gamma_k, 2, stats::median)
        } else if (moment == "sd") {
          gamma_hat[[k]] <- apply(gamma_k, 2, stats::sd)
        } else if (!is.na(as.numeric(moment))){
          q <- as.numeric(moment)
          if (q >= 0 && q <= 1) {
            gamma_hat[[k]] <- apply(gamma_k, 2, stats::quantile, probs = q)
            moment_label <- paste0(q, "_quantile")
          } else {
            warning_count <- warning_count + 1
            warning(paste(moment, " is not a valid value for the moment. The posterior summary is not
                    calculated for it."))
            valid_moment <- FALSE
            break
          }
        } else {
          warning_count <- warning_count + 1
          warning(paste(moment, " is not a valid value for the moment. The posterior summary is not
                    calculated for it."))
          valid_moment <- FALSE
          break
        }
      } else {
        if (moment >= 0 && moment <= 1) {
          gamma_hat[[k]] <- apply(gamma_k, 2, stats::quantile, probs = moment)
          moment_label <- paste0(moment, "_quantile")
        } else {
          warning_count <- warning_count + 1
          warning(paste(moment, " is not a valid value for the moment. The posterior summary is not
                    calculated for it."))
          valid_moment <- FALSE
          break
        }
      }
      f_hat[[k]] <- basis %*% gamma_hat[[k]]
    }
    if (!valid_moment) next

    outputs[[moment_label]] <- list(
      gamma = gamma_hat,
      f     = f_hat
    )
  }

  if (warning_count == number_moments) {
    stop("No valid moments were provided.")
  }

  return(outputs)

}

#' Relative Mean Integrated Square Error
#'
#' @description
#' `r_mise()` computes the relative mean integrated square error (R-MISE) between two functions.
#'
#' @details
#' The user may provide a single reference `true` function, and multiple fitted functions `fit`, as long
#' as the coordinates at which they are observed coincide. The opposite is not possible: one cannot provide
#' multiple `true` functions, but provide a single `fit` function.
#'
#'
#' @param true A vector or matrix with observed true values of a function. If a matrix, each row corresponds
#' to a different function.
#' @param fit A vector or matrix with obersved fitted values of a function. If a matrix, each row corresponds
#' to a different function.
#'
#' @return The R-MISE value. If `fit` is a matrix, then the output is a vector of R-MISE values.
#' @export
#'
r_mise <- function(true, fit) {
  if (is.matrix(true)) {
    N <- nrow(true)
    n <- ncol(true)

    if (!is.matrix(fit)) {
      stop("The dimensions of 'fit' do not match those of 'true'.")
    }

    if (nrow(fit) != N | ncol(fit) != n) {
      stop("The dimensions of 'fit' do not match those of 'true'.")
    }

    mise <- apply((true - fit) ^ 2, 1, sum)
    rmise <- mise / apply(true ^ 2, 1, sum)
  } else {
    n <- length(true)

    if (is.matrix(fit)) {
      if (ncol(fit) != n) {
        stop("The number of columns in 'fit' must match the length of 'true'.")
      }
      mise <- apply((true - fit)^2, 1, sum)
      rmise <- mise / apply(true ^ 2, 1, sum)

    } else {
        if (length(fit) != n) {
          stop("The length of 'fit' must match the length of 'true'.")
        }
      mise <- sum((true - fit) ^ 2)
      rmise <- mise / sum(true ^ 2)
    }
  }
  return(rmise)
}


#' Posterior estimates for time-transformation functions
#'
#' @description
#' `posterior_tt()` obtains the desired posterior summary statistics for quantities related to the time-
#' transformation functions. If the sample provided does not including rescaling for the second
#' feature, the posterior summary is given for the individual time-transformation B-spline coefficients,
#' the individual time-transformation functions, and the stochastic time points. When rescaling for the second feature,
#' posterior summaries of the rescaling parameter, as well as the time-transformation functions and stochastic time
#' points associated with feature 2 are also computed.
#'
#' @details
#' Note that this function is designed for two-feature mixed membership models with optional time-transformation
#' rescaling of the second feature, as used in the EEG case study (it does not generalize to K features).
#'
#' @inheritParams posterior_regression
#'
#' @return
#' When there is no rescaling of the time-transformations in the second feature, returns a list containing:
#'
#' * `phi` : a matrix for which each row corresponds to the posterior sample summary for an individual time-transformation
#' B-spline coefficient
#' * `tt_functions`: a matrix for which each row corresponds to the posterior sample summary for an individual
#' time-transformation function
#' * `stochastic_time`: a matrix for which each row corresponds to the posterior sample summary of the individual
#' stochastic time points
#'
#' If there is rescaling, the list additionally contains:
#' * `rho`: the posterior sample summary of the time-transformation rescaling parameter
#' * `tt_functions2` : a matrix for which each row corresponds to the posterior sample summary for an individual
#' time-transformation function associated with the second feature
#' * `stochastic_time`: a matrix for which each row corresponds to the posterior sample summary of the individual
#' stochastic time points associated with the second feature
#'
#' @export
#'
posterior_tt <- function(crMM_samples, h, degree = 3, intercept = FALSE,
                         moment = "mean", t_boundary, eval_points) {

  if (!inherits(crMM_samples, "crMM_Obj")) {
    stop("'crMM_samples' must be an object of class 'crMM_Obj'.")
  }

  if (intercept == FALSE) {
    df <- h + degree
  } else {
    df <- h + degree + 1
  }

  N <- ncol(crMM_samples$c)
  n <- ncol(crMM_samples$stochastic_time) / N

  if (typeof(moment) == "character") {
    if (moment == "mean") {
      phi <- apply(crMM_samples$phi, 2, mean)
      phi <- matrix(phi, nrow = N, ncol = df, byrow = T)

      tt <- apply(crMM_samples$stochastic_time, 2, mean)
      tt <- matrix(tt, nrow = N, ncol = n, byrow = T)

    } else if (moment == "median") {
      phi <- apply(crMM_samples$phi, 2, stats::median)
      phi <- matrix(phi, nrow = N, ncol = df, byrow = T)

      tt <- apply(crMM_samples$stochastic_time, 2, stats::median)
      tt <- matrix(tt, nrow = N, ncol = n, byrow = T)
    } else if (moment == "sd") {
      phi <- apply(crMM_samples$phi, 2, stats::sd)
      phi <- matrix(phi, nrow = N, ncol = df, byrow = T)

      tt <- apply(crMM_samples$stochastic_time, 2, stats::sd)
      tt <- matrix(tt, nrow = N, ncol = n, byrow = T)
    } else {
      stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
    }
  } else {
    if (moment >= 0 & moment <= 1) {
      phi <- apply(crMM_samples$phi, 2, stats::quantile, probs = moment)
      phi <- matrix(phi, nrow = N, ncol = df, byrow = T)

      tt <- apply(crMM_samples$stochastic_time, 2, stats::quantile, probs = moment)
      tt <- matrix(tt, nrow = N, ncol = n, byrow = T)
    } else {
      stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
    }
  }

  t1 <- min(t_boundary)
  tn <- max(t_boundary)

  knots <- seq(t1, tn, length.out = h + 2)[2:(h + 1)]
  eval_t <- seq(t1, tn, length.out = eval_points)
  basis <- splines::bs(x = eval_t, knots = knots, degree = degree, intercept = intercept)

  tt_functions <- basis %*% t(phi)

  if (!is.null(crMM_samples$rho)) {

    scaled_phi <- sweep(crMM_samples$phi, 1, crMM_samples$rho, FUN = "*")

    if (typeof(moment) == "character") {
      if (moment == "mean") {
        scaled_phi <- apply(scaled_phi, 2, mean)
        scaled_phi <- matrix(scaled_phi, nrow = N, ncol = df, byrow = T)

        rho <- apply(crMM_samples$rho, 2, mean)

        tt2 <- apply(crMM_samples$stochastic_time2, 2, mean)
        tt2 <- matrix(tt2, nrow = N, ncol = n, byrow = T)
      } else if (moment == "median") {
        scaled_phi <- apply(scaled_phi, 2, stats::median)
        scaled_phi <- matrix(scaled_phi, nrow = N, ncol = df, byrow = T)

        rho <- apply(crMM_samples$rho, 2, stats::median)

        tt2 <- apply(crMM_samples$stochastic_time2, 2, stats::median)
        tt2 <- matrix(tt2, nrow = N, ncol = n, byrow = T)
      } else if (moment == "sd") {
        scaled_phi <- apply(scaled_phi, 2, stats::sd)
        scaled_phi <- matrix(scaled_phi, nrow = N, ncol = df, byrow = T)

        rho <- apply(crMM_samples$rho, 2, stats::sd)

        tt2 <- apply(crMM_samples$stochastic_time2, 2, stats::sd)
        tt2 <- matrix(tt2, nrow = N, ncol = n, byrow = T)
      } else {
        stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
      }
    } else {
      if (moment >= 0 & moment <= 1) {
        scaled_phi <- apply(scaled_phi, 2, stats::quantile, probs = moment)
        scaled_phi <- matrix(scaled_phi, nrow = N, ncol = df, byrow = T)

        rho <- apply(crMM_samples$rho, 2, stats::quantile, probs = moment)

        tt2 <- apply(crMM_samples$stochastic_time2, 2, stats::quantile, probs = moment)
        tt2 <- matrix(tt2, nrow = N, ncol = n, byrow = T)
      } else {
        stop(paste(moment, " is not a valid value for the moment. The posterior summary cannot be
                    calculated for it."))
      }
    }
    tt_functions2 <- basis %*% t(scaled_phi) - rho * eval_t + eval_t
    return(list(phi = phi,
                rho = rho,
                tt_functions = tt_functions,
                tt_functions2 = tt_functions2,
                stochastic_time = tt,
                stochastic_time2 = tt2))
  } else {
    return(list(phi = phi,
                tt_functions = tt_functions,
                stochastic_time = tt))
  }
}



