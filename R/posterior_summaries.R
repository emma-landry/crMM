#' Posterior estimates of the shape functions
#'
#' @description
#' `posterior_shape()` obtains posterior summaries for the shape functions. It calculates the desired summary
#' statistic/ moments of the posterior MCMM samples for the B-spline coefficients of the two features,
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
#' `'mean'`, `'median'`, `'sd'`, or values between 0 and 1 for quantiles. Default is `c('mean')`.
#' @param eval_points Number of points at which to evaluated the function. Default is `1000`.
#'
#' @return A list containing a list for each of the inputted moments, that consists of:
#'
#' * `gamma1` : the posterior summary of the sample for the B-spline coefficient associated with feature 1.
#' * `gamma2` : the posterior summary of the sample for the B-spline coefficient associated with feature 2.
#' * `f1` : the evaluated function associated with feature 1 using the estimate `gamma1`.
#' * `f2` : the evaluated function associated with feature 1 using the estimate `gamma2`.
#'
#' @export
#'
posterior_shape <- function(crMM_samples, t, p, degree = 3, intercept = FALSE,
                            moments = c('mean'), eval_points = 1000) {
  if (!inherits(crMM_samples, "crMM_Obj")) {
    stop("`crMM_samples` must be an object of class `crMM_Obj`.")
  }

  df <- ncol(crMM_samples$gamma1)

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

  for (i in 1:number_moments) {
    moment <- moments[i]
    warning_count <- 0
    if (typeof(moment) == "character") {
      if (moment == "mean") {
        gamma1 <- apply(crMM_samples$gamma1, 2, mean)
        gamma2 <- apply(crMM_samples$gamma2, 2, mean)
      } else if (moment == "median") {
        gamma1 <- apply(crMM_samples$gamma1, 2, median)
        gamma2 <- apply(crMM_samples$gamma2, 2, median)
      } else if (moment == "sd") {
        gamma1 <- apply(crMM_samples$gamma1, 2, sd)
        gamma2 <- apply(crMM_samples$gamma2, 2, sd)
      } else {
        warning_count <- warning_count + 1
        warning(paste(moment, " is not a valid value for the moment. The posterior summary is not
                    calculated for it."))
        next
      }
    } else {
      if (moment >= 0 & moment <= 1) {
        gamma1 <- apply(crMM_samples$gamma1, 2, stats::quantile, probs = moment)
        gamma2 <- apply(crMM_samples$gamma2, 2, stats::quantile, probs = moment)
        moment <- paste(moment, "quantile")
      } else {
        warning_count <- warning_count + 1
        warning(paste(moment, " is not a valid value for the moment. The posterior summary is not
                    calculated for it."))
        next
      }
    }

    f1 <- basis %*% gamma1
    f2 <- basis %*% gamma2

    outputs$moment <- list(gamma1 = gamma1,
                           gamma2 = gamma2,
                           f1 = f1,
                           f2 = f2)
  }

  if (warning_count == number_moments) {
    stop("No valid moments were provided.")
  }

  return(outputs)

}
