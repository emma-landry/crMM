BsplinePrecision <- function(df = NULL, knots = NULL, degree = 3, intercept = FALSE) {
  if (intercept == FALSE) {
    p <- degree
  } else {
    p <- degree + 1
  }

  if (!is.null(df)) {
    if (df <= p ){
      stop("'df' needs to be greater than degree (or degree + 1 if 'intercept' = T).")
    }
  } else {
    if (is.null(knots)) {
      stop("An input needs to be provided for either 'knots' or 'df'.")
    } else {
      h <- length(knots)
      df <- h + p
    }
  }

  Precision <- matrix(0, nrow = df, ncol = df)
  diag(Precision) <- c(rep(1, df - 1), 1 / 2)
  diag(Precision[1:(df - 1), 2:df]) <- rep(-1, df - 1)
  Precision <- Precision + t(Precision)

  return(Precision)
}

BsplineCov <- function(df = NULL, knots = NULL, degree = 3, intercept = FALSE) {
  if (intercept == FALSE) {
    p <- degree
  } else {
    p <- degree + 1
  }

  if (!is.null(df)) {
    if (df <= p ){
      stop("'df' needs to be greater than degree (or degree + 1 if 'intercept' = T).")
    }
  } else {
    if (is.null(knots)) {
      stop("An input needs to be provided for either 'knots' or 'df'.")
    } else {
      h <- length(knots)
      df <- h + p
    }
  }

  Covariance <- matrix(1, nrow = df, ncol = df)

  for (i in 2:df) {
    Covariance[i, i:df] <- i
    Covariance[i:df, i] <- i
  }

  return(Covariance)
}


#' Compute the Relative Integrated Mean Squared Error
#'
#' @description
#' `compute_rmse()` computes the Relative Integrated Mean Sqaured Error (R-MISE) between
#' a reference function and its estimate, evaluated at the same input points.
#'
#' @param true Discrete evaluations of the true function
#' @param fit Discrete evaluations of the estimated function
#'
#' @return
#' The R-MISE value
#'
#' @export
#'
compute_rmse <- function(true, fit) {
  n <- length(true)

  if (length(fit) != n) {
    stop("'true' and 'fit' must have the same length.")
  }

  rmse <- sum((true- fit) ^ 2)
  rmse <- rmse / sum(true ^ 2)
  return(rmse)
}
