identityTT <- function(Boundary.knots, knots = NULL, df = NULL, degree = 3, intercept = FALSE) {
  if (length(Boundary.knots) != 2) {
    stop("'Boundary.knots' needs to have length 2.")
  }

  t1 <- min(Boundary.knots)
  tn <- max(Boundary.knots)

  if (intercept == FALSE) {
    p <- degree
  } else {
    p <- degree + 1
  }

  if (!is.null(knots)) {
    if (min(knots) <= t1) {
      stop("At least one element in 'knots' is less than or equal to the lower boundary knot.")
    }

    if (max(knots) <= tn) {
      stop("At least one element in 'knots' is greater than or equal to the upper boundary knot.")
    }

    h <- length(knots)

    if (h <1) {
      stop("'knots' needs to be of length greater than or equal to 1.")
    }

    padded_knots <- c(rep(t1, p), knots, rep(tn, p))
  } else {
    if (is.null(df)){
      stop("An input needs to be provided for either 'knots' or 'df'.")
    } else {
      if (df <= p ){
        stop("'df' needs to be greater than degree (or degree + 1 if 'intercept' = T).")
      }

      h <- df - p
      knots <- seq(t1, tn, length.out = h + 2)[2:(h + 1)]
      padded_knots <- c(rep(t1, p), knots, rep(tn, p))
    }
  }

  Ups <- rep(t1, h + p)

  for (i in (1:(h + p -1))) {
    Ups[i + 1] <- (padded_knots[i + p] - padded_knots[i + 1]) / (p - 1) + Ups[i]
  }
  return(Ups)
}

jupp <- function(x) {
  if (is.unsorted(x)) {
    stop("The entries in 'x' need to be increasing.")
  }
  p <- length(x)
  y <- x
  y[2:(p - 1)] <- log(x[3:p] - x[2:(p - 1)]) - log(x[2:(p - 1)] - x[1:(p-2)])
  return(y)
}

juppinv <- function(y) {
  p <- length(y)
  x <- y
  a <- c(1, exp(cumsum(y[2:(p-1)])))
  h <- (y[p] - y[1]) * a / sum(a[1:(p - 1)])
  x[2:(p-1)] <- y[1] + cumsum(h[1:(p - 2)])
  return(x)
}

