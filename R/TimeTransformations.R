identityTT <- function(Boundary.knots, knots = NULL, df = NULL, degree = 3, intercept = FALSE){
  t1 <- Boundary.knots[1]
  tn <- Boundary.knots[2]

  if(intercept == FALSE){
    p <- degree
  } else{
    p <- degree + 1
  }

  if (!is.null(knots)){
    h <- length(knots)
    padded_knots <- c(rep(t1, p), knots, rep(tn, p))
  } else {
    if (is.null(df)){
      stop("An input needs to be provided for either 'knots' or 'df'.")
    }else{
      h <- df - p
      knots <- seq(t1, tn, length.out = h + 2)[2:(h + 1)]
      padded_knots <- c(rep(t1, p), knots, rep(tn, p))
    }
  }

  Ups <- rep(t1, h + p)

  for (i in (1: (h + p -1))){
    Ups[i + 1] <- (padded_knots[i + p] - padded_knots[i + 1]) / (p - 1) + Ups[i]
  }
  return(Ups)
}


