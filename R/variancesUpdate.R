var_cUpdate <- function(c, a_c, b_c) {
  N <- length(c)
  var_c <- 1 / stats::rgamma(n = 1, shape = a_c + N / 2, rate = b_c + sum(c ^ 2) / 2)
  return(var_c)
}
