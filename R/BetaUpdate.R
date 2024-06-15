BetaUpdate <- function(phi, X, Upsilon, B0, V0, var_phi, U, it_num){
  Q <- ncol(phi)
  N <- nrow(phi)
  l <- ncol(X)

  if(length(Upsilon) != Q) {
    stop("The length of 'Upsilon' must match the number of columns in 'phi'.")
  }

  if (nrow(X) != N) {
    stop("The number of rows in 'X' must match the number of rows in 'phi'.")
  }

  if (nrow(B0) != l) {
    stop("The number of rows in 'B0' must match the number of columns in 'X'.")
  }

  if (ncol(B0) != Q - 2) {
    stop("The number of columns in 'B0' must be 'ncol(phi) - 2'.")
  }

  if (ncol(V0) != l | nrow(V0) != l) {
    stop("'V0' must be a square matrix of dimension matching the number of columns in 'X'.")
  }

  if (ncol(U) != l | nrow(U) != l) {
    stop("'U' must be a square matrix of dimension matching the number of columns in 'X'.")
  }

  jupp_mean <- jupp(Upsilon)[-c(1, Q)]
  phiT <- t(phi)

  #if (it_num > 10000) {
  if (TRUE %in% is.na(phiT)) {
    cat("The transposed phi are", "\n")
    print(phiT)
  }

  jupp_phiT <- apply(phiT, 2, jupp)[-c(1, Q),]

  #if (it_num > 10000) {
  if (TRUE %in% is.na(jupp_phiT)) {
    cat("The transposed phi after Jupp are", "\n")
    print(jupp_phiT)
  }

  H <- t(jupp_phiT)
  posterior_mean <- U %*% (t(X) %*% (H - jupp_mean) + V0 %*% B0)
  B <- matrixNormal::rmatnorm(s = 1, M = posterior_mean, U = U, V = var_phi * diag(Q - 2))
  return(B)
}
