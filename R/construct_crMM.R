construct_crMM <- function(gamma1, gamma2, c, variance, pi1, fit_sample, fit, fit2, ...){
  output <- list(gamma1 = gamma1,
                 gamma2 = gamma2,
                 c = c,
                 variance = variance,
                 pi1 = pi1,
                 fit_sample = fit_sample,
                 fit = fit,
                 fit2 = fit2,
                 ...)
  class(output) <-  "crMM_Obj"
  return(output)
}
