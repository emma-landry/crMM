construct_crMM <- function(c, variance, fit_sample, fit, fit2, ...){
  output <- list(c = c,
                 variance = variance,
                 fit_sample = fit_sample,
                 fit = fit,
                 fit2 = fit2,
                 ...)
  class(output) <-  "crMM_Obj"
  return(output)
}
