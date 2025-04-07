fit.clogit.model <- function(outcome, exposures, covariates, id, dataset,
                      run.wald=TRUE, run.aic=TRUE, run.bic=TRUE) {
  
  ### construct model call
  mod_call <- paste(c(
    "clogit(",
    outcome,
    " ~ ",
    ifelse(!is.null(covariates),
           paste(c(exposures, covariates, paste0("strata(", id, ")")), collapse=" + "),
           paste(c(exposures, paste0("strata(", id, ")")), collapse=" + ")),
    ", data=dataset",
    ", model=FALSE)"
  ), collapse="")
  
  ### fit model
  mod <- eval(parse(text=mod_call))
  
  ### extract model inferences
  coeff <- as.data.frame(summary(mod)$coefficients)
  
  ### extract VCOV
  vcov <- as.data.frame(vcov(mod))
  
  ### extract Wald tests
  if(run.wald==TRUE) {wald <- as.data.frame(Anova(mod, test.statistic="Wald"))} else {wald <- NULL}
  
  ### extract information criteria
  if(run.aic==TRUE) {aic <- AIC(mod)} else {aic <- NULL}
  if(run.bic==TRUE) {bic <- BIC(mod)} else {bic <- NULL}
  
  return(list(mod=mod, vcov=vcov, coeff=coeff, wald=wald, aic=aic, bic=bic))
  
}
