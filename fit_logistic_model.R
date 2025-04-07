fit.logistic.model <- function(outcome, exposures, covariates, dataset, robust=TRUE,
                      run.lr=TRUE, run.vif=TRUE, run.aic=TRUE, run.bic=TRUE, run.auc=TRUE) {
  
  ### construct model call
  mod_call <- paste(c(
    "glm(",
    outcome,
    " ~ ",
    ifelse(!is.null(covariates),
           paste(c(exposures, covariates), collapse=" + "),
           paste(exposures, collapse=" + ")),
    ", data=dataset",
    ", family=binomial",
    ", model=FALSE)"
  ), collapse="")
  
  ### fit model
  mod <- eval(parse(text=mod_call))
  
  ### extract VCOV
  if(robust==TRUE) {vcov_obj <- vcovCL(mod, cluster=~participant_id)}
  if(robust==FALSE) {vcov_obj <- vcov(mod)}
  vcov <- as.data.frame(vcov_obj)
  
  ### extract model inferences
  if(robust==TRUE) {coeff <- as.data.frame(coeftest(mod, vcov=vcov_obj)[,])}
  if(robust==FALSE) {coeff <- as.data.frame(summary(mod)$coefficients)}
  
  ### extract LR tests
  if(run.lr==TRUE) {lr <- as.data.frame(Anova(mod))} else {lr <- NULL}
  
  ### extracts VIFs
  if(run.vif==TRUE) {
    if(!is.null(covariates)) {
      if(class(vif(mod))=="numeric") {
        vif <- as.data.frame(vif(mod))
        colnames(vif) <- "VIF"
      } else {
        vif <- as.data.frame(vif(mod)[,3])
        colnames(vif) <- "GVIF^(1/(2*Df))"
      }
    } else {vif <- NA}
  } else {vif <- NULL}
  
  ### extract information criteria
  if(run.aic==TRUE) {aic <- AIC(mod)} else {aic <- NULL}
  if(run.bic==TRUE) {bic <- BIC(mod)} else {bic <- NULL}
  
  ### extract AUROC
  if(run.auc==TRUE) {
    auc <- as.numeric(auc(roc(response=mod$y,
                              predictor=mod$fitted.values,
                              levels=0:1, direction="<")))
  } else {auc <- NULL}
  
  return(list(mod=mod, vcov=vcov_obj, coeff=coeff, lr=lr, vif=vif, aic=aic, bic=bic, auc=auc))
  
}
