library(splines)
library(survival)
library(lmtest)
library(car)

source("filepath\\fit_clogit_model.R")

### set file paths
infile = "filepath\\dataset_filtered.RData"
out_dir = "filepath"

### define time-varying covariates to control for in models
covs = c(
  "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
  "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
  "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
)

### define variables to interact with exposure
int_vars = c(
  "age2grp",
  "age3grp",
  "age4grp",
  "sex_visit0",
  "health_status_visit0",
  "white_visit0",
  "imd_quintile_visit0",
  "remote_collection",
  "reinfected"
)

### assign names to interaction variables to use in output file names
int_names = c(
  "age2grp",
  "age3grp",
  "age4grp",
  "sex",
  "health",
  "white",
  "imd",
  "mode",
  "reinfected"
)

### include interaction between calendar time and effect modifier?
include_time_int <- c(
  FALSE,
  FALSE,
  FALSE,
  FALSE,
  FALSE,
  TRUE,
  TRUE,
  FALSE,
  FALSE
)

### read in dataset
load(infile)

### rename dataset
orig <- dat
rm(dat); gc()

### create four-, thee- and two-group versions of age
orig$age4grp <- orig$age10

orig$age3grp <- orig$age10
orig$age3grp[orig$age3grp %in% c("(15,24]", "(24,34]")] <- "(15,34]"

orig$age2grp <- ifelse(orig$age_at_visit>=50, "(49,64]", "(15,49]")

### loop over each interaction variable
for(i in 1:length(int_vars)) {
  
  ### reset dataset
  dat <- orig
  
  ### subset dataset to exclude missing values of interaction variable
  if(int_vars[i] %in% c("work_sector", "work_sector_coarse")) {dat <- dat[dat$work_sector!="Unknown",]}
  if(int_vars[i]=="soc_major") {dat <- dat[dat$soc_major!="Unknown",]}
  
  ### define covariate list for baseline model (i.e. adjusted model without interaction with exposure)
  time_int <- ifelse(include_time_int[i]==TRUE,
                     paste0("ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) : ", int_vars[i]),
                     NA)
  
  main_effect <- int_vars[i]
  
  covs_baseline <- c(covs,
                     time_int,
                     main_effect)
  
  covs_baseline <- as.character(na.omit(covs_baseline))
  
  ### fit baseline model
  mod_baseline <- fit.clogit.model(outcome="notlooking",
                                   exposures="exposure",
                                   covariates=covs_baseline,
                                   id="participant_id",
                                   dataset=dat,
                                   run.wald=FALSE,
                                   run.aic=FALSE,
                                   run.bic=FALSE)
  
  ### add interaction with exposure to baseline model
  exposure_int <- paste0("exposure : ", int_vars[i])
  
  covs_int <- c(covs_baseline,
                exposure_int)
  
  ### fit adjusted model including interaction with exposure
  mod_int <- fit.clogit.model(outcome="notlooking",
                              exposures="exposure",
                              covariates=covs_int,
                              id="participant_id",
                              dataset=dat,
                              run.wald=FALSE,
                              run.aic=FALSE,
                              run.bic=FALSE)
  
  ### save and write out LR test results
  mod_lr <- as.data.frame(lrtest(mod_int$mod, mod_baseline$mod))
  write.csv(mod_lr, file=paste0(out_dir, "\\lr_", int_names[i], ".csv"), row.names=FALSE)
  
  ### save and write out model coefficients
  mod_coeff <- mod_int$coeff
  write.csv(mod_coeff, file=paste0(out_dir, "\\coeffs_", int_names[i], ".csv"))
  
  ### save and write out variance-covariance matrix
  mod_vcov <- mod_int$vcov
  write.csv(mod_vcov, file=paste0(out_dir, "\\vcov_", int_names[i], ".csv"))
  
  ### restrict coeffs and VCOV to main effects and interactions of interest
  mod_coeff <- mod_coeff[grep("exposure", rownames(mod_coeff)),]
  mod_vcov <- mod_vcov[grep("exposure", rownames(mod_vcov)),
                       grep("exposure", colnames(mod_vcov))]
  
  ### count number of levels in model for main and interaction variables
  n_main_levels <- nlevels(as.factor(dat[["exposure"]])) - 1
  n_int_levels <- nlevels(as.factor(dat[[int_vars[i]]])) - 1
  
  ### set up empty lists to store ORs and CIs
  or_list <- NULL
  lcl_list <- NULL
  ucl_list <- NULL
  
  ### calculate ORs and CIs for main effects
  main_coeff <- mod_coeff[1:n_main_levels, 1]
  main_var <- diag(as.matrix(mod_vcov))[1:n_main_levels]
  main_se <- sqrt(main_var)
  
  or_list[[1]] <- exp(main_coeff)
  lcl_list[[1]] <- exp(main_coeff - 1.96 * main_se)
  ucl_list[[1]] <- exp(main_coeff + 1.96 * main_se)
  
  ### calculate ORs and CIs for each level of interaction variable
  for(j in 1:n_int_levels) {
    
    int_coeff <- main_coeff +
      mod_coeff[(j*n_main_levels+1):((j+1)*n_main_levels), 1]
    
    int_var <- main_var +
      diag(as.matrix(mod_vcov))[(j*n_main_levels+1):((j+1)*n_main_levels)] +
      2 * diag(as.matrix(mod_vcov[(j*n_main_levels+1):((j+1)*n_main_levels), 1:n_main_levels]))
    
    int_se <- sqrt(int_var)
    
    or_list[[j+1]] <- exp(int_coeff)
    lcl_list[[j+1]] <- exp(int_coeff - 1.96 * int_se)
    ucl_list[[j+1]] <- exp(int_coeff + 1.96 * int_se)
    
  }
  
  ### combine ORs and CIs into single data.frame
  or_df <- data.frame(main_var=rep(levels(as.factor(dat[["exposure"]]))[-1], n_int_levels+1),
                      int_var=rep(levels(as.factor(dat[[int_vars[i]]])), each=n_main_levels),
                      or=unlist(or_list),
                      lcl=unlist(lcl_list),
                      ucl=unlist(ucl_list))
  
  ### write out ORs and CIs
  write.csv(or_df, file=paste0(out_dir, "\\or_", int_names[i], ".csv"), row.names=FALSE)
  
  ### print iteration when complete
  print(i)
  
}
