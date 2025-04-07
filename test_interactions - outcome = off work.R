library(splines)
library(survival)
library(lmtest)
library(car)
library(zoo)

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
  "sex_visit0",
  "health_conditions_visit0",
  "white_visit0",
  "imd_quintile_visit0",
  "remote_collection",
  "reinfected",
  "work_sector",
  "work_sector_coarse",
  "soc_major",
  "self_employed"
)

### assign names to interaction variables to use in output file names
int_names = c(
  "age2grp",
  "age3grp",
  "sex",
  "health",
  "white",
  "imd",
  "mode",
  "reinfected",
  "sector",
  "sector_coarse",
  "occupation",
  "selfemp"
)

### include interaction between calendar time and effect modifier?
include_time_int <- c(
  FALSE,
  FALSE,
  FALSE,
  FALSE,
  TRUE,
  TRUE,
  FALSE,
  FALSE,
  TRUE,
  TRUE,
  TRUE,
  TRUE
)

### read in dataset
load(infile)

### rename dataset
orig <- dat
rm(dat); gc()

### restrict sample to visits when participants were in employment
orig <- orig[orig$working==1,]

### restrict sample to visits from 1 Oct 2021 onwards (end of furlough)
orig <- orig[orig$visit_date>=as.numeric(as.Date("2021-10-01")),]

### create four-, thee- and two-group versions of age
orig$age4grp <- orig$age10

orig$age3grp <- orig$age10
orig$age3grp[orig$age3grp %in% c("(15,24]", "(24,34]")] <- "(15,34]"

orig$age2grp <- ifelse(orig$age_at_visit>=50, "(49,64]", "(15,49]")

################################################################################

### sort dataset by participant ID, visit date, and visit ID
orig <- orig[with(orig, order(participant_id, visit_date, -xtfrm(visit_id))),]

### create function for LOCF
locf.fun <- function(x) {na.locf(x, na.rm=FALSE)}

### replace 'Unknown' occupation category with NA
orig$soc_major[orig$soc_major=="Unknown"] <- NA

### impute missing values of work sector and occupation group
orig$work_sector <- ave(orig$work_sector, orig$participant_id, FUN=locf.fun)
orig$soc_major <- ave(orig$soc_major, orig$participant_id, FUN=locf.fun)

### reset unknown occupation from NA to 'Unknown'
orig$soc_major[is.na(orig$soc_major)] <- "Unknown"

################################################################################

### assign NA work sector to 'unknown' group
orig$work_sector[is.na(orig$work_sector)] <- "Unknown"

### group social care in with healthcare
orig$work_sector[orig$work_sector %in% c("2", "3")] <- "02_03"

### group smaller work sectors in with 'other'
orig$work_sector[orig$work_sector %in% c("7", "8", "9", "10", "13", "14", "15")] <- "Other"

### append leading 0 to single-digit work sectors
orig$work_sector <- ifelse(nchar(orig$work_sector)==1, paste0("0", orig$work_sector), orig$work_sector)

### derive coarse version of work sector
orig$work_sector_coarse <- "0_other"
orig$work_sector_coarse[orig$work_sector=="01"] <- "1_education"
orig$work_sector_coarse[orig$work_sector=="02_03"] <- "2_health_social_care"
orig$work_sector_coarse[orig$work_sector=="Unknown"] <- "3_unknown"

### group SOC major 0 in with 'unknown'
orig$soc_major[orig$soc_major=="0"] <- "Unknown"

### group together SOC major 8 and 9
orig$soc_major[orig$soc_major %in% c("8", "9")] <- "8_9"

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
  mod_baseline <- fit.clogit.model(outcome="offwork",
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
  mod_int <- fit.clogit.model(outcome="offwork",
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
