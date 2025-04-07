source("filepath\\cov_dist_cont.R")
source("filepath\\cov_dist_cat.R")

infile = "filepath\\dataset.RData"
out_dir = "filepath"

########################### DATASET FILTERING ###########################

### read in dataset
load(infile)

### drop weekly visits and monthly visits where the LC question was not answered
dat <- dat[!(dat$visit_num %in% 1:3) &
             !is.na(dat$long_covid_have_symptoms) &
             !is.na(dat$lc_response_date) &
             dat$visit_date >= dat$lc_response_date,]

### drop participants whose first positive swab was before 11/11/2020 (12 weeks
### before LC question was introduced on 03/02/2021)
dat <- dat[!(!is.na(dat$first_swab_date) &
               dat$first_swab_date < as.numeric(as.Date("2020-11-11"))),]

### drop participants whose first positive swab was on or before CIS enrolment
dat <- dat[!(!is.na(dat$first_swab_date) &
               dat$first_swab_date <= dat$visit0_date),]

### drop participants with confirmed or suspected infection more than 14 days
### before first positive swab
dat <- dat[!((!is.na(dat$infection_date) & is.na(dat$first_swab_date)) |
               (!is.na(dat$infection_date) &  dat$first_swab_date - dat$infection_date > 14)),]

### restrict dataset to follow-up visits when participants were aged 16-64
dat <- dat[!is.na(dat$age_at_visit) & dat$age_at_visit>=16 & dat$age_at_visit<=64,]

### drop follow-up visits when participants were students
dat <- dat[dat$student==0,]

### drop follow-up visits with unknown work status
dat <- dat[!is.na(dat$work_status),]

### drop participants with missing covariate info
dat <- dat[!is.na(dat$age_at_visit),]
dat <- dat[!is.na(dat$age_visit0),]
dat <- dat[!is.na(dat$sex_visit0),]
dat <- dat[!is.na(dat$white_visit0),]
dat <- dat[!is.na(dat$gor9d_visit0),]
dat <- dat[!is.na(dat$imd_quintile_visit0),]
dat <- dat[dat$imd_quintile_visit0!="-999",]
dat <- dat[!is.na(dat$health_status_visit0),]

################################# DERIVATIONS #################################

### recalculate last visit date after applying study selection criteria
dat$last_visit_date <- ave(dat$visit_date, dat$participant_id, FUN=max)

### set infection dates to NA if later than last visit during follow-up period
dat$infection_date[dat$infection_date > dat$last_visit_date] <- NA
dat$first_swab_date[dat$first_swab_date > dat$last_visit_date] <- NA

### re-derive ever-infected flag
dat$ever_infected <- as.numeric(!is.na(dat$first_swab_date))

### re-derive flags indicating if ever had LC
dat$ever_lc_any <- ave(dat$lc_any, dat$participant_id, FUN=max)
dat$ever_lc_lim <- ave(dat$lc_lim, dat$participant_id, FUN=max)

### re-derive calendar time of first positive swab (days since 24 Jan 2020)
dat$calendar_time_infection <- dat$first_swab_date - as.numeric(as.Date("2020-01-24"))

### re-derive month and year of first positive swab
dat$yyyymm_infection <- format(as.Date(dat$first_swab_date, origin="1970-01-01"), "%Y-%m")

### re-derive month of first positive swab
dat$month_infection <- as.character(substr(dat$yyyymm_infection, 6, 7))

### re-derive time since first positive swab at each visit
dat$time_since_infection <- dat$visit_date - dat$first_swab_date

### re-derive time period of first positive swab
dat$variant_period_infection <- NA

dat$variant_period_infection[!is.na(dat$first_swab_date) &
                               dat$first_swab_date <= as.numeric(as.Date("2020-11-15"))] <- "Pre-Alpha"

dat$variant_period_infection[!is.na(dat$first_swab_date) &
                               dat$first_swab_date >= as.numeric(as.Date("2020-11-16")) &
                               dat$first_swab_date <= as.numeric(as.Date("2021-05-16"))] <- "Alpha"

dat$variant_period_infection[!is.na(dat$first_swab_date) &
                               dat$first_swab_date >= as.numeric(as.Date("2021-05-17")) &
                               dat$first_swab_date <= as.numeric(as.Date("2021-12-19"))] <- "Delta"

dat$variant_period_infection[!is.na(dat$first_swab_date) &
                               dat$first_swab_date >= as.numeric(as.Date("2021-12-20"))] <- "Omicron"

### re-derive vaccination status at first positive swab
dat$vacc_status_infection <- NA

dat$vacc_status_infection[!is.na(dat$first_swab_date) &
                            (is.na(dat$covid_vaccine_date1) |
                               ((dat$first_swab_date - dat$covid_vaccine_date1) <= 14))] <- "0"

dat$vacc_status_infection[!is.na(dat$first_swab_date) &
                            !is.na(dat$covid_vaccine_date1) &
                            ((dat$first_swab_date - dat$covid_vaccine_date1) >= 14)] <- "1"

dat$vacc_status_infection[!is.na(dat$first_swab_date) &
                            !is.na(dat$covid_vaccine_date2) &
                            ((dat$first_swab_date - dat$covid_vaccine_date2) >= 14)] <- "2"

dat$vacc_status_infection[!is.na(dat$first_swab_date) &
                            !is.na(dat$covid_vaccine_date3) &
                            ((dat$first_swab_date - dat$covid_vaccine_date3) >= 14)] <- "3"

dat$vacc_status_infection[!is.na(dat$first_swab_date) &
                            !is.na(dat$covid_vaccine_date4) &
                            ((dat$first_swab_date - dat$covid_vaccine_date4) >= 14)] <- "4"

### re-derive alternative vaccination statuses at first positive swab
dat$vacc_status_infection2 <- ifelse(dat$vacc_status_infection %in% c("1","2","3","4"), "1+",
                                     dat$vacc_status_infection)

dat$vacc_status_infection3 <- ifelse(dat$vacc_status_infection %in% c("2","3","4"), "2+",
                                     dat$vacc_status_infection)

### derive new variant period variable
dat$variant_period_infection2 <- NA
dat$variant_period_infection2[is.na(dat$variant_period_infection)] <- "0_uninfected"
dat$variant_period_infection2[!is.na(dat$variant_period_infection) & 
                                dat$variant_period_infection!="Omicron"] <- "1_pre_omicron"
dat$variant_period_infection2[!is.na(dat$variant_period_infection) & 
                                dat$variant_period_infection=="Omicron"] <- "2_omicron"

### derive new vaccination status variable
dat$vacc_status_infection4 <- NA
dat$vacc_status_infection4[is.na(dat$vacc_status_infection2)] <- "0_uninfected"
dat$vacc_status_infection4[!is.na(dat$vacc_status_infection2) &
                             dat$vacc_status_infection2=="0"] <- "1_unvaccinated"
dat$vacc_status_infection4[!is.na(dat$vacc_status_infection2) &
                             dat$vacc_status_infection2=="1+"] <- "2_vaccinated"

### derive variant period by vaccination status variable
dat$var_vacc <- NA
dat$var_vacc[dat$variant_period_infection2=="0_uninfected"] <- "0_uninfected"
dat$var_vacc[dat$variant_period_infection2=="1_pre_omicron" &
               dat$vacc_status_infection4=="1_unvaccinated"] <- "1_pre_omicron_unvaccinated"
dat$var_vacc[dat$variant_period_infection2=="1_pre_omicron" &
               dat$vacc_status_infection4=="2_vaccinated"] <- "2_pre_omicron_vaccinated"
dat$var_vacc[dat$variant_period_infection2=="2_omicron" &
               dat$vacc_status_infection4=="1_unvaccinated"] <- "3_omicron_unvaccinated"
dat$var_vacc[dat$variant_period_infection2=="2_omicron" &
               dat$vacc_status_infection4=="2_vaccinated"] <- "4_omicron_vaccinated"

### derive time since last dose at infection
dat$time_since_last_dose_infection <- ifelse(dat$vacc_status_infection=="1", dat$first_swab_date - dat$covid_vaccine_date1,
                                             ifelse(dat$vacc_status_infection=="2", dat$first_swab_date - dat$covid_vaccine_date2,
                                                    ifelse(dat$vacc_status_infection=="3", dat$first_swab_date - dat$covid_vaccine_date3,
                                                           ifelse(dat$vacc_status_infection=="4", dat$first_swab_date - dat$covid_vaccine_date4, NA))))

### derive categorical time since last dose at infection
dat$time_since_last_dose_infection_cat <- NA
dat$time_since_last_dose_infection_cat[is.na(dat$vacc_status_infection)] <- "0_uninfected"
dat$time_since_last_dose_infection_cat[!is.na(dat$vacc_status_infection) &
                                         dat$vacc_status_infection=="0"] <- "1_unvaccinated"
dat$time_since_last_dose_infection_cat[!is.na(dat$time_since_last_dose_infection) &
                                         dat$time_since_last_dose_infection>=14 &
                                         dat$time_since_last_dose_infection<=90] <- "2_14-90days"
dat$time_since_last_dose_infection_cat[!is.na(dat$time_since_last_dose_infection) &
                                         dat$time_since_last_dose_infection>=91 &
                                         dat$time_since_last_dose_infection<=180] <- "3_91-180days"
dat$time_since_last_dose_infection_cat[!is.na(dat$time_since_last_dose_infection) &
                                         dat$time_since_last_dose_infection>=181] <- "4_ge181days"

### derive categorical time since last dose at infection - version 2
dat$time_since_last_dose_infection_cat2 <- NA
dat$time_since_last_dose_infection_cat2[dat$ever_infected==0] <- "0_uninfected"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          (is.na(dat$covid_vaccine_date1) | (dat$first_swab_date < dat$covid_vaccine_date1))] <- "1_unvaccinated"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          dat$first_swab_date >= dat$covid_vaccine_date1 &
                                          (is.na(dat$covid_vaccine_date2) | (dat$first_swab_date <= (dat$covid_vaccine_date2 + 13)))] <- "2_dose1"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date2) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date2 + 14) &
                                          dat$first_swab_date <= (dat$covid_vaccine_date2 + 90) &
                                          (is.na(dat$covid_vaccine_date3) | (dat$first_swab_date <= (dat$covid_vaccine_date3 + 13)))] <- "3_dose2_14-90days"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date2) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date2 + 91) &
                                          dat$first_swab_date <= (dat$covid_vaccine_date2 + 180) &
                                          (is.na(dat$covid_vaccine_date3) | (dat$first_swab_date <= (dat$covid_vaccine_date3 + 13)))] <- "4_dose2_91-180days"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date2) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date2 + 181) &
                                          (is.na(dat$covid_vaccine_date3) | (dat$first_swab_date <= (dat$covid_vaccine_date3 + 13)))] <- "5_dose2_ge181days"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date3) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date3 + 14) &
                                          dat$first_swab_date <= (dat$covid_vaccine_date3 + 90) &
                                          (is.na(dat$covid_vaccine_date4) | (dat$first_swab_date <= (dat$covid_vaccine_date4 + 13)))] <- "6_dose3-4_14-90days"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date4) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date4 + 14) &
                                          dat$first_swab_date <= (dat$covid_vaccine_date4 + 90)] <- "6_dose3-4_14-90days"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date3) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date3 + 91) &
                                          (is.na(dat$covid_vaccine_date4) | (dat$first_swab_date <= (dat$covid_vaccine_date4 + 13)))] <- "7_dose3-4_ge91days"
dat$time_since_last_dose_infection_cat2[is.na(dat$time_since_last_dose_infection_cat2) &
                                          !is.na(dat$covid_vaccine_date4) &
                                          dat$first_swab_date >= (dat$covid_vaccine_date4 + 91)] <- "7_dose3-4_ge91days"

### derive variant period by time since last dose variable
dat$var_vacc2 <- NA
dat$var_vacc2[dat$variant_period_infection2=="0_uninfected"] <- "0_uninfected"
dat$var_vacc2[dat$variant_period_infection2=="1_pre_omicron" &
                is.na(dat$time_since_last_dose_infection)] <- "1_pre_omicron_unvaccinated"
dat$var_vacc2[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection>=14 & dat$time_since_last_dose_infection<=90] <- "2_pre_omicron_14-90days"
dat$var_vacc2[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection>=91 & dat$time_since_last_dose_infection<=180] <- "3_pre_omicron_91-180days"
dat$var_vacc2[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection>=181] <- "4_pre_omicron_ge181days"
dat$var_vacc2[dat$variant_period_infection2=="2_omicron" &
                is.na(dat$time_since_last_dose_infection)] <- "5_omicron_unvaccinated"
dat$var_vacc2[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection>=14 & dat$time_since_last_dose_infection<=90] <- "6_omicron_14-90days"
dat$var_vacc2[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection>=91 & dat$time_since_last_dose_infection<=180] <- "7_omicron_91-180days"
dat$var_vacc2[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection>=181] <- "8_omicron_ge181days"

### derive variant period by time since last dose variable - version 2
dat$var_vacc3 <- NA
dat$var_vacc3[dat$ever_infected==0] <- "0_uninfected"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="1_unvaccinated"] <- "1_pre_omicron_unvaccinated"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="2_dose1"] <- "2_pre_omicron_dose1"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="3_dose2_14-90days"] <- "3_pre_omicron_dose2_14-90days"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="4_dose2_91-180days"] <- "4_pre_omicron_dose2_91-180days"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="5_dose2_ge181days"] <- "5_pre_omicron_dose2_ge181days"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="6_dose3-4_14-90days"] <- "6_pre_omicron_dose3-4_14-90days"
dat$var_vacc3[dat$variant_period_infection2=="1_pre_omicron" &
                dat$time_since_last_dose_infection_cat2=="7_dose3-4_ge91days"] <- "7_pre_omicron_dose3-4_ge91days"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="1_unvaccinated"] <- "8_omicron_unvaccinated"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="2_dose1"] <- "9_omicron_dose1"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="3_dose2_14-90days"] <- "10_omicron_dose2_14-90days"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="4_dose2_91-180days"] <- "11_omicron_dose2_91-180days"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="5_dose2_ge181days"] <- "12_omicron_dose2_ge181days"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="6_dose3-4_14-90days"] <- "13_omicron_dose3-4_14-90days"
dat$var_vacc3[dat$variant_period_infection2=="2_omicron" &
                dat$time_since_last_dose_infection_cat2=="7_dose3-4_ge91days"] <- "14_omicron_dose3-4_ge91days"

### derive employment group at baseline
dat$work_group_visit0 <- NA
dat$work_group_visit0[dat$work_status_visit0 %in% c("1", "2", "3", "4")] <- "1_working"
dat$work_group_visit0[dat$work_status_visit0 %in% c("5")] <- "2_unemployed"
dat$work_group_visit0[dat$work_status_visit0 %in% c("6")] <- "3_notlooking"
dat$work_group_visit0[dat$work_status_visit0 %in% c("7")] <- "4_retired"
dat$work_group_visit0[dat$work_status_visit0 %in% c("8", "9", "10", "11", "12")] <- "5_student"

### re-factor work sector variable
dat$work_sector_visit0[dat$work_sector_visit0=="99"] <- "Not_working"
dat$work_sector_visit0[dat$work_sector_visit0=="999"] <- "Not_applicable"
dat$work_sector_visit0[is.na(dat$work_sector_visit0)] <- "Unknown"

### derive ever-not-looking flag
dat$ever_notlooking <- ave(dat$notlooking, dat$participant_id, FUN=max)

### derive ever-remote-response flag
dat$ever_remote_collection <- ave(dat$remote_collection, dat$participant_id, FUN=max)

### restrict dataset to participants with a positive swab during follow-up
dat_lc_vs_infected <- dat[dat$ever_infected==1,]

### restrict dataset to participants who ever had LC or without a positive swab during follow-up
dat_lc_vs_uninfected <- dat[dat$ever_lc_any==1 | dat$ever_infected==0,]

### create person-level datasets
dat_pers_lc_vs_nolc <- dat[!duplicated(dat$participant_id),]
dat_pers_lc_vs_infected <- dat_lc_vs_infected[!duplicated(dat_lc_vs_infected$participant_id),]
dat_pers_lc_vs_uninfected <- dat_lc_vs_uninfected[!duplicated(dat_lc_vs_uninfected$participant_id),]

############################# EVER-LC VS NEVER-LC #############################

### covariate distributions - continuous variables
cov_dist_cont1 <- cov.dist.cont(
  vars = c("age_visit0"),
  dataset = dat_pers_lc_vs_nolc,
  exposure = "ever_lc_any"
)

write.csv(cov_dist_cont1,
          file=paste0(out_dir, "\\cov_dist_cont_lc_vs_nolc.csv"),
          row.names=FALSE)

### covariate distributions - categorical variables
cov_dist_cat1 <- cov.dist.cat(
  vars = c("age10_visit0", "sex_visit0", "white_visit0", "gor9d_visit0",
           "imd_quintile_visit0", "health_conditions_visit0", "health_status_visit0",
           "work_group_visit0",
           "vacc_status_infection", "vacc_status_infection4",
           "time_since_last_dose_infection_cat", "time_since_last_dose_infection_cat2",
           "variant_period_infection", "variant_period_infection2",
           "var_vacc", "var_vacc2", "var_vacc3", "ever_remote_collection"),
  dataset = dat_pers_lc_vs_nolc,
  exposure = "ever_lc_any"
)

write.csv(cov_dist_cat1,
          file=paste0(out_dir, "\\cov_dist_cat_lc_vs_nolc.csv"),
          row.names=FALSE)

### covariate distributions - employment variables
cov_dist_emp1 <- cov.dist.cat(
  vars = c("work_sector_visit0", "soc_major_visit0", "self_employed_visit0"),
  dataset = dat_pers_lc_vs_nolc[!is.na(dat_pers_lc_vs_nolc$work_group_visit0) & dat_pers_lc_vs_nolc$work_group_visit0=="1_working",],
  exposure = "ever_lc_any"
)

write.csv(cov_dist_emp1,
          file=paste0(out_dir, "\\cov_dist_emp_lc_vs_nolc.csv"),
          row.names=FALSE)

######################## EVER-LC VS INFECTED WITHOUT LC ########################

### covariate distributions - continuous variables
cov_dist_cont2 <- cov.dist.cont(
  vars = c("age_visit0"),
  dataset = dat_pers_lc_vs_infected,
  exposure = "ever_lc_any"
)

write.csv(cov_dist_cont2,
          file=paste0(out_dir, "\\cov_dist_cont_lc_vs_infected.csv"),
          row.names=FALSE)

### covariate distributions - categorical variables
cov_dist_cat2 <- cov.dist.cat(
  vars = c("age10_visit0", "sex_visit0", "white_visit0", "gor9d_visit0",
           "imd_quintile_visit0", "health_conditions_visit0", "health_status_visit0",
           "work_group_visit0",
           "vacc_status_infection", "vacc_status_infection4",
           "time_since_last_dose_infection_cat", "time_since_last_dose_infection_cat2",
           "variant_period_infection", "variant_period_infection2",
           "var_vacc", "var_vacc2", "var_vacc3", "ever_remote_collection"),
  dataset = dat_pers_lc_vs_infected,
  exposure = "ever_lc_any"
)

write.csv(cov_dist_cat2,
          file=paste0(out_dir, "\\cov_dist_cat_lc_vs_infected.csv"),
          row.names=FALSE)

### covariate distributions - employment variables
cov_dist_emp2 <- cov.dist.cat(
  vars = c("work_sector_visit0", "soc_major_visit0", "self_employed_visit0"),
  dataset = dat_pers_lc_vs_infected[!is.na(dat_pers_lc_vs_infected$work_group_visit0) & dat_pers_lc_vs_infected$work_group_visit0=="1_working",],
  exposure = "ever_lc_any"
)

write.csv(cov_dist_emp2,
          file=paste0(out_dir, "\\cov_dist_emp_lc_vs_infected.csv"),
          row.names=FALSE)

############################ EVER-LC VS UNINFECTED ############################

### covariate distributions - continuous variables
cov_dist_cont3 <- cov.dist.cont(
  vars = c("age_visit0"),
  dataset = dat_pers_lc_vs_uninfected,
  exposure = "ever_lc_any"
)

write.csv(cov_dist_cont3,
          file=paste0(out_dir, "\\cov_dist_cont_lc_vs_uninfected.csv"),
          row.names=FALSE)

### covariate distributions - categorical variables
cov_dist_cat3 <- cov.dist.cat(
  vars = c("age10_visit0", "sex_visit0", "white_visit0", "gor9d_visit0",
           "imd_quintile_visit0", "health_conditions_visit0", "health_status_visit0",
           "work_group_visit0",
           "ever_remote_collection"),
  dataset = dat_pers_lc_vs_uninfected,
  exposure = "ever_lc_any"
)

write.csv(cov_dist_cat3,
          file=paste0(out_dir, "\\cov_dist_cat_lc_vs_uninfected.csv"),
          row.names=FALSE)

### covariate distributions - employment variables
cov_dist_emp3 <- cov.dist.cat(
  vars = c("work_sector_visit0", "soc_major_visit0", "self_employed_visit0"),
  dataset = dat_pers_lc_vs_uninfected[!is.na(dat_pers_lc_vs_uninfected$work_group_visit0) & dat_pers_lc_vs_uninfected$work_group_visit0=="1_working",],
  exposure = "ever_lc_any"
)

write.csv(cov_dist_emp3,
          file=paste0(out_dir, "\\cov_dist_emp_lc_vs_uninfected.csv"),
          row.names=FALSE)

############################# FOLLOW-UP TIME STATS #############################

### calculate visit-level time-since-infection stats by exposoure group
time_since_infection_stats <- data.frame(
  group = levels(as.factor(dat$grp2)),
  min = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=min)$x,
  q1 = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=function(x) quantile(x,0.25))$x,
  median = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=median)$x,
  mean = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=mean)$x,
  q3 = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=function(x) quantile(x,0.75))$x,
  max = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=max)$x,
  sd = aggregate(dat$time_since_infection[dat$ever_infected==1], by=list(dat$grp2[dat$ever_infected==1]), FUN=sd)$x
)

write.csv(time_since_infection_stats,
          file=paste(out_dir, "\\time_since_infection_stats.csv", sep=""),
          row.names=FALSE)

### calculate person-level follow-up time from first CIS visit - all participants
dat_pers_lc_vs_nolc$futime_visit0 <- dat_pers_lc_vs_nolc$last_visit_date - dat_pers_lc_vs_nolc$visit0_date

### calculate follow-up time stats from first CIS visit - all participants
futime_stats_visit0_all <- c(as.numeric(summary(dat_pers_lc_vs_nolc$futime_visit0)), sd(dat_pers_lc_vs_nolc$futime_visit0))

### calculate follow-up time stats from first CIS visit - ever-LC
futime_stats_visit0_ever_lc <- c(as.numeric(summary(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_lc_any==1])), sd(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_lc_any==1]))

### calculate follow-up time stats from first CIS visit - never-LC
futime_stats_visit0_never_lc <- c(as.numeric(summary(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_lc_any==0])), sd(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_lc_any==0]))

### calculate follow-up time stats from first CIS visit - infected without LC
futime_stats_visit0_infected_nolc <- c(as.numeric(summary(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_infected==1 & dat_pers_lc_vs_nolc$ever_lc_any==0])), sd(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_infected==1 & dat_pers_lc_vs_nolc$ever_lc_any==0]))

### calculate follow-up time stats from first CIS visit - never infected
futime_stats_visit0_uninfected <- c(as.numeric(summary(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_infected==0])), sd(dat_pers_lc_vs_nolc$futime_visit0[dat_pers_lc_vs_nolc$ever_infected==0]))

### combine into single data.frame
futime_stats <- data.frame(
  all = futime_stats_visit0_all,
  ever_lc = futime_stats_visit0_ever_lc,
  never_lc = futime_stats_visit0_never_lc,
  infected_nolc = futime_stats_visit0_infected_nolc,
  uninfected = futime_stats_visit0_uninfected
)

rownames(futime_stats) <- c("min", "q1", "median", "mean", "q3", "max", "sd")

write.csv(futime_stats, file=paste(out_dir, "\\futime_from_visit0_stats.csv", sep=""))
