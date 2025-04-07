library(splines)
library(survival)
library(sandwich)
library(lmtest)
library(car)
library(pROC)
library(ggplot2)

########################### SOURCE FUNCTIONS ###########################

source("filepath\\fit_logistic_model.R")
source("filepath\\fit_clogit_model.R")
source("filepath\\modelling_logistic.R")
source("filepath\\modelling_clogit.R")
source("filepath\\modelling_clogit_lt40w.R")
source("filepath\\modelling_clogit_lt52w.R")
source("filepath\\modelling_clogit_ge12w.R")
source("filepath\\modelling_clogit_grp4.R")

########################### SET GLOBAL PARAMETERS ###########################

infile = "filepath\\dataset.RData"
root = "filepath\\Modelling"
chart_colours = c("#7fbf76", "#1b7837", "#762a83", "#af8dc3")

########################### DATASET FILTERING ###########################

### read in dataset
load(infile)

### initial sample waterfall
n_visits <- nrow(dat)
names(n_visits) <- "Initial sample"
n_people <- nrow(dat[!duplicated(dat$participant_id),])
names(n_people) <- names(n_visits)

### drop weekly visits and monthly visits where the LC question was not answered
dat <- dat[!(dat$visit_num %in% 1:3) &
             !is.na(dat$long_covid_have_symptoms) &
             !is.na(dat$lc_response_date) &
             dat$visit_date >= dat$lc_response_date,]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Monthly follow-up visits where LC question was answered"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### drop participants whose first positive swab was before 11/11/2020 (12 weeks
### before LC question was introduced on 03/02/2021)
dat <- dat[!(!is.na(dat$first_swab_date) &
               dat$first_swab_date < as.numeric(as.Date("2020-11-11"))),]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Participants without a positive swab, or first positive swab from 11/11/2020"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### drop participants whose first positive swab was on or before CIS enrolment
dat <- dat[!(!is.na(dat$first_swab_date) &
               dat$first_swab_date <= dat$visit0_date),]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Participants without a positive swab, or first positive swab after CIS enrolment visit"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### drop participants with confirmed or suspected infection more than 14 days
### before first positive swab
dat <- dat[!((!is.na(dat$infection_date) & is.na(dat$first_swab_date)) |
               (!is.na(dat$infection_date) &  dat$first_swab_date - dat$infection_date > 14)),]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Participants without confirmed or suspected infection >14 days before first positive swab"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### restrict dataset to follow-up visits when participants were aged 16-64
dat <- dat[!is.na(dat$age_at_visit) & dat$age_at_visit>=16 & dat$age_at_visit<=64,]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Follow-up visits where participants were aged 16-64"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### drop follow-up visits when participants were students
dat <- dat[dat$student==0,]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Follow-up visits where participants were not students"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### drop follow-up visits with unknown work status
dat <- dat[!is.na(dat$work_status),]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Follow-up visits where participants had a known employment status"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### drop participants with missing covariate info
dat <- dat[!is.na(dat$age_at_visit),]
dat <- dat[!is.na(dat$age_visit0),]
dat <- dat[!is.na(dat$sex_visit0),]
dat <- dat[!is.na(dat$white_visit0),]
dat <- dat[!is.na(dat$gor9d_visit0),]
dat <- dat[!is.na(dat$imd_quintile_visit0),]
dat <- dat[dat$imd_quintile_visit0!="-999",]
dat <- dat[!is.na(dat$health_status_visit0),]

n_visits <- c(n_visits, nrow(dat))
names(n_visits)[length(n_visits)] <- "Participants with complete covariate info"
n_people <- c(n_people, nrow(dat[!duplicated(dat$participant_id),]))
names(n_people) <- names(n_visits)

### save sample waterfall to working directory
sample_waterfall <- as.data.frame(cbind(n_visits, n_people))
colnames(sample_waterfall) <- c("Visits", "Participants")
write.csv(sample_waterfall, file=paste0(root, "\\sample_waterfall.csv"))

########################### VARIABLE DERIVATIONS ###########################

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

### derive exposure variable - infection/LC status by time since infection
dat$exposure <- NA

dat$exposure[is.na(dat$time_since_infection) | dat$time_since_infection<0] <- "00_uninfected"

dat$exposure[dat$time_since_infection>=0 & dat$time_since_infection<84] <- "01_lt12w"

dat$exposure[dat$time_since_infection>=84 & dat$time_since_infection<126 & dat$grp2=="3_previous_infection_ge12w"] <- "02_12-18w_nolc"
dat$exposure[dat$time_since_infection>=126 & dat$time_since_infection<168 & dat$grp2=="3_previous_infection_ge12w"] <- "03_18-24w_nolc"
dat$exposure[dat$time_since_infection>=168 & dat$time_since_infection<210 & dat$grp2=="3_previous_infection_ge12w"] <- "04_24-30w_nolc"
dat$exposure[dat$time_since_infection>=210 & dat$time_since_infection<280 & dat$grp2=="3_previous_infection_ge12w"] <- "05_30-40w_nolc"
dat$exposure[dat$time_since_infection>=280 & dat$time_since_infection<364 & dat$grp2=="3_previous_infection_ge12w"] <- "06_40-52w_nolc"
dat$exposure[dat$time_since_infection>=364 & dat$grp2=="3_previous_infection_ge12w"] <- "07_ge52w_nolc"

dat$exposure[dat$time_since_infection>=84 & dat$time_since_infection<126 & dat$grp2=="4_long_covid_any"] <- "08_12-18w_currlc"
dat$exposure[dat$time_since_infection>=126 & dat$time_since_infection<168 & dat$grp2=="4_long_covid_any"] <- "09_18-24w_currlc"
dat$exposure[dat$time_since_infection>=168 & dat$time_since_infection<210 & dat$grp2=="4_long_covid_any"] <- "10_24-30w_currlc"
dat$exposure[dat$time_since_infection>=210 & dat$time_since_infection<280 & dat$grp2=="4_long_covid_any"] <- "11_30-40w_currlc"
dat$exposure[dat$time_since_infection>=280 & dat$time_since_infection<364 & dat$grp2=="4_long_covid_any"] <- "12_40-52w_currlc"
dat$exposure[dat$time_since_infection>=364 & dat$grp2=="4_long_covid_any"] <- "13_ge52w_currlc"

dat$exposure[dat$time_since_infection>=84 & dat$time_since_infection<126 & dat$grp2=="5_previous_long_covid_any"] <- "14_12-18w_prevlc"
dat$exposure[dat$time_since_infection>=126 & dat$time_since_infection<168 & dat$grp2=="5_previous_long_covid_any"] <- "15_18-24w_prevlc"
dat$exposure[dat$time_since_infection>=168 & dat$time_since_infection<210 & dat$grp2=="5_previous_long_covid_any"] <- "16_24-30w_prevlc"
dat$exposure[dat$time_since_infection>=210 & dat$time_since_infection<280 & dat$grp2=="5_previous_long_covid_any"] <- "17_30-40w_prevlc"
dat$exposure[dat$time_since_infection>=280 & dat$time_since_infection<364 & dat$grp2=="5_previous_long_covid_any"] <- "18_40-52w_prevlc"
dat$exposure[dat$time_since_infection>=364 & dat$grp2=="5_previous_long_covid_any"] <- "19_ge52w_prevlc"

### group uninfected and infected without LC (i.e. no LC vs. current LC vs. previous LC)
dat$grp4 <- NA
dat$grp4[!(dat$grp2 %in% c("4_long_covid_any", "5_previous_long_covid_any"))] <- "1_without_lc"
dat$grp4[dat$grp2=="4_long_covid_any"] <- "2_current_lc"
dat$grp4[dat$grp2=="5_previous_long_covid_any"] <- "3_previous_lc"

### derive alternative exposure variable - LC (but not infection) status by time since infection
dat$exposure2 <- NA

dat$exposure2[dat$grp4=="1_without_lc"] <- "00_nolc"

dat$exposure2[dat$time_since_infection>=84 & dat$time_since_infection<126 & dat$grp4=="2_current_lc"] <- "01_12-18w_currlc"
dat$exposure2[dat$time_since_infection>=126 & dat$time_since_infection<168 & dat$grp4=="2_current_lc"] <- "02_18-24w_currlc"
dat$exposure2[dat$time_since_infection>=168 & dat$time_since_infection<210 & dat$grp4=="2_current_lc"] <- "03_24-30w_currlc"
dat$exposure2[dat$time_since_infection>=210 & dat$time_since_infection<280 & dat$grp4=="2_current_lc"] <- "04_30-40w_currlc"
dat$exposure2[dat$time_since_infection>=280 & dat$time_since_infection<364 & dat$grp4=="2_current_lc"] <- "05_40-52w_currlc"
dat$exposure2[dat$time_since_infection>=364 & dat$grp4=="2_current_lc"] <- "06_ge52w_currlc"

dat$exposure2[dat$time_since_infection>=84 & dat$time_since_infection<126 & dat$grp4=="3_previous_lc"] <- "07_12-18w_prevlc"
dat$exposure2[dat$time_since_infection>=126 & dat$time_since_infection<168 & dat$grp4=="3_previous_lc"] <- "08_18-24w_prevlc"
dat$exposure2[dat$time_since_infection>=168 & dat$time_since_infection<210 & dat$grp4=="3_previous_lc"] <- "09_24-30w_prevlc"
dat$exposure2[dat$time_since_infection>=210 & dat$time_since_infection<280 & dat$grp4=="3_previous_lc"] <- "10_30-40w_prevlc"
dat$exposure2[dat$time_since_infection>=280 & dat$time_since_infection<364 & dat$grp4=="3_previous_lc"] <- "11_40-52w_prevlc"
dat$exposure2[dat$time_since_infection>=364 & dat$grp4=="3_previous_lc"] <- "12_ge52w_prevlc"

### save dataset to working directory
save(dat, file=paste0(root, "\\dataset_filtered.RData"))

######################## OUTCOME = NOT LOOKING FOR WORK ########################

### main analysis
modelling.clogit(
  dset = dat,
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\_Main analysis"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### increasing interior knots in splines from 1 to 2
modelling.clogit(
  dset = dat,
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - two interior knots"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=3, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=3, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=3, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=3, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### increasing interior knots in splines from 1 to 3
modelling.clogit(
  dset = dat,
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - three interior knots"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=4, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=4, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=4, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=4, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### participants aged <50 only
modelling.clogit(
  dset = dat[dat$age_at_visit < 50,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - participants aged under 50 only"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### participants aged 50+ only
modelling.clogit(
  dset = dat[dat$age_at_visit >= 50,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - participants aged 50+ only"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### excluding retired participants
modelling.clogit(
  dset = dat[dat$retired==0,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - excluding retired participants"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### ever-infected study participants only
modelling.clogit(
  dset = dat[dat$ever_infected==1,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - ever-infected participants"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### pre-Omicron infections, unvaccinated when infected
modelling.clogit(
  dset = dat[dat$var_vacc %in% c("0_uninfected", "1_pre_omicron_unvaccinated"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - pre-Omicron, unvaccinated"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### pre-Omicron infections, vaccinated when infected
modelling.clogit(
  dset = dat[dat$var_vacc %in% c("0_uninfected", "2_pre_omicron_vaccinated"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - pre-Omicron, vaccinated"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### Omicron infections, vaccinated when infected
### note: follow-up restricted to <40 weeks post-infection (max possible for Omicron)
modelling.clogit.lt40w(
  dset = dat[dat$var_vacc %in% c("0_uninfected", "4_omicron_vaccinated") &
               (is.na(dat$time_since_infection) | dat$time_since_infection<280),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - Omicron, vaccinated"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### unvaccinated when infected
modelling.clogit(
  dset = dat[dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "1_unvaccinated"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - unvaccinated"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### vaccinated 14-90 days ago when infected
modelling.clogit(
  dset = dat[dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "2_14-90days"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - vaccinated 14-90 days ago"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### vaccinated 91-180 days ago when infected
modelling.clogit(
  dset = dat[dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "3_91-180days"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - vaccinated 91-180 days ago"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### vaccinated 181+ days ago when infected
### note: follow-up restricted to <52 weeks post-infection (max possible in 181+ group)
modelling.clogit.lt52w(
  dset = dat[dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "4_ge181days") &
               (is.na(dat$time_since_infection) | dat$time_since_infection<364),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - vaccinated 181+ days ago"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### face-to-face data collection only (i.e. excluding remote responses)
modelling.clogit(
  dset = dat[dat$remote_collection==0,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - face-to-face data collection only"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### excluding visits after reinfection
modelling.clogit(
  dset = dat[dat$reinfected==0,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - excluding visits after reinfection"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### controlling for time-varying rather than baseline health status
modelling.clogit(
  dset = dat,
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\CLR - controlling for time-varying health status"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### binary logistic regression, all eligible participants
modelling.logistic(
  dset = dat,
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\BLR - all participants"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0",
    "white_visit0",
    "gor9d_visit0",
    "imd_quintile_visit0",
    "health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### binary logistic regression, ever-infected participants only
modelling.logistic(
  dset = dat[dat$ever_infected==1,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = not looking\\BLR - ever-infected participants"),
  outcome = "notlooking",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0",
    "white_visit0",
    "gor9d_visit0",
    "imd_quintile_visit0",
    "health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

######################### OUTCOME = LONG-TERM ABSENCE #########################

### NOTE: sample restricted to visits when participants were in employment, from
### 1 Oct 2021 onwards (end of furlough)

### main analysis
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\_Main analysis"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### increasing interior knots in splines from 1 to 2
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - two interior knots"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=3, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=3, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=3, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=3, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### increasing interior knots in splines from 1 to 3
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - three interior knots"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=4, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=4, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=4, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=4, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### participants aged <50 only
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$age_at_visit < 50,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - participants aged under 50 only"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### participants aged 50+ only
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$age_at_visit >= 50,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - participants aged 50+ only"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### ever-infected study participants only
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$ever_infected==1,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - ever-infected participants"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### unvaccinated when infected
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "1_unvaccinated"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - unvaccinated"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### vaccinated 14-90 days ago when infected
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "2_14-90days"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - vaccinated 14-90 days ago"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### vaccinated 91-180 days ago when infected
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "3_91-180days"),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - vaccinated 91-180 days ago"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### vaccinated 181+ days ago when infected
### note: follow-up restricted to <52 weeks post-infection (max possible in 181+ group)
modelling.clogit.lt52w(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$time_since_last_dose_infection_cat %in% c("0_uninfected", "4_ge181days") &
               (is.na(dat$time_since_infection) | dat$time_since_infection<364),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - vaccinated 181+ days ago"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### controlling for time-varying rather than baseline health status
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - controlling for time-varying health status"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### face-to-face data collection only (i.e. excluding remote responses)
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$remote_collection==0,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - face-to-face data collection only"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### excluding visits after reinfection
modelling.clogit(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$reinfected==0,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\CLR - excluding visits after reinfection"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### binary logistic regression, all eligible participants
modelling.logistic(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")),],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\BLR - all participants"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0",
    "white_visit0",
    "gor9d_visit0",
    "imd_quintile_visit0",
    "health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)

### binary logistic regression, ever-infected participants only
modelling.logistic(
  dset = dat[dat$working==1 &
               dat$visit_date>=as.numeric(as.Date("2021-10-01")) &
               dat$ever_infected==1,],
  chart_cols = chart_colours,
  out_dir = paste0(root, "\\Outcome = off work\\BLR - ever-infected participants"),
  outcome = "offwork",
  covs = c(
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * ns(age_at_visit, df=2, Boundary.knots=quantile(age_at_visit, c(0.1, 0.9)))",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * sex_visit0",
    "ns(calendar_time_visit, df=2, Boundary.knots=quantile(calendar_time_visit, c(0.1, 0.9))) * health_status_visit0",
    "white_visit0",
    "gor9d_visit0",
    "imd_quintile_visit0",
    "health_status_visit0"
  ),
  run_mod1 = TRUE,
  run_mod2 = TRUE
)
