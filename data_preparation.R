library(haven)
library(sqldf)

out_dir = "filepath"
visit_dataset_date = "20221010"
vacc_dataset_date = "20221010"
cutoff_date = "2022-09-30"

############################### READ IN DATASETS ###############################

### read in visit-level data
input_file <- paste0("filepath/data_participant_clean.dta")

vars_of_interest <- c(
  "participant_id",
  "hh_id_fake",
  "visit_id",
  "visit_num",
  "visit_date",
  "dataset",
  "d_survey_mode_preference",
  "age_at_visit",
  "sex",
  "ethnicityg",
  "country",
  "gor9d",
  "cis20_samp",
  "imd_samp",
  "work_status_v1",
  "work_sector",
  "work_direct_contact_patients_etc",
  "gold_code",
  "health_conditions",
  "health_conditions_impact",
  "result_mk",
  "result_combined",
  "ct_mean",
  "ctpattern",
  "covid_test_swab_pos_first_date",
  "covid_test_blood_pos_first_date",
  "covid_date",
  "long_covid_have_symptoms",
  "reduce_activities_long_covid"
)

dat <- haven::read_dta(input_file, col_select=all_of(vars_of_interest))
dat <- zap_labels(dat)

### drop visits after cut-off date
dat <- dat[dat$visit_date <= cutoff_date,]

### sort dataset by participant ID, visit date, and visit ID
dat <- dat[with(dat, order(participant_id, visit_date, -xtfrm(visit_id))),]

### read in vaccination data
dat_vacc <- haven::read_dta(paste0("filepath/data_participant_vaccination.dta"),
                            col_select = c("participant_id",
                                           "covid_vaccine_date1",
                                           "covid_vaccine_date2",
                                           "covid_vaccine_date3",
                                           "covid_vaccine_date4",
                                           "covid_vaccine_type1",
                                           "covid_vaccine_type2",
                                           "covid_vaccine_type3",
                                           "covid_vaccine_type4"))

dat_vacc <- zap_labels(dat_vacc)

### join vaccination data
dat <- sqldf("
  select
    a.*, 
    b.covid_vaccine_date1, b.covid_vaccine_date2, b.covid_vaccine_date3, b.covid_vaccine_date4,
    b.covid_vaccine_type1, b.covid_vaccine_type2, b.covid_vaccine_type3, b.covid_vaccine_type4
    from dat as a
    left join dat_vacc as b
    on a.participant_id = b.participant_id
")

rm(dat_vacc); gc()

### convert all dates to numeric
dat$visit_date <- as.numeric(dat$visit_date)
dat$covid_date <- as.numeric(dat$covid_date)
dat$covid_test_swab_pos_first_date <- as.numeric(dat$covid_test_swab_pos_first_date)
dat$covid_test_blood_pos_first_date <- as.numeric(dat$covid_test_blood_pos_first_date)
dat$covid_vaccine_date1 <- as.numeric(dat$covid_vaccine_date1)
dat$covid_vaccine_date2 <- as.numeric(dat$covid_vaccine_date2)
dat$covid_vaccine_date3 <- as.numeric(dat$covid_vaccine_date3)
dat$covid_vaccine_date4 <- as.numeric(dat$covid_vaccine_date4)

############################## VACCINATION STATUS ##############################

### find last visit date
dat$last_visit_date <- ave(dat$visit_date, dat$participant_id, FUN=max)

### set vaccination dates to NA if later than last visit during follow-up period
dat$covid_vaccine_date1[!is.na(dat$covid_vaccine_date1) & dat$covid_vaccine_date1 > dat$last_visit_date] <- NA
dat$covid_vaccine_date2[!is.na(dat$covid_vaccine_date2) & dat$covid_vaccine_date2 > dat$last_visit_date] <- NA
dat$covid_vaccine_date3[!is.na(dat$covid_vaccine_date3) & dat$covid_vaccine_date3 > dat$last_visit_date] <- NA
dat$covid_vaccine_date4[!is.na(dat$covid_vaccine_date4) & dat$covid_vaccine_date4 > dat$last_visit_date] <- NA
dat$covid_vaccine_type1[is.na(dat$covid_vaccine_date1)] <- NA
dat$covid_vaccine_type2[is.na(dat$covid_vaccine_date2)] <- NA
dat$covid_vaccine_type3[is.na(dat$covid_vaccine_date3)] <- NA
dat$covid_vaccine_type4[is.na(dat$covid_vaccine_date4)] <- NA

### derive vaccination status at each visit
dat$vacc_status_visit <- "0"

dat$vacc_status_visit[!is.na(dat$covid_vaccine_date1) &
                        ((dat$visit_date - dat$covid_vaccine_date1) >= 14)] <- "1"

dat$vacc_status_visit[!is.na(dat$covid_vaccine_date2) &
                        ((dat$visit_date - dat$covid_vaccine_date2) >= 14)] <- "2"

dat$vacc_status_visit[!is.na(dat$covid_vaccine_date3) &
                        ((dat$visit_date - dat$covid_vaccine_date3) >= 14)] <- "3"

dat$vacc_status_visit[!is.na(dat$covid_vaccine_date4) &
                        ((dat$visit_date - dat$covid_vaccine_date4) >= 14)] <- "4"

### derive time since last dose at each visit
dat$time_since_last_dose <- ifelse(dat$vacc_status_visit=="0", NA,
                                   ifelse(dat$vacc_status_visit=="1", dat$visit_date - dat$covid_vaccine_date1,
                                          ifelse(dat$vacc_status_visit=="2", dat$visit_date - dat$covid_vaccine_date2,
                                                 ifelse(dat$vacc_status_visit=="3", dat$visit_date - dat$covid_vaccine_date3,
                                                        ifelse(dat$vacc_status_visit=="4", dat$visit_date - dat$covid_vaccine_date4, -999)))))

############################### EARLIEST DATES ###############################

### find date of earliest CIS visit
dat <- sqldf("
  select a.*, b.visit0_date
  from dat as a
  left join(
    select participant_id, min(visit_date) as visit0_date
    from dat
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
")


### find date of earliest positive PCR test during study
dat <- sqldf("
  select a.*, b.swab_study_date
  from dat as a
  left join(
    select participant_id, min(visit_date) as swab_study_date
    from dat
    where result_mk=1
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
")

### find earliest positive blood test during the study
dat <- sqldf("
  select a.*, b.blood_study_date
  from dat as a
  left join(
    select participant_id, min(visit_date) as blood_study_date
    from dat
    where result_combined=1
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
")

### find dates of earliest positive tests outside of study
dat <- sqldf("
  select a.*, b.swab_nonstudy_date, b.blood_nonstudy_date
  from dat as a
  left join(
    select
      participant_id,
      min(covid_test_swab_pos_first_date) as swab_nonstudy_date,
      min(covid_test_blood_pos_first_date) as blood_nonstudy_date
    from dat
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
")

### find date when participants first thought they had COVID
### Note: restricting this to visits from 15 Oct 2020, as evidence that a large
### proportion of participants change their mind re. when they first thought
### they had COVID on follow-up visits after this (i.e. people thought they had
### COVID in the first wave, then during the second wave, when knowledge of the
### virus had improved, they amended their first COVID date to later, or got rid
### of it altogether)
dat <- sqldf("
  select a.*, b.think_covid_date
  from dat as a
  left join(
    select
      participant_id,
      min(covid_date) as think_covid_date
    from dat
    where visit_date >= 18550
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
")

### find date of first response to LC question
dat <- sqldf("
  select a.*, b.lc_response_date
  from dat as a
  left join(
    select participant_id, min(visit_date) as lc_response_date
    from dat
    where long_covid_have_symptoms is not NULL
    group by participant_id
  ) as b
  on a.participant_id = b.participant_id
")

### set LC responses before 3 Feb 2021 (when the question was implemented) to NA
dat$lc_response_date[dat$lc_response_date < as.numeric(as.Date("2021-02-03"))] <- NA
dat$long_covid_have_symptoms[is.na(dat$lc_response_date)] <- NA

### set positive antibody results after first vaccination to NA
dat$blood_study_date[dat$blood_study_date >= dat$covid_vaccine_date1] <- NA
dat$blood_nonstudy_date[dat$blood_nonstudy_date >= dat$covid_vaccine_date1] <- NA

### COVID-19 arrived in the UK on 24 Jan 2020 - set any infection dates before this to NA
dat$swab_nonstudy_date[!is.na(dat$swab_nonstudy_date) &
                         (dat$swab_nonstudy_date < as.numeric(as.Date("2020-01-24")))] <- NA

dat$blood_nonstudy_date[!is.na(dat$blood_nonstudy_date) &
                          (dat$blood_nonstudy_date < as.numeric(as.Date("2020-01-24")))] <- NA

dat$think_covid_date[!is.na(dat$think_covid_date) &
                       (dat$think_covid_date < as.numeric(as.Date("2020-01-24")))] <- NA

### find infection date based on any suspected or confirmed infection
dat$infection_date <- pmin(dat$swab_study_date,
                           dat$swab_nonstudy_date,
                           dat$blood_study_date,
                           dat$blood_nonstudy_date,
                           dat$think_covid_date,
                           na.rm=TRUE)

### find date of first positive swab
dat$first_swab_date <- pmin(dat$swab_study_date,
                            dat$swab_nonstudy_date,
                            na.rm=TRUE)

### set infection dates to NA if later than last visit during follow-up period
dat$infection_date[dat$infection_date > dat$last_visit_date] <- NA
dat$first_swab_date[dat$first_swab_date > dat$last_visit_date] <- NA

############################## LONG COVID STATUS ##############################

### if no response to LC question then set activity limitation to NA
dat$reduce_activities_long_covid[is.na(dat$long_covid_have_symptoms)] <- NA

### if LC = no then set activity limitation to 0
dat$reduce_activities_long_covid[!is.na(dat$long_covid_have_symptoms) &
                                   dat$long_covid_have_symptoms==0] <- 0

### if LC = yes but activity limitation is missing then set activity limitation to none
dat$reduce_activities_long_covid[!is.na(dat$long_covid_have_symptoms) &
                                   dat$long_covid_have_symptoms==1 &
                                   is.na(dat$reduce_activities_long_covid)] <- 4

### derive flag for LC of any severity
dat$lc_any <- ifelse(!is.na(dat$long_covid_have_symptoms) &
                       dat$long_covid_have_symptoms==1, 1, 0)

### derive flag for activity-limiting LC
dat$lc_lim <- ifelse(!is.na(dat$reduce_activities_long_covid) &
                       dat$reduce_activities_long_covid %in% 5:6, 1, 0)

### set LC flags to 0 if before first positive swab or within first 12 weeks
dat$lc_any[is.na(dat$first_swab_date) |
             (!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) < 84)] <- 0

dat$lc_lim[is.na(dat$first_swab_date) |
             (!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) < 84)] <- 0

### flag all visits from first occurrence of LC=0 after at least one visit where LC=1
prev.lc <- function(lc_var) {
  lc_diff <- c(0, diff(lc_var))
  first0 <- which(lc_diff==-1)[1]
  out_var <- rep(0, length(lc_diff))
  if(!is.na(first0)) {out_var[first0:length(out_var)] <- 1}
  return(out_var)
}

dat$lc_any_prev <- ave(dat$lc_any, dat$participant_id, FUN=prev.lc)
dat$lc_lim_prev <- ave(dat$lc_lim, dat$participant_id, FUN=prev.lc)

############################## SOCIO-DEMOGRAPHICS ##############################

### define 10-year age-band variable
age_breaks <- c(2, 15, 24, 34, 49, 64, max(dat$age_at_visit, na.rm=TRUE))

dat$age10 <- cut(dat$age_at_visit,
                 breaks=age_breaks,
                 include.lowest=TRUE)

dat$age10 <- as.character(dat$age10)

### convert sex to character
dat$sex <- as.character(dat$sex)

### define white/non-white variable
dat$white <- as.character(ifelse(dat$ethnicityg==1, 1, 0))

### convert ethnic group to character
dat$ethnicityg <- as.character(dat$ethnicityg)

### define IMD quintile
dat$imd_quintile_eng <- ifelse(dat$imd_samp>(0*32844/5) &
                                 dat$imd_samp<=(1*32844/5), 1,
                               ifelse(dat$imd_samp>(1*32844/5) &
                                        dat$imd_samp<=(2*32844/5), 2,
                                      ifelse(dat$imd_samp>(2*32844/5) &
                                               dat$imd_samp<=(3*32844/5), 3,
                                             ifelse(dat$imd_samp>(3*32844/5) &
                                                      dat$imd_samp<=(4*32844/5), 4,
                                                    ifelse(dat$imd_samp>(4*32844/5) &
                                                             dat$imd_samp<=(5*32844/5), 5, -999)))))

dat$imd_quintile_wal <- ifelse(dat$imd_samp>(0*1909/5) &
                                 dat$imd_samp<=(1*1909/5), 1,
                               ifelse(dat$imd_samp>(1*1909/5) &
                                        dat$imd_samp<=(2*1909/5), 2,
                                      ifelse(dat$imd_samp>(2*1909/5) &
                                               dat$imd_samp<=(3*1909/5), 3,
                                             ifelse(dat$imd_samp>(3*1909/5) &
                                                      dat$imd_samp<=(4*1909/5), 4,
                                                    ifelse(dat$imd_samp>(4*1909/5) &
                                                             dat$imd_samp<=(5*1909/5), 5, -999)))))

dat$imd_quintile_sco <- ifelse(dat$imd_samp>(0*6976/5) &
                                 dat$imd_samp<=(1*6976/5), 1,
                               ifelse(dat$imd_samp>(1*6976/5) &
                                        dat$imd_samp<=(2*6976/5), 2,
                                      ifelse(dat$imd_samp>(2*6976/5) &
                                               dat$imd_samp<=(3*6976/5), 3,
                                             ifelse(dat$imd_samp>(3*6976/5) &
                                                      dat$imd_samp<=(4*6976/5), 4,
                                                    ifelse(dat$imd_samp>(4*6976/5) &
                                                             dat$imd_samp<=(5*6976/5), 5, -999)))))

dat$imd_quintile_ni <- ifelse(dat$imd_samp>(0*890/5) &
                                dat$imd_samp<=(1*890/5), 1,
                              ifelse(dat$imd_samp>(1*890/5) &
                                       dat$imd_samp<=(2*890/5), 2,
                                     ifelse(dat$imd_samp>(2*890/5) &
                                              dat$imd_samp<=(3*890/5), 3,
                                            ifelse(dat$imd_samp>(3*890/5) &
                                                     dat$imd_samp<=(4*890/5), 4,
                                                   ifelse(dat$imd_samp>(4*890/5) &
                                                            dat$imd_samp<=(5*890/5), 5, -999)))))

dat$imd_quintile <- ifelse(dat$country==0, dat$imd_quintile_eng,
                           ifelse(dat$country==1, dat$imd_quintile_wal,
                                  ifelse(dat$country==3, dat$imd_quintile_sco,
                                         ifelse(dat$country==2, dat$imd_quintile_ni, -999))))

dat$imd_quintile <- as.character(dat$imd_quintile)

### convert GOR to character
dat$gor9d <- as.character(dat$gor9d)

### convert CIS area to character
dat$cis20_samp <- as.character(dat$cis20_samp)

### convert country to character
dat$country <- as.character(dat$country)

### impute missing health conditions variable
dat$health_conditions[is.na(dat$health_conditions)] <- 0

### define health/disability status variable
dat$health_status <- 0
dat$health_status[dat$health_conditions==1 & is.na(dat$health_conditions_impact)] <- 1
dat$health_status[dat$health_conditions==1 & dat$health_conditions_impact==0] <- 1
dat$health_status[dat$health_conditions==1 & dat$health_conditions_impact==1] <- 2
dat$health_status[dat$health_conditions==1 & dat$health_conditions_impact==2] <- 3

### convert health conditions and health/disability status to character
dat$health_conditions <- as.character(dat$health_conditions)
dat$health_status <- as.character(dat$health_status)

### rename work status
dat$work_status <- dat$work_status_v1
dat$work_status_v1 <- NULL

### set work sector to 'unknown' if work sector is 'not working' but work status is 'working'
dat$work_sector[dat$work_sector==99 & dat$work_status %in% 1:4] <- NA

### set work sector to 'not working' if work status is 'not working'
dat$work_sector[dat$work_status %in% 5:6] <- 99

### set work sector to 'not applicable' if retired or in education
dat$work_sector <- ifelse(dat$work_status %in% 7:12, 999, dat$work_sector)

### only allow workers to be patient-facing if they work in health or social care
dat$healthsocialcare_pf <- ifelse(!is.na(dat$work_direct_contact_patients_etc) &
                                    dat$work_direct_contact_patients_etc==1 &
                                    !is.na(dat$work_sector) &
                                    dat$work_sector %in% 2:3, 1, 0)

### flag if self-employed
dat$self_employed <- ifelse(!is.na(dat$work_status) & dat$work_status %in% 3:4, 1, 0)

### derive SOC major group
dat$soc_major <- substr(dat$gold_code, 1, 1)

### set SOC major group to 'unknown' if SOC code is uncodeable or missing
dat$soc_major[dat$soc_major %in% c("u", "", "-")] <- "Unknown"

### set SOC major group to 'not working if work status is 'not working'
dat$soc_major[dat$work_status %in% 5:6] <- "Not_working"

### set SOC to 'not applicable' if retired or in education
dat$soc_major[dat$work_status %in% 7:12] <- "Not_applicable"

### convert employment variables to character
dat$work_status <- as.character(dat$work_status)
dat$work_sector <- as.character(dat$work_sector)
dat$healthsocialcare_pf <- as.character(dat$healthsocialcare_pf)
dat$self_employed <- as.character(dat$self_employed)

### get person-level variables at first CIS visit (i.e. enrolment)
dat <- sqldf("
  select
    a.*,
    b.age_at_visit as age_visit0,
    b.age10 as age10_visit0,
    b.sex as sex_visit0,
    b.ethnicityg as ethnicityg_visit0,
    b.white as white_visit0,
    b.imd_quintile as imd_quintile_visit0,
    b.gor9d as gor9d_visit0,
    b.cis20_samp as cis20_samp_visit0,
    b.country as country_visit0,
    b.health_conditions as health_conditions_visit0,
    b.health_status as health_status_visit0,
    b.work_status as work_status_visit0,
    b.work_sector as work_sector_visit0,
    b.healthsocialcare_pf as healthsocialcare_pf_visit0,
    b.self_employed as self_employed_visit0,
    b.soc_major as soc_major_visit0
  from dat as a
  left join(
    select *
    from dat
    where visit_date = visit0_date
  ) as b
  on a.participant_id = b.participant_id
")

### get person-level variables at first LC response
dat <- sqldf("
  select
    a.*,
    b.age_at_visit as age_lc0,
    b.age10 as age10_lc0,
    b.sex as sex_lc0,
    b.ethnicityg as ethnicityg_lc0,
    b.white as white_lc0,
    b.imd_quintile as imd_quintile_lc0,
    b.gor9d as gor9d_lc0,
    b.cis20_samp as cis20_samp_lc0,
    b.country as country_lc0,
    b.health_conditions as health_conditions_lc0,
    b.health_status as health_status_lc0,
    b.work_status as work_status_lc0,
    b.work_sector as work_sector_lc0,
    b.healthsocialcare_pf as healthsocialcare_pf_lc0,
    b.self_employed as self_employed_lc0,
    b.soc_major as soc_major_lc0
  from dat as a
  left join(
    select *
    from dat
    where visit_date = lc_response_date
  ) as b
  on a.participant_id = b.participant_id
")

############################### EMPLOYMENT FLAGS ###############################

### employment flags (visit-level)
dat$working <- ifelse(dat$work_status %in% c("1", "2", "3", "4"), 1, 0)
dat$offwork <- ifelse(dat$work_status %in% c("2", "4"), 1, 0)
dat$nonworking <- ifelse(dat$work_status %in% c("5", "6", "7",
                                                "8", "9", "10", "11", "12"), 1, 0)
dat$nonworking2 <- ifelse(dat$work_status %in% c("5", "6"), 1, 0)
dat$unemployed <- ifelse(dat$work_status %in% c("5"), 1, 0)
dat$notlooking <- ifelse(dat$work_status %in% c("6"), 1, 0)
dat$retired <- ifelse(dat$work_status %in% c("7"), 1, 0)
dat$student <- ifelse(dat$work_status %in% c("8", "9", "10", "11", "12"), 1, 0)
dat$inactive <- ifelse(dat$work_status %in% c("6", "7",
                                              "8", "9", "10", "11", "12"), 1, 0)

############################### EXPOSURE GROUPS ###############################

### define exposure group setup 1 - uninfected, previous infection
dat$grp1 <- "1_uninfected"
dat$grp1[!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) >= 0 & (dat$visit_date - dat$first_swab_date) < 84] <- "2_previous_infection_lt12w"
dat$grp1[!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) >= 84] <- "3_previous_infection_ge12w"

### define exposure group setup 2 - uninfected, previous infection,
### current LC of any severity, previous LC of any severity
dat$grp2 <- "1_uninfected"
dat$grp2[!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) >= 0 & (dat$visit_date - dat$first_swab_date) < 84] <- "2_previous_infection_lt12w"
dat$grp2[!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) >= 84] <- "3_previous_infection_ge12w"
dat$grp2[dat$lc_any==1] <- "4_long_covid_any"
dat$grp2[dat$lc_any==0 & dat$lc_any_prev==1] <- "5_previous_long_covid_any"

### define exposure group setup 3 - as above, but with activity-limiting LC
dat$grp3 <- "1_uninfected"
dat$grp3[!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) >= 0 & (dat$visit_date - dat$first_swab_date) < 84] <- "2_previous_infection_lt12w"
dat$grp3[!is.na(dat$first_swab_date) & (dat$visit_date - dat$first_swab_date) >= 84] <- "3_previous_infection_ge12w"
dat$grp3[dat$lc_lim==1] <- "4_long_covid_lim"
dat$grp3[dat$lc_lim==0 & dat$lc_lim_prev==1] <- "5_previous_long_covid_lim"

############################ PERSON-LEVEL VARIABLES ############################

### flag if ever received positive swab
dat$ever_infected <- as.numeric(!is.na(dat$first_swab_date))

### flag if ever had LC
dat$ever_lc_any <- ave(dat$lc_any, dat$participant_id, FUN=max)
dat$ever_lc_lim <- ave(dat$lc_lim, dat$participant_id, FUN=max)

### derive calendar time of first positive swab (days since 24 Jan 2020)
dat$calendar_time_infection <- dat$first_swab_date - as.numeric(as.Date("2020-01-24"))

### derive month and year of first positive swab
dat$yyyymm_infection <- format(as.Date(dat$first_swab_date, origin="1970-01-01"), "%Y-%m")

### derive month of first positive swab
dat$month_infection <- as.character(substr(dat$yyyymm_infection, 6, 7))

### derive time period of first positive swab
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

### derive vaccination status at first positive swab
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

### derive alternative vaccination statuses at first positive swab
dat$vacc_status_infection2 <- ifelse(dat$vacc_status_infection %in% c("1","2","3","4"), "1+",
                                     dat$vacc_status_infection)

dat$vacc_status_infection3 <- ifelse(dat$vacc_status_infection %in% c("2","3","4"), "2+",
                                     dat$vacc_status_infection)

############################ VISIT-LEVEL VARIABLES ############################

### derive calendar time of each visit (days since 24 Jan 2020)
dat$calendar_time_visit <- dat$visit_date - as.numeric(as.Date("2020-01-24"))

### derive month and year of each visit
dat$yyyymm_visit <- format(as.Date(dat$visit_date, origin="1970-01-01"), "%Y-%m")

### derive month of each visit
dat$month_visit <- as.character(substr(dat$yyyymm_visit, 6, 7))

### derive time since first positive swab at each visit
dat$time_since_infection <- dat$visit_date - dat$first_swab_date

### flag visits after first positive swab
dat$post_infection <- ifelse(!is.na(dat$time_since_infection) & dat$time_since_infection >= 0, 1, 0)

### derive variables relating to data collection mode
dat$remote_collection <- ifelse(dat$dataset==3, 1, 0)

dat$collection_mode <- NA
dat$collection_mode[is.na(dat$d_survey_mode_preference)] <- "1_studyworker"
dat$collection_mode[!is.na(dat$d_survey_mode_preference) &
                      dat$d_survey_mode_preference==0] <- "2_online"
dat$collection_mode[!is.na(dat$d_survey_mode_preference) &
                      dat$d_survey_mode_preference==1] <- "3_telephone"

############################ IDENTIFY REINFECTIONS ############################

### flag negative PCR tests after first positive swab
dat$neg_pcr <- ifelse(!is.na(dat$result_mk) & dat$result_mk==0, 1, 0)
dat$neg_after_pos <- ifelse(dat$post_infection==1, dat$neg_pcr, 0)

### count number of successive negatives; counter restarts after each positive
successive.negs <- function(x) {return(ifelse(x==0, 0, sequence(rle(x)$lengths)))}
dat$successive_negs <- ave(dat$neg_after_pos, dat$participant_id, FUN=successive.negs)

### convert to max number of previous successive negatives up until each visit
dat$successive_negs <- ave(dat$successive_negs, dat$participant_id, FUN=cummax)

### identify reinfections using CIS PCR tests (rules taken from CIS publication)
dat$reinfection <- 0

dat$reinfection[!is.na(dat$result_mk) & dat$result_mk==1 &
                  !is.na(dat$time_since_infection) & dat$time_since_infection>=120 &
                  dat$successive_negs>=1] <- 1

dat$reinfection[!is.na(dat$result_mk) & dat$result_mk==1 &
                  !is.na(dat$time_since_infection) & dat$time_since_infection>=90 &
                  dat$successive_negs>=2] <- 1

dat$reinfection[dat$visit_date >= as.numeric(as.Date("2021-12-20")) &
                  !is.na(dat$result_mk) & dat$result_mk==1 &
                  !is.na(dat$time_since_infection) & dat$time_since_infection>=90 &
                  dat$successive_negs>=1] <- 1

dat$reinfection[!is.na(dat$result_mk) & dat$result_mk==1 &
                  !is.na(dat$time_since_infection) & dat$time_since_infection>=60 &
                  dat$successive_negs>=3] <- 1

dat$reinfection[!is.na(dat$result_mk) & dat$result_mk==1 &
                  !is.na(dat$time_since_infection) & dat$time_since_infection>0 &
                  dat$successive_negs>=4] <- 1

### also identify reinfections using positive swabs outside of the CIS
dat$reinfection[!is.na(dat$covid_test_swab_pos_first_date) &
                  dat$covid_test_swab_pos_first_date <= dat$visit_date &
                  !is.na(dat$first_swab_date) &
                  dat$covid_test_swab_pos_first_date - dat$first_swab_date >= 120 &
                  dat$successive_negs>=1] <- 1

dat$reinfection[!is.na(dat$covid_test_swab_pos_first_date) &
                  dat$covid_test_swab_pos_first_date <= dat$visit_date &
                  !is.na(dat$first_swab_date) &
                  dat$covid_test_swab_pos_first_date - dat$first_swab_date >= 90 &
                  dat$successive_negs>=2] <- 1

dat$reinfection[dat$visit_date >= as.numeric(as.Date("2021-12-20")) &
                  !is.na(dat$covid_test_swab_pos_first_date) &
                  dat$covid_test_swab_pos_first_date <= dat$visit_date &
                  !is.na(dat$first_swab_date) &
                  dat$covid_test_swab_pos_first_date - dat$first_swab_date >= 90 &
                  dat$successive_negs>=1] <- 1

dat$reinfection[!is.na(dat$covid_test_swab_pos_first_date) &
                  dat$covid_test_swab_pos_first_date <= dat$visit_date &
                  !is.na(dat$first_swab_date) &
                  dat$covid_test_swab_pos_first_date - dat$first_swab_date >= 60 &
                  dat$successive_negs>=3] <- 1

dat$reinfection[!is.na(dat$covid_test_swab_pos_first_date) &
                  dat$covid_test_swab_pos_first_date <= dat$visit_date &
                  !is.na(dat$first_swab_date) &
                  dat$covid_test_swab_pos_first_date - dat$first_swab_date > 0 &
                  dat$successive_negs>=4] <- 1

### finally, carry forward reinfection flag from first visit after reinfection
dat$reinfected <- ave(dat$reinfection, dat$participant_id, FUN=cummax)

############################## WRITE OUT DATASETS ##############################

### save dataset to working directory
save(dat, file=paste0(out_dir, "\\dataset.RData"))
