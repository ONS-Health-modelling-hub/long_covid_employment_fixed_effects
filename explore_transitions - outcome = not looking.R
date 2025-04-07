library(sqldf)

infile = "filepath\\dataset_filtered.RData"
out_dir = "filepath"

################################## DATA PREP ##################################

### read in dataset
load(infile)

### derive pre vs. post Long Covid flag
dat$post_lc <- ifelse(dat$lc_any==1 | dat$lc_any_prev==1, 1, 0)

### re-derive ever Long Covid flag
dat$ever_lc_any <- ave(dat$lc_any, dat$participant_id, FUN=max)

### identify each participant's first visit during follow-up
dat$first_visit_date <- ave(dat$visit_date, dat$participant_id, FUN=min)
dat$first_visit <- ifelse(dat$visit_date==dat$first_visit_date, 1, 0)

### identify if each participant was working or inactive at their previous visit
dat$prev_participant_id <- c("X", dat$participant_id[1:(nrow(dat)-1)])

dat$prev_working <- c(0, dat$working[1:(nrow(dat)-1)])
dat$prev_working <- dat$prev_working * (dat$participant_id==dat$prev_participant_id)

dat$prev_notlooking <- c(0, dat$notlooking[1:(nrow(dat)-1)])
dat$prev_notlooking <- dat$prev_notlooking * (dat$participant_id==dat$prev_participant_id)

### flag transitions to inactive
dat$transition_inactive <- 0
dat$transition_inactive[dat$first_visit==0 &
                          dat$prev_notlooking==0 &
                          dat$notlooking==1] <- 1

### flag transitions to working
dat$transition_working <- 0
dat$transition_working[dat$first_visit==0 &
                         dat$prev_working==0 &
                         dat$working==1] <- 1

### flag if participants ever transitioned to inactive or working
dat$ever_transition_inactive <- ave(dat$transition_inactive, dat$participant_id, FUN=max)
dat$ever_transition_working <- ave(dat$transition_working, dat$participant_id, FUN=max)

### flag instances of participants transitioning to inactivity whilst having Long Covid
dat$transition_inactive_lc <- dat$transition_inactive * dat$lc_any
dat$ever_transition_inactive_lc <- ave(dat$transition_inactive_lc, dat$participant_id, FUN=max)

### flag instances of participants transitioning to working whilst having Long Covid
### or previously having Long Covid
dat$transition_working_lc <- ifelse(dat$transition_working==1 &
                                      dat$grp2 %in% c("4_long_covid_any", "5_previous_long_covid_any"), 1, 0)
dat$ever_transition_working_lc <- ave(dat$transition_working_lc, dat$participant_id, FUN=max)

### find date when participants first transitioned to inactivity whilst having Long Covid
dat <- sqldf("
  select a.*, b.transition_inactive_lc_date
  from dat as a
  left join(
    select participant_id, min(visit_date) as transition_inactive_lc_date
    from dat
    where transition_inactive_lc=1
    group by participant_id
  ) as b
  on a.participant_id=b.participant_id
")

### flag instances of participants transitioning back to working having previously
### transitioned to inactive whilst having Long Covid
dat$transition_back_working_lc <- ifelse(dat$ever_transition_inactive_lc==1 &
                                           dat$transition_working_lc==1 &
                                           dat$visit_date >= dat$transition_inactive_lc_date, 1, 0)
dat$ever_transition_back_working_lc <- ave(dat$transition_back_working_lc, dat$participant_id, FUN=max)

### find date when participants first transitioning back to working having previously
### transitioned to inactive whilst having Long Covid
dat <- sqldf("
  select a.*, b.transition_back_working_lc_date
  from dat as a
  left join(
    select participant_id, min(visit_date) as transition_back_working_lc_date
    from dat
    where transition_back_working_lc=1
    group by participant_id
  ) as b
  on a.participant_id=b.participant_id
")

### subset dataset to visits when participants were working at the previous visit
### (and removing first visits during follow-up)
#dat_working <- dat[dat$first_visit==0 & dat$prev_working==1,]

### subset dataset to visits when participants were inactive at the previous visit
### (and removing first visits during follow-up)
#dat_inactive <- dat[dat$first_visit==0 & dat$prev_notlooking==1,]

########################### SUMMARIES: BACK TO WORK ###########################

### find percentage of who participants who transitioned back to working having
### previously transitioned to inactive whilst having Long Covid
transitions_back_to_working_lc <- sqldf("
  select
    count(*) as n_participants,
    sum(transition) as n_transitioned,
    sum(transition) / count(*) * 100 as pct_transitioned
  from(
    select participant_id, max(ever_transition_back_working_lc) as transition
    from dat
    where ever_transition_inactive_lc=1
    group by participant_id
  )
")

write.csv(transitions_back_to_working_lc,
          file=paste0(out_dir, "\\transitions_back_to_working_lc.csv"),
          row.names=FALSE)

############################ SUMMARIES: INACTIVITY ############################

### find percentage of who participants who were ever inactive according to whether
### they ever had Long Covid
ever_inactive_by_ever_lc <- sqldf("
  select
    ever_lc_any,
    count(*) as n_participants,
    sum(notlooking) as n_notlooking,
    sum(notlooking) / count(*) * 100 as pct_notlooking
  from(
    select participant_id, ever_lc_any, max(notlooking) as notlooking
    from dat
    group by participant_id, ever_lc_any
  )
  group by ever_lc_any
")

write.csv(ever_inactive_by_ever_lc,
          file=paste0(out_dir, "\\ever_inactive_by_ever_lc.csv"),
          row.names=FALSE)

### find percentage of who participants who were ever inactive while in each exposure state
ever_inactive <- sqldf("
  select
    grp2,
    count(*) as n_participants,
    sum(notlooking) as n_notlooking,
    sum(notlooking) / count(*) * 100 as pct_notlooking
  from(
    select participant_id, grp2, max(notlooking) as notlooking
    from dat
    group by participant_id, grp2
  )
  group by grp2
")

write.csv(ever_inactive,
          file=paste0(out_dir, "\\ever_inactive.csv"),
          row.names=FALSE)

### find percentage of visits when participants were inactive while in each exposure state
visits_inactive <- sqldf("
  select
    grp2,
    count(*) as n_visits,
    sum(notlooking) as n_notlooking,
    sum(notlooking) / count(*) * 100 as pct_notlooking
  from dat
  group by grp2
")

write.csv(visits_inactive,
          file=paste0(out_dir, "\\visits_inactive.csv"),
          row.names=FALSE)

### find percentage of who participants who were ever inactive before vs. after Long Covid
ever_inactive_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_participants,
    sum(notlooking) as n_notlooking,
    sum(notlooking) / count(*) * 100 as pct_notlooking
  from(
    select participant_id, post_lc, max(notlooking) as notlooking
    from dat
    group by participant_id, post_lc
  )
  group by post_lc
")

write.csv(ever_inactive_coarse,
          file=paste0(out_dir, "\\ever_inactive_coarse.csv"),
          row.names=FALSE)

### find percentage of visits when participants were inactive before vs. after Long Covid
visits_inactive_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_visits,
    sum(notlooking) as n_notlooking,
    sum(notlooking) / count(*) * 100 as pct_notlooking
  from dat
  group by post_lc
")

write.csv(visits_inactive_coarse,
          file=paste0(out_dir, "\\visits_inactive_coarse.csv"),
          row.names=FALSE)

### find percentage of who participants who ever transitioned to inactivity
### according to whether they ever had Long Covid
ever_transition_to_inactive_by_ever_lc <- sqldf("
  select
    ever_lc_any,
    count(*) as n_participants,
    sum(ever_transition) as n_ever_transition,
    sum(ever_transition) / count(*) * 100 as pct_ever_transition
  from(
    select participant_id, ever_lc_any, max(ever_transition_inactive) as ever_transition
    from dat
    group by participant_id, ever_lc_any
  )
  group by ever_lc_any
")

write.csv(ever_transition_to_inactive_by_ever_lc,
          file=paste0(out_dir, "\\ever_transition_to_inactive_by_ever_lc.csv"),
          row.names=FALSE)

### find percentage of who participants who transitioned to inactivity
### while in each exposure state
transitions_to_inactivity <- sqldf("
  select
    grp2,
    count(*) as n_participants,
    sum(transition) as n_transitioned,
    sum(transition) / count(*) * 100 as pct_transitioned
  from(
    select participant_id, grp2, max(transition_inactive) as transition
    from dat
    group by participant_id, grp2
  )
  group by grp2
")

write.csv(transitions_to_inactivity,
          file=paste0(out_dir, "\\transitions_to_inactivity.csv"),
          row.names=FALSE)

### find percentage of who participants who transitioned to inactivity
### before vs. after Long Covid
transitions_to_inactivity_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_participants,
    sum(transition) as n_transitioned,
    sum(transition) / count(*) * 100 as pct_transitioned
  from(
    select participant_id, post_lc, max(transition_inactive) as transition
    from dat
    group by participant_id, post_lc
  )
  group by post_lc
")

write.csv(transitions_to_inactivity_coarse,
          file=paste0(out_dir, "\\transitions_to_inactivity_coarse.csv"),
          row.names=FALSE)

############################## SUMMARIES: WORKING ##############################

### find percentage of who participants who were ever working according to whether
### they ever had Long Covid
ever_working_by_ever_lc <- sqldf("
  select
    ever_lc_any,
    count(*) as n_participants,
    sum(working) as n_working,
    sum(working) / count(*) * 100 as pct_working
  from(
    select participant_id, ever_lc_any, max(working) as working
    from dat
    group by participant_id, ever_lc_any
  )
  group by ever_lc_any
")

write.csv(ever_working_by_ever_lc,
          file=paste0(out_dir, "\\ever_working_by_ever_lc.csv"),
          row.names=FALSE)

### find percentage of who participants who were ever working while in each exposure state
ever_working <- sqldf("
  select
    grp2,
    count(*) as n_participants,
    sum(working) as n_working,
    sum(working) / count(*) * 100 as pct_working
  from(
    select participant_id, grp2, max(working) as working
    from dat
    group by participant_id, grp2
  )
  group by grp2
")

write.csv(ever_working,
          file=paste0(out_dir, "\\ever_working.csv"),
          row.names=FALSE)

### find percentage of visits when participants were working while in each exposure state
visits_working <- sqldf("
  select
    grp2,
    count(*) as n_visits,
    sum(working) as n_working,
    sum(working) / count(*) * 100 as pct_working
  from dat
  group by grp2
")

write.csv(visits_working,
          file=paste0(out_dir, "\\visits_working.csv"),
          row.names=FALSE)

### find percentage of who participants who were ever working before vs. after Long Covid
ever_working_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_participants,
    sum(working) as n_working,
    sum(working) / count(*) * 100 as pct_working
  from(
    select participant_id, post_lc, max(working) as working
    from dat
    group by participant_id, post_lc
  )
  group by post_lc
")

write.csv(ever_working_coarse,
          file=paste0(out_dir, "\\ever_working_coarse.csv"),
          row.names=FALSE)

### find percentage of visits when participants were working before vs. after Long Covid
visits_working_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_visits,
    sum(working) as n_working,
    sum(working) / count(*) * 100 as pct_working
  from dat
  group by post_lc
")

write.csv(visits_working_coarse,
          file=paste0(out_dir, "\\visits_working_coarse.csv"),
          row.names=FALSE)

### find percentage of who participants who ever transitioned to working
### according to whether they ever had Long Covid
ever_transition_to_working_by_ever_lc <- sqldf("
  select
    ever_lc_any,
    count(*) as n_participants,
    sum(ever_transition) as n_ever_transition,
    sum(ever_transition) / count(*) * 100 as pct_ever_transition
  from(
    select participant_id, ever_lc_any, max(ever_transition_working) as ever_transition
    from dat
    group by participant_id, ever_lc_any
  )
  group by ever_lc_any
")

write.csv(ever_transition_to_working_by_ever_lc,
          file=paste0(out_dir, "\\ever_transition_to_working_by_ever_lc.csv"),
          row.names=FALSE)

### find percentage of who participants who transitioned to working
### while in each exposure state
transitions_to_working <- sqldf("
  select
    grp2,
    count(*) as n_participants,
    sum(transition) as n_transitioned,
    sum(transition) / count(*) * 100 as pct_transitioned
  from(
    select participant_id, grp2, max(transition_working) as transition
    from dat
    group by participant_id, grp2
  )
  group by grp2
")

write.csv(transitions_to_working,
          file=paste0(out_dir, "\\transitions_to_working.csv"),
          row.names=FALSE)

### find percentage of who participants who transitioned to working
### before vs. after Long Covid
transitions_to_working_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_participants,
    sum(transition) as n_transitioned,
    sum(transition) / count(*) * 100 as pct_transitioned
  from(
    select participant_id, post_lc, max(transition_working) as transition
    from dat
    group by participant_id, post_lc
  )
  group by post_lc
")

write.csv(transitions_to_working_coarse,
          file=paste0(out_dir, "\\transitions_to_working_coarse.csv"),
          row.names=FALSE)

############################## SUMMARIES: RETIRED ##############################

### find percentage of visits when participants were retired while in each exposure state
visits_retired <- sqldf("
  select
    grp2,
    count(*) as n_visits,
    sum(retired) as n_retired,
    sum(retired) / count(*) * 100 as pct_retired
  from dat
  group by grp2
")

write.csv(visits_retired,
          file=paste0(out_dir, "\\visits_retired.csv"),
          row.names=FALSE)

### find percentage of visits when participants were retired before vs. after Long Covid
visits_retired_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_visits,
    sum(retired) as n_retired,
    sum(retired) / count(*) * 100 as pct_retired
  from dat
  group by post_lc
")

write.csv(visits_retired_coarse,
          file=paste0(out_dir, "\\visits_retired_coarse.csv"),
          row.names=FALSE)

############################ SUMMARIES: UNEMPLOYED ############################

### find percentage of visits when participants were unemployed while in each exposure state
visits_unemployed <- sqldf("
  select
    grp2,
    count(*) as n_visits,
    sum(unemployed) as n_unemployed,
    sum(unemployed) / count(*) * 100 as pct_unemployed
  from dat
  group by grp2
")

write.csv(visits_unemployed,
          file=paste0(out_dir, "\\visits_unemployed.csv"),
          row.names=FALSE)

### find percentage of visits when participants were unemployed before vs. after Long Covid
visits_unemployed_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_visits,
    sum(unemployed) as n_unemployed,
    sum(unemployed) / count(*) * 100 as pct_unemployed
  from dat
  group by post_lc
")

write.csv(visits_unemployed_coarse,
          file=paste0(out_dir, "\\visits_unemployed_coarse.csv"),
          row.names=FALSE)

######################### SUMMARIES: LONG-TERM ABSENCE #########################

### find percentage of who participants who were ever absent according to whether
### they ever had Long Covid (whilst in work, visits from 1 Oct 2021 only)
ever_offwork_by_ever_lc <- sqldf("
  select
    ever_lc_any,
    count(*) as n_participants,
    sum(offwork) as n_offwork,
    sum(offwork) / count(*) * 100 as pct_offwork
  from(
    select participant_id, ever_lc_any, max(offwork) as offwork
    from dat
    where working=1 and visit_date>=18901
    group by participant_id, ever_lc_any
  )
  group by ever_lc_any
")

write.csv(ever_offwork_by_ever_lc,
          file=paste0(out_dir, "\\ever_offwork_by_ever_lc.csv"),
          row.names=FALSE)

### find percentage of visits when participants were absent while in each exposure state
visits_offwork <- sqldf("
  select
    grp2,
    count(*) as n_visits,
    sum(offwork) as n_offwork,
    sum(offwork) / count(*) * 100 as pct_offwork
  from dat
  where working=1 and visit_date>=18901
  group by grp2
")

write.csv(visits_offwork,
          file=paste0(out_dir, "\\visits_offwork.csv"),
          row.names=FALSE)

### find percentage of visits when participants were absent before vs. after Long Covid
visits_offwork_coarse <- sqldf("
  select
    post_lc,
    count(*) as n_visits,
    sum(offwork) as n_offwork,
    sum(offwork) / count(*) * 100 as pct_offwork
  from dat
  where working=1 and visit_date>=18901
  group by post_lc
")

write.csv(visits_offwork_coarse,
          file=paste0(out_dir, "\\visits_offwork_coarse.csv"),
          row.names=FALSE)
