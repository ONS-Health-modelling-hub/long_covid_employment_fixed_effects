library(sqldf)
library(purrr)
library(reshape2)
library(ggplot2)

set.seed(230111)
n_iter = 10000
totals_file = "filepath\\Descriptive plots\\Rates by time since infection among LC\\estimates.csv"
mod_file = "filepath\\Modelling\\Outcome = not looking\\_Main analysis\\or_mod2.csv"
out_dir = "filepath"
out_name_ext = "age16-64"

#set.seed(230112)
#n_iter = 10000
#totals_file = "filepath\\Descriptive plots\\Rates by time since infection among LC\\estimates_50plus.csv"
#mod_file = "filepath\\Modelling\\Outcome = not looking\\Heterogeneous effects\\or_age2grp.csv"
#out_dir = "filepath"
#out_name_ext = "age50-64"

### read in estimates totals and modelled ORs
totals <- read.csv(totals_file)
mod <- read.csv(mod_file)

### drop totals for <12 weeks from infection
totals <- totals[totals$time_since_infection_cat!="00_lt12w",]
totals$time_since_infection_cat <- droplevels(totals$time_since_infection_cat)

### if input model doesn't include interactions:
### restrict ORs to those relating to current LC from adjusted model
if(out_name_ext=="age16-64") {
  mod <- mod[mod$model=="Adjusted",]
  mod <- mod[grep("current Long Covid", mod$exposure),]
}

### if input model includes interactions:
### restrict ORs to those relating to current LC for participants aged 50-64
if(out_name_ext=="age50-64") {
  mod <- mod[mod$int_var==1,]
  mod <- mod[grep("currlc", mod$main_var),]
}

### calculate coefficients and SE from ORs and CIs
mod$coeff <- log(mod$or)
mod$se_coeff <- (log(mod$ucl) - mod$coeff) / 1.96

### add time since infection variable to model coefficients
mod$time_since_infection_cat <- levels(totals$time_since_infection_cat)

### restrict model coefficients dataset to variables of interest
mod <- mod[,c("time_since_infection_cat", "coeff", "se_coeff")]

### join model coefficients onto totals dataset
totals <- merge(x=totals, y=mod, all.x=TRUE, all.y=FALSE, by="time_since_infection_cat")

### sort totals dataset by month and time since infection
totals <- totals[with(totals, order(month, time_since_infection_cat)),]
rownames(totals) <- NULL

### calculate ORs from model coefficients
totals$or <- exp(totals$coeff)

### calculate probability and odds of inactivity
totals$p_inactivity <- totals$y_wtd / totals$n_wtd
totals$odds_inactivity <- totals$p_inactivity / (1-totals$p_inactivity)

### calculate odds and probability of inactivity had participants not been
### infected with SARS-CoV-2 and developed LC (i.e. the counterfactual)
totals$odds_inactivity_counterfactual <- totals$odds_inactivity / totals$or
totals$p_inactivity_counterfactual <- totals$odds_inactivity_counterfactual / (1+totals$odds_inactivity_counterfactual)

### calculate total inactivity had participants not been infected with
### SARS-CoV-2 and developed LC (i.e. the counterfactual)
totals$y_wtd_counterfactual <- totals$n_wtd * totals$p_inactivity_counterfactual

### calculate inactivity "attributable to" LC
totals$inactivity_due_to_lc <- totals$y_wtd - totals$y_wtd_counterfactual

### sum observed and counterfactual inactivity across time-since-infection
### categories within months
results_df <- sqldf("
  select
    month,
    sum(y_wtd) as observed_est,
    sum(y_wtd_counterfactual) as counterfactual_est,
    sum(inactivity_due_to_lc) as attributable_est
  from totals
  group by month
")

############################# CONFIDENCE INTERVALS #############################

### set up empty lists to store simulated totals
observed_sim_list <- as.list(NULL)
counterfactual_sim_list <- as.list(NULL)
attributable_sim_list <- as.list(NULL)

### loop over iterations
for(i in 1:n_iter) {
  
  ### simulate inactivity totals
  y_sim <- unlist(map2(.x=totals$y_wtd, .y=totals$y_se, .f=function(x,y) rnorm(n=1, x, y)))
  
  ### simulate ORs
  or_sim <- exp(unlist(map2(.x=totals$coeff, .y=totals$se_coeff, .f=function(x,y) rnorm(n=1, x, y))))
  
  ### calculate probability and odds of inactivity
  p_inactivity_sim <- y_sim / totals$n_wtd
  odds_inactivity_sim <- p_inactivity_sim / (1-p_inactivity_sim)
  
  ### calculate odds and probability of inactivity had participants not been
  ### infected with SARS-CoV-2 and developed LC (i.e. the counterfactual)
  odds_inactivity_counterfactual_sim <- odds_inactivity_sim / or_sim
  p_inactivity_counterfactual_sim <- odds_inactivity_counterfactual_sim / (1+odds_inactivity_counterfactual_sim)
  
  ### calculate total inactivity had participants not been infected with
  ### SARS-CoV-2 and developed LC (i.e. the counterfactual)
  y_counterfactual_sim <- totals$n_wtd * p_inactivity_counterfactual_sim
  
  ### calculate inactivity "attributable to" LC
  inactivity_due_to_lc_sim <- y_sim - y_counterfactual_sim
  
  ### collate variables into data.frame
  totals_sim <- data.frame(
    time_since_infection_cat = totals$time_since_infection_cat,
    month = totals$month,
    y_sim = y_sim,
    y_counterfactual_sim = y_counterfactual_sim,
    inactivity_due_to_lc_sim = inactivity_due_to_lc_sim
  )
  
  ### sum observed and counterfactual inactivity across time-since-infection
  ### categories within months
  results_sim <- sqldf("
    select
      month,
      sum(y_sim) as observed,
      sum(y_counterfactual_sim) as counterfactual,
      sum(inactivity_due_to_lc_sim) as attributable
    from totals_sim
   group by month
  ")
  
  ### store observed, counterfactual and attributable totals in separate lists
  observed_sim_list[[i]] <- results_sim$observed
  counterfactual_sim_list[[i]] <- results_sim$counterfactual
  attributable_sim_list[[i]] <- results_sim$attributable
  
}

### coerce lists to data.frames
observed_sim_df <- as.data.frame(observed_sim_list)
counterfactual_sim_df <- as.data.frame(counterfactual_sim_list)
attributable_sim_df <- as.data.frame(attributable_sim_list)

### find standard deviations of the sampling distributions (i.e. estimated SEs)
observed_se <- apply(observed_sim_df, MARGIN=1, FUN=sd)
counterfactual_se <- apply(counterfactual_sim_df, MARGIN=1, FUN=sd)
attributable_se <- apply(attributable_sim_df, MARGIN=1, FUN=sd)

### calculate 95% confidence intervals around point estimates
results_df$observed_lcl <- results_df$observed_est - 1.96 * observed_se
results_df$counterfactual_lcl <- results_df$counterfactual_est - 1.96 * counterfactual_se
results_df$attributable_lcl <- results_df$attributable_est - 1.96 * attributable_se

results_df$observed_ucl <- results_df$observed_est + 1.96 * observed_se
results_df$counterfactual_ucl <- results_df$counterfactual_est + 1.96 * counterfactual_se
results_df$attributable_ucl <- results_df$attributable_est + 1.96 * attributable_se

### find 95% confidence intervals around point estimates as 2.5th and 97.5th
### percentiles of the empirical sampling distribution
results_df$observed_lclemp <- apply(observed_sim_df, MARGIN=1, FUN=function(x) quantile(x, 0.025))
results_df$counterfactual_lclemp <- apply(counterfactual_sim_df, MARGIN=1, FUN=function(x) quantile(x, 0.025))
results_df$attributable_lclemp <- apply(attributable_sim_df, MARGIN=1, FUN=function(x) quantile(x, 0.025))

results_df$observed_uclemp <- apply(observed_sim_df, MARGIN=1, FUN=function(x) quantile(x, 0.975))
results_df$counterfactual_uclemp <- apply(counterfactual_sim_df, MARGIN=1, FUN=function(x) quantile(x, 0.975))
results_df$attributable_uclemp <- apply(attributable_sim_df, MARGIN=1, FUN=function(x) quantile(x, 0.975))

### truncate lower confidence limits at 0
results_df$observed_lcl[results_df$observed_lcl < 0] <- 0
results_df$counterfactual_lcl[results_df$counterfactual_lcl < 0] <- 0
results_df$attributable_lcl[results_df$attributable_lcl < 0] <- 0

results_df$observed_lclemp[results_df$observed_lclemp < 0] <- 0
results_df$counterfactual_lclemp[results_df$counterfactual_lclemp < 0] <- 0
results_df$attributable_lclemp[results_df$attributable_lclemp < 0] <- 0

### round estimates to nearest integer
results_df <- round(results_df, 0)

### re-order columns
results_df <- results_df[,c(
  "month",
  "observed_est",
  "observed_lcl",
  "observed_ucl",
  "observed_lclemp",
  "observed_uclemp",
  "counterfactual_est",
  "counterfactual_lcl",
  "counterfactual_ucl",
  "counterfactual_lclemp",
  "counterfactual_uclemp",
  "attributable_est",
  "attributable_lcl",
  "attributable_ucl",
  "attributable_lclemp",
  "attributable_uclemp"
)]

### write out results to working directory
write.csv(results_df,
          file=paste0(out_dir, "\\estimates_", out_name_ext, ".csv", sep=""),
          row.names=FALSE)

##################################### PLOT #####################################

### convert dataset to long format
chart_data <- melt(results_df, id.vars="month")

### derive component and statistic indicators
chart_data$variable <- as.character(chart_data$variable)
chart_data$component <- t(as.data.frame(strsplit(chart_data$variable, "_")))[,1]
chart_data$statistic <- t(as.data.frame(strsplit(chart_data$variable, "_")))[,2]
chart_data$variable <- NULL

### give the counterfactual component the CIs relating to the total (so that the
### chart displays correctly)
observed_data <- chart_data[chart_data$component=="observed",]
observed_data$component <- NULL
colnames(observed_data)[colnames(observed_data)=="value"] <- "value_observed"

chart_data <- merge(
  x=chart_data,
  y=observed_data,
  on=c("month", "statistic"),
  all.x=TRUE,
  all.y=FALSE
)

chart_data$value <- ifelse(chart_data$component=="counterfactual" &
                             chart_data$statistic %in% c("lcl", "ucl", "lclemp", "uclemp"),
                           chart_data$value_observed, chart_data$value)

chart_data$value_observed <- NULL

### convert dataset to wide format
chart_data <- dcast(chart_data, month + component ~ statistic, value.var="value")

### drop overall totals, keeping the two components
chart_data <- chart_data[chart_data$component!="observed",]

### re-label components
chart_data$component <- as.factor(chart_data$component)
chart_data$component <- factor(chart_data$component, labels=c("Attributable", "Not attributable"))
chart_data$component <- relevel(chart_data$component, ref="Not attributable")

### produce plot - CIs calculated as y +/- 1.96 * SE
results_plot <- ggplot(chart_data, aes(x=as.factor(month), y=est/1000, fill=component)) +
  geom_bar(stat="identity", position="stack") +
  geom_errorbar(aes(ymin=lcl/1000, ymax=ucl/1000), size=0.3, width=0.2) +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)), limits=c(0,NA)) +
  scale_fill_brewer(palette="Blues") +
  ylab("Thousands") +
  theme(
    axis.title.x=element_blank(),
    axis.text.x = element_text(size=8, colour="black", face="plain", angle=90, vjust=0.5),
    axis.ticks.x=element_line(size=0.5, colour="black"),
    axis.line.x=element_blank(),
    axis.title.y=element_text(size=9, colour="black", face="bold"),
    axis.text.y = element_text(size=8, colour="black", face="plain"),
    axis.ticks.y=element_line(size=0.5, colour="black"),
    axis.line.y=element_blank(),
    panel.border=element_rect(size=0.5, colour="black", fill=NA),
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.spacing=unit(0.5, "cm"),
    strip.background=element_blank(),
    strip.text=element_text(size=9, colour="black", face="bold"),
    legend.title=element_blank(),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.justification="center",
    legend.text=element_text(size=9, colour="black", face="plain"),
    plot.margin=margin(0.2, 0.2, 0.2, 0.2, unit="cm")
  )

ggsave(filename=paste0(out_dir, "\\estimates_", out_name_ext, ".jpg"),
       plot=results_plot,
       width=14, height=12, units="cm")

### produce plot - CIs calculated as 2.5th and 97.5th percentiles of empirical distribution
results_plot <- ggplot(chart_data, aes(x=as.factor(month), y=est/1000, fill=component)) +
  geom_bar(stat="identity", position="stack") +
  geom_errorbar(aes(ymin=lclemp/1000, ymax=uclemp/1000), size=0.3, width=0.2) +
  scale_y_continuous(expand=expansion(mult=c(0,0.02)), limits=c(0,NA)) +
  scale_fill_brewer(palette="Blues") +
  ylab("Thousands") +
  theme(
    axis.title.x=element_blank(),
    axis.text.x = element_text(size=8, colour="black", face="plain", angle=90, vjust=0.5),
    axis.ticks.x=element_line(size=0.5, colour="black"),
    axis.line.x=element_blank(),
    axis.title.y=element_text(size=9, colour="black", face="bold"),
    axis.text.y = element_text(size=8, colour="black", face="plain"),
    axis.ticks.y=element_line(size=0.5, colour="black"),
    axis.line.y=element_blank(),
    panel.border=element_rect(size=0.5, colour="black", fill=NA),
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.spacing=unit(0.5, "cm"),
    strip.background=element_blank(),
    strip.text=element_text(size=9, colour="black", face="bold"),
    legend.title=element_blank(),
    legend.position="bottom",
    legend.direction="horizontal",
    legend.justification="center",
    legend.text=element_text(size=9, colour="black", face="plain"),
    plot.margin=margin(0.2, 0.2, 0.2, 0.2, unit="cm")
  )

ggsave(filename=paste0(out_dir, "\\estimates_", out_name_ext, "_empirical_ci.jpg"),
       plot=results_plot,
       width=14, height=12, units="cm")
