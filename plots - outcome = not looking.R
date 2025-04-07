library(ggplot2)

### list input directories containing OR data
in_dirs = c(
  "filepath\\Heterogeneous effects",
  "filepath\\Heterogeneous effects\\Time since last vaccine dose",
  "filepath\\Heterogeneous effects\\Time since last vaccine dose v2",
  "filepath\\Heterogeneous effects\\Variant period",
  "filepath\\Heterogeneous effects\\Variant x time since last vaccine dose"
)

### define output directory to save plots
out_dir = "filepath"

### loop over input directories
for(j in 1:length(in_dirs)) {
  
  ### only carry on with this iteration if the directory contains OR datasets
  if(length(list.files(in_dirs[j], pattern="or_")) > 0) {
    
    ### list all OR datasets in input directory
    file_list <- paste0(in_dirs[j], "\\", list.files(in_dirs[j], pattern="or_"))
    file_list <- file_list[grep(".csv", file_list)]
    if(length(grep("coeffs_", file_list))>0) {file_list <- file_list[-grep("coeffs_", file_list)]}
    if(length(grep("vcov_", file_list))>0) {file_list <- file_list[-grep("vcov_", file_list)]}
    if(length(grep("lr_", file_list))>0) {file_list <- file_list[-grep("lr_", file_list)]}
    
    ### list all interaction variables
    int_names <- t(as.data.frame(strsplit(file_list, "\\\\")))
    int_names <- int_names[,ncol(int_names)]
    int_names <- gsub(".csv", "", int_names)
    int_names <- substr(int_names, 4, nchar(int_names))
    int_names <- unique(int_names)
    
    ### loop over interaction variables
    for(i in 1:length(file_list)) {
      
      ### read in odds ratios
      dat <- read.csv(file_list[i])
      
      ### restrict dataset to ORs relating to current Long Covid
      dat <- dat[grep("_currlc", dat$main_var),]
      
      ### convert interaction variable to factor
      dat$int_var <- as.factor(dat$int_var)
      
      ### assign labels to interaction variable
      if(int_names[i]=="age2grp") {dat$int_var <- factor(dat$int_var, labels=c("16-49 years", "50-64 years"))}
      if(int_names[i]=="age3grp") {dat$int_var <- factor(dat$int_var, labels=c("16-34 years", "35-49 years", "50-64 years"))}
      if(int_names[i]=="age4grp") {dat$int_var <- factor(dat$int_var, labels=c("16-24 years", "25-34 years", "35-49 years", "50-64 years"))}
      if(int_names[i]=="health") {dat$int_var <- factor(dat$int_var, labels=c("No health conditions", "Not limited", "Limited a little", "Limited a lot"))}
      if(int_names[i]=="imd") {dat$int_var <- factor(dat$int_var, labels=c("1 (most deprived)", "2", "3", "4", "5 (least deprived)"))}
      if(int_names[i]=="mode") {dat$int_var <- factor(dat$int_var, labels=c("Face-to-face", "Remote"))}
      if(int_names[i]=="reinfected") {dat$int_var <- factor(dat$int_var, labels=c("Not reinfected", "Reinfected"))}
      if(int_names[i]=="sex") {dat$int_var <- factor(dat$int_var, labels=c("Male", "Female"))}
      if(int_names[i]=="white") {dat$int_var <- factor(dat$int_var, labels=c("Non-white", "White"))}
      if(int_names[i]=="variant") {dat$int_var <- factor(dat$int_var, labels=c("Pre-Omicron", "Omicron"))}
      if(int_names[i]=="time_since_last_dose") {dat$int_var <- factor(dat$int_var,
                                                                      labels=c("Unvaccinated", "14-90 days",
                                                                               "91-180 days", ">180 days"))}
      if(int_names[i]=="time_since_last_dose2") {dat$int_var <- factor(dat$int_var,
                                                                       labels=c("Unvaccinated", "First dose",
                                                                                "14-90 days after second dose",
                                                                                "91-180 days after second dose",
                                                                                ">180 days after second dose",
                                                                                "14-90 days after third/fourth dose",
                                                                                ">90 days after third/fourth dose"))}
      if(int_names[i]=="variant_by_time") {dat$int_var <- factor(dat$int_var,
                                                                 labels=c("Unvaccinated", "Pre-Omicron, 14-90 days",
                                                                          "Pre-Omicron, 91-180 days", "Pre-Omicron, >180 days",
                                                                          "Omicron, 14-90 days", "Omicron, 91-180 days",
                                                                          "Omicron, >180 days"))}
      if(int_names[i]=="sector") {dat$int_var <- factor(dat$int_var,
                                                        labels=c("Education", "Health or social care", "Transport",
                                                                 "Retail", "Hospitality", "Manufacturing or construction",
                                                                 "Civil service or local government", "Other"))}
      if(int_names[i]=="sector_coarse") {dat$int_var <- factor(dat$int_var,
                                                               labels=c("Other", "Education", "Health or social care"))}
      if(int_names[i]=="occupation") {dat$int_var <- factor(dat$int_var,
                                                            labels=c("Managers, directors and senior officials",
                                                                     "Professional occupations",
                                                                     "Associate professional and technical occupations",
                                                                     "Administrative and secretarial occupations",
                                                                     "Skilled trades occupations",
                                                                     "Caring, leisure and other service occupations",
                                                                     "Sales and customer service occupations",
                                                                     "Process, plant and machine operatives; and elementary occupations"))}
      if(int_names[i]=="selfemp") {dat$int_var <- factor(dat$int_var,
                                                         labels=c("Employee", "Self-employed"))}
      
      ### reverse order of interaction variable levels
      dat$int_var <- factor(dat$int_var, levels=levels(dat$int_var)[nlevels(dat$int_var):1])
      
      ### derive variable for time since infection
      dat$time <- NA
      dat$time[grep("12-18w", dat$main_var)] <- 1
      dat$time[grep("18-24w", dat$main_var)] <- 2
      dat$time[grep("24-30w", dat$main_var)] <- 3
      dat$time[grep("30-40w", dat$main_var)] <- 4
      dat$time[grep("40-52w", dat$main_var)] <- 5
      dat$time[grep("ge52w", dat$main_var)] <- 6
      
      ### convert time since infection to factor and assign labels
      time_labs <- c("12 to <18 weeks", "18 to <24 weeks",
                     "24 to <30 weeks", "30 to <40 weeks",
                     "40 to <52 weeks", "52+ weeks")
      
      dat$time <- factor(dat$time, labels=time_labs[1:max(dat$time)])
      
      ### reverse order of time since infection levels
      dat$time <- factor(dat$time, levels=levels(dat$time)[nlevels(dat$time):1])
      
      ### set rows for non-defined estimates to NA
      dat$or[dat$or>1000 | dat$lcl==0 | dat$ucl==Inf] <- NA
      dat$lcl[is.na(dat$or)] <- NA
      dat$ucl[is.na(dat$or)] <- NA
      
      ### labels for OR plot
      dat$label <- paste0(
        format(round(dat$or, 2), nsmall=2),
        " (",
        format(round(dat$lcl, 2), nsmall=2),
        " to ",
        format(round(dat$ucl, 2), nsmall=2),
        ")"
      )
      
      ### set labels for non-defined estimates to blank
      dat$label[is.na(dat$or)] <- ""
      
      ### determine number of rows for legend
      legend_rows <- ifelse(nlevels(dat$int_var)<=3, 1,
                            ifelse(nlevels(dat$int_var)<=6, 2, 3))
      
      ### plot ORs
      or_plot <- ggplot(dat, aes(x=or, y=time, group=int_var, colour=int_var, label=label)) +
        geom_vline(xintercept=1, linetype="dashed", colour="grey40", size=0.3) +
        geom_point(position=position_dodge(0.7), size=1.5) +
        geom_errorbarh(aes(xmin=lcl, xmax=ucl, group=int_var), position=position_dodge(0.7), height=0, size=0.4) +
        geom_text(aes(x=max(ucl,na.rm=TRUE)), hjust=0, vjust=0.5, position=position_dodge(0.7), size=2.7) +
        coord_trans(x="log") +
        scale_x_continuous(expand=expansion(mult=c(0.02,0.35))) +
        scale_colour_brewer(palette="Dark2") +
        xlab("Odds ratio (log scale)") +
        theme(
          axis.title.x=element_text(size=9, colour="black", face="bold"),
          axis.text.x = element_text(size=8, colour="black", face="plain"),
          axis.ticks.x=element_line(size=0.5, colour="black"),
          axis.line.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(size=8, colour="black", face="plain"),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          panel.border=element_rect(size=0.5, colour="black", fill=NA),
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(size=9, colour="black", face="bold"),
          legend.title=element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.justification="center",
          legend.text=element_text(size=9, colour="black", face="plain")
        ) +
        guides(colour=guide_legend(reverse=TRUE, nrow=legend_rows, byrow=TRUE))
      
      ggsave(filename=paste0(out_dir, "\\or_", int_names[i], ".jpg"),
             plot=or_plot,
             width=12, height=20, units="cm")
      
    }
    
  }
  
}
