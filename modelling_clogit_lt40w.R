modelling.clogit.lt40w <- function(dset, chart_cols, out_dir, outcome,
                                   covs, run_mod1, run_mod2) {
  
  #################### MODEL 1: IGNORING TIME SINCE INFECTION ####################
  
  if(run_mod1==TRUE) {
    
    ### fit unadjusted model
    mod1_unadj <- fit.clogit.model(outcome=outcome,
                                   exposures="grp2",
                                   covariates=NULL,
                                   id="participant_id",
                                   dataset=dset)
    
    ### fit adjusted model
    mod1_adj <- fit.clogit.model(outcome=outcome,
                                 exposures="grp2",
                                 covariates=covs,
                                 id="participant_id",
                                 dataset=dset,
                                 run.wald=FALSE)
    
    ### fit adjusted model - baseline
    mod1_adj_baseline <- fit.clogit.model(outcome=outcome,
                                          exposures="1",
                                          covariates=covs,
                                          id="participant_id",
                                          dataset=dset,
                                          run.wald=FALSE)
    
    ### write out model outputs
    write.csv(mod1_unadj$coeff, file=paste0(out_dir, "\\coeffs_unadj_mod1.csv"))
    write.csv(mod1_unadj$wald, file=paste0(out_dir, "\\wald_unadj_mod1.csv"))
    
    write.csv(mod1_adj$coeff, file=paste0(out_dir, "\\coeffs_adj_mod1.csv"))
    write.csv(lrtest(mod1_adj$mod, mod1_adj_baseline$mod)[5][2,],
              file=paste0(out_dir, "\\wald_adj_mod1.csv"))
    
    ### model metrics
    mod1_metrics <- data.frame(model = c("unadjusted", "adjusted"),
                               aic = c(mod1_unadj$aic, mod1_adj$aic),
                               bic = c(mod1_unadj$bic, mod1_adj$bic))
    
    write.csv(mod1_metrics, file=paste0(out_dir, "\\metrics_mod1.csv"), row.names=FALSE)
    
    ### coeffs for unadjusted model
    or1_unadj <- mod1_unadj$coeff[grep("grp2", rownames(mod1_unadj$coeff)),]
    or1_unadj$model <- "1_unadjusted"
    or1_unadj$grp2 <- rownames(or1_unadj)
    
    ### coeffs for adjusted model
    or1_adj <- mod1_adj$coeff[grep("grp2", rownames(mod1_adj$coeff)),]
    or1_adj$model <- "2_adjusted"
    or1_adj$grp2 <- rownames(or1_adj)
    
    ### stack coeffs for unadjusted and adjusted models
    or1 <- rbind(or1_unadj, or1_adj)
    
    ### calculate ORs
    or1$or <- exp(or1[,1])
    or1$lcl <- exp(or1[,1] - 1.96*or1[,3])
    or1$ucl <- exp(or1[,1] + 1.96*or1[,3])
    
    ### re-label variables
    or1$model <- factor(or1$model,
                        labels=c("Unadjusted", "Adjusted"))
    
    or1$grp2 <- factor(or1$grp2,
                       labels=c("Infected <12 weeks ago", "Infected 12+ weeks ago",
                                "Current Long Covid", "Previous Long Covid"))
    
    or1$grp2 <- factor(or1$grp2, levels=levels(or1$grp2)[nlevels(or1$grp2):1])
    
    ### write out ORs
    write.csv(or1[,6:10], file=paste0(out_dir, "\\or_mod1.csv"), row.names=FALSE)
    
    ### labels for OR plot
    or1$label <- paste0(
      format(round(or1$or, 2), nsmall=2),
      " (",
      format(round(or1$lcl, 2), nsmall=2),
      " to ",
      format(round(or1$ucl, 2), nsmall=2),
      ")"
    )
    
    ### set rows for non-defined estimates to NA
    or1$or[or1$lcl==0 | or1$ucl==Inf] <- NA
    or1$lcl[is.na(or1$or)] <- NA
    or1$ucl[is.na(or1$or)] <- NA
    or1$label[is.na(or1$or)] <- ""
    
    ### plot ORs - faceted on unadjusted vs. adjusted
    or1_plot_facet <- ggplot(or1, aes(x=or, y=grp2, label=label)) +
      geom_vline(xintercept=1, linetype="dashed", colour="grey40", size=0.3) +
      geom_point(size=2.0, colour=rep(chart_cols,2)) +
      geom_errorbarh(aes(xmin=lcl, xmax=ucl), height=0, size=0.5, colour=rep(chart_cols,2)) +
      geom_text(aes(x=or), hjust=0.5, vjust=-1, size=2.7, colour=rep(chart_cols,2)) +
      facet_wrap(~model, nrow=1) +
      coord_trans(x="log") +
      scale_x_continuous(expand=expansion(mult=c(0.2,0.2))) +
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
      )
    
    ggsave(filename=paste0(out_dir, "\\or_mod1_faceted.jpg"),
           plot=or1_plot_facet,
           width=18, height=11, units="cm")
    
    ### plot ORs - adjusted estimates only
    or1_plot_adj <- ggplot(or1[or1$model=="Adjusted",], aes(x=or, y=grp2, label=label)) +
      geom_vline(xintercept=1, linetype="dashed", colour="grey40", size=0.3) +
      geom_point(size=2.0, colour=chart_cols) +
      geom_errorbarh(aes(xmin=lcl, xmax=ucl), height=0, size=0.5, colour=chart_cols) +
      geom_text(aes(x=or), hjust=0.5, vjust=-1, size=2.7, colour=chart_cols) +
      coord_trans(x="log") +
      scale_x_continuous(expand=expansion(mult=c(0.1,0.1))) +
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
      )
    
    ggsave(filename=paste0(out_dir, "\\or_mod1_adjusted.jpg"),
           plot=or1_plot_adj,
           width=12, height=10, units="cm")
    
  }
  
  ################# MODEL 2: ACCOUNTING FOR TIME SINCE INFECTION #################
  
  if(run_mod2==TRUE) {
    
    ### fit unadjusted model
    mod2_unadj <- fit.clogit.model(outcome=outcome,
                                   exposures="exposure",
                                   covariates=NULL,
                                   id="participant_id",
                                   dataset=dset)
    
    ### fit adjusted model
    mod2_adj <- fit.clogit.model(outcome=outcome,
                                 exposures="exposure",
                                 covariates=covs,
                                 id="participant_id",
                                 dataset=dset,
                                 run.wald=FALSE)
    
    ### fit adjusted model - baseline
    mod2_adj_baseline <- fit.clogit.model(outcome=outcome,
                                          exposures="1",
                                          covariates=covs,
                                          id="participant_id",
                                          dataset=dset,
                                          run.wald=FALSE)
    
    ### write out model outputs
    write.csv(mod2_unadj$coeff, file=paste0(out_dir, "\\coeffs_unadj_mod2.csv"))
    write.csv(mod2_unadj$wald, file=paste0(out_dir, "\\wald_unadj_mod2.csv"))
    
    write.csv(mod2_adj$coeff, file=paste0(out_dir, "\\coeffs_adj_mod2.csv"))
    write.csv(lrtest(mod2_adj$mod, mod2_adj_baseline$mod)[5][2,],
              file=paste0(out_dir, "\\wald_adj_mod2.csv"))
    
    ### model metrics
    mod2_metrics <- data.frame(model = c("unadjusted", "adjusted"),
                               aic = c(mod2_unadj$aic, mod2_adj$aic),
                               bic = c(mod2_unadj$bic, mod2_adj$bic))
    
    write.csv(mod2_metrics, file=paste0(out_dir, "\\metrics_mod2.csv"), row.names=FALSE)
    
    ### coeffs for unadjusted model
    or2_unadj <- mod2_unadj$coeff[grep("exposure", rownames(mod2_unadj$coeff)),]
    or2_unadj$model <- "1_unadjusted"
    or2_unadj$exposure <- rownames(or2_unadj)
    
    ### coeffs for adjusted model
    or2_adj <- mod2_adj$coeff[grep("exposure", rownames(mod2_adj$coeff)),]
    or2_adj$model <- "2_adjusted"
    or2_adj$exposure <- rownames(or2_adj)
    
    ### stack coeffs for unadjusted and adjusted models
    or2 <- rbind(or2_unadj, or2_adj)
    
    ### calculate ORs
    or2$or <- exp(or2[,1])
    or2$lcl <- exp(or2[,1] - 1.96*or2[,3])
    or2$ucl <- exp(or2[,1] + 1.96*or2[,3])
    
    ### derive variable for time since infection
    or2$time <- NA
    or2$time[grep("lt12w", or2$exposure)] <- 1
    or2$time[grep("12-18w", or2$exposure)] <- 2
    or2$time[grep("18-24w", or2$exposure)] <- 3
    or2$time[grep("24-30w", or2$exposure)] <- 4
    or2$time[grep("30-40w", or2$exposure)] <- 5
    
    ### derive variable for infection/LC status
    or2$status <- NA
    or2$status[grep("lt12w", or2$exposure)] <- 1
    or2$status[grep("nolc", or2$exposure)] <- 2
    or2$status[grep("currlc", or2$exposure)] <- 3
    or2$status[grep("prevlc", or2$exposure)] <- 4
    
    ### re-label variables
    or2$model <- factor(or2$model,
                        labels=c("Unadjusted", "Adjusted"))
    
    or2$time <- factor(or2$time,
                       labels=c("<12 weeks", "12 to <18 weeks",
                                "18 to <24 weeks", "24 to <30 weeks",
                                "30 to <40 weeks"))
    
    or2$time <- factor(or2$time, levels=levels(or2$time)[nlevels(or2$time):1])
    
    or2$status <- factor(or2$status,
                         labels=c("Infected <12 weeks ago", "Infected 12+ weeks ago",
                                  "Current Long Covid", "Previous Long Covid"))
    
    or2$status <- factor(or2$status, levels=levels(or2$status)[nlevels(or2$status):1])
    
    or2$exposure <- factor(or2$exposure,
                           labels=c("<12 weeks",
                                    "12 to <18 weeks, previously infected", "18 to <24 weeks, previously infected",
                                    "24 to <30 weeks, previously infected", "30 to <40 weeks, previously infected",
                                    "12 to <18 weeks, current Long Covid", "18 to <24 weeks, current Long Covid",
                                    "24 to <30 weeks, current Long Covid", "30 to <40 weeks, current Long Covid",
                                    "12 to <18 weeks, previous Long Covid", "18 to <24 weeks, previous Long Covid",
                                    "24 to <30 weeks, previous Long Covid", "30 to <40 weeks, previous Long Covid"))
    
    or2$exposure <- factor(or2$exposure, levels=levels(or2$exposure)[nlevels(or2$exposure):1])
    
    ### write out ORs
    write.csv(or2[,6:10], file=paste0(out_dir, "\\or_mod2.csv"), row.names=FALSE)
    
    ### labels for OR plot
    or2$label <- paste0(
      format(round(or2$or, 2), nsmall=2),
      " (",
      format(round(or2$lcl, 2), nsmall=2),
      " to ",
      format(round(or2$ucl, 2), nsmall=2),
      ")"
    )
    
    ### set rows for non-defined estimates to NA
    or2$or[or2$lcl==0 | or2$ucl==Inf] <- NA
    or2$lcl[is.na(or2$or)] <- NA
    or2$ucl[is.na(or2$or)] <- NA
    or2$label[is.na(or2$or)] <- ""
    
    ### plot ORs - faceted on unadjusted vs. adjusted
    or2_plot_facet <- ggplot(or2, aes(x=or, y=time, group=status, colour=status, label=label)) +
      geom_vline(xintercept=1, linetype="dashed", colour="grey40", size=0.3) +
      geom_point(position=position_dodge(0.7), size=1.5) +
      geom_errorbarh(aes(xmin=lcl, xmax=ucl, group=status), position=position_dodge(0.7), height=0, size=0.4) +
      geom_text(aes(x=max(ucl,na.rm=TRUE)), hjust=0, vjust=0.5, position=position_dodge(0.7), size=2.7) +
      facet_wrap(~model, nrow=1) +
      coord_trans(x="log") +
      scale_x_continuous(expand=expansion(mult=c(0.02,0.5))) +
      scale_colour_manual(values=chart_cols[4:1]) +
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
      guides(colour=guide_legend(reverse=TRUE, nrow=2, byrow=TRUE))
    
    ggsave(filename=paste0(out_dir, "\\or_mod2_faceted.jpg"),
           plot=or2_plot_facet,
           width=18, height=20, units="cm")
    
    ### plot ORs - adjusted estimates only
    or2_plot_adj <- ggplot(or2[or2$model=="Adjusted",], aes(x=or, y=time, group=status, colour=status, label=label)) +
      geom_vline(xintercept=1, linetype="dashed", colour="grey40", size=0.3) +
      geom_point(position=position_dodge(0.7), size=1.5) +
      geom_errorbarh(aes(xmin=lcl, xmax=ucl, group=status), position=position_dodge(0.7), height=0, size=0.4) +
      geom_text(aes(x=max(ucl,na.rm=TRUE)), hjust=0, vjust=0.5, position=position_dodge(0.7), size=2.7) +
      coord_trans(x="log") +
      scale_x_continuous(expand=expansion(mult=c(0.02,0.35))) +
      scale_colour_manual(values=chart_cols[4:1]) +
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
      guides(colour=guide_legend(reverse=TRUE, nrow=2, byrow=TRUE))
    
    ggsave(filename=paste0(out_dir, "\\or_mod2_adjusted.jpg"),
           plot=or2_plot_adj,
           width=12, height=20, units="cm")
    
  }
  
}
