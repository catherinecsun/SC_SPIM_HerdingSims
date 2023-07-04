#run
# script3_processingSPIMresults.R
# script3_processingSCresults.R
library(tidyr)

#combine SPIM and SC results together
colnames(SPIMresults)
colnames(SCresults)

colnames(SPIMresults)[which(colnames(SPIMresults)%in%c("X50."))]<-"Median"
SPIM_and_SC_results<-rbind.fill(cbind(SPIMresults,
                                      Model=SPIMresults$PID),
                                cbind(SCresults,Model="SC"))

#combine coverage results for SPIM and SC together
colnames(coverages_calc_SC)
coverages_calc_SC$aggregation<-as.numeric(levels(coverages_calc_SC$aggregation))[coverages_calc_SC$aggregation]
colnames(coverages_calc_SPIM)
coverages_calc<-rbind.fill(cbind(coverages_calc_SPIM,
                                 Model=paste0("SPIM: ",coverages_calc_SPIM$PID)),
                           cbind(coverages_calc_SC,Model="SC"))

##### assign variables of interest for plots #####
param<-c("N","sigma")
PID<-c("antlerssexcollarcoat", "sexcollarcoat", 
       "sexcoat","sexcollar",
       "sex","collar")

temp<-SPIM_and_SC_results[SPIM_and_SC_results$param%in%param,]
temp<-temp[temp$Model=="SC"|temp$PID%in%PID,]
temp$Model[!temp$Model=="SC"]<-paste0("SPIM: ",temp$Model[!temp$Model=="SC"])
temp$Model<-factor(temp$Model,
                   levels=c("SC",unique(temp$Model)[!startsWith(unique(temp$Model),"SC")]))
###########summaries##########
#RB
summaryRB<-temp%>%
  group_by(Model,param,p0,aggregation,cohesion)%>%
  summarize(RB_min=min(RB),
            RB_mean=mean(RB),
            RB_max=max(RB),
            RB_SD=sd(RB))

summaryRB$RB<-paste0(round(summaryRB$RB_min,2),"-",round(summaryRB$RB_max,2)," ; ",
                     round(summaryRB$RB_mean,2)," (",round(summaryRB$RB_SD,2),")")

summaryRB_N_wide<-summaryRB[summaryRB$param=="N",][,c(1,3,4,5,10)]
summaryRB_N_wide<-pivot_wider(summaryRB_N_wide,names_from=cohesion,values_from = RB)
write.csv(summaryRB_N_wide,"summaryRB_N_wide.csv")

summaryRB_Sigma_wide<-summaryRB[summaryRB$param=="sigma",][,c(1,3,4,5,10)]
summaryRB_Sigma_wide<-pivot_wider(summaryRB_Sigma_wide,names_from=cohesion,values_from = RB)
write.csv(summaryRB_Sigma_wide,"summaryRB_Sigma_wide.csv")

#CV
summaryCV<-temp%>%
  group_by(Model,param,p0,aggregation,cohesion)%>%
  summarize(CV_min=min(CoV),
            CV_mean=mean(CoV),
            CV_max=max(CoV),
            CV_SD=sd(CoV))

summaryCV$CV<-paste0(round(summaryCV$CV_min,2),"-",round(summaryCV$CV_max,2)," ; ",
                     round(summaryCV$CV_mean,2)," (",round(summaryCV$CV_SD,2),")")

summaryCV_N_wide<-summaryCV[summaryCV$param=="N",][,c(1,3,4,5,10)]
summaryCV_N_wide<-pivot_wider(summaryCV_N_wide,names_from=cohesion,values_from = CV)
write.csv(summaryCV_N_wide,"summaryCV_N_wide.csv")

summaryCV_Sigma_wide<-summaryCV[summaryCV$param=="sigma",][,c(1,3,4,5,10)]
summaryCV_Sigma_wide<-pivot_wider(summaryCV_Sigma_wide,names_from=cohesion,values_from = CV)
write.csv(summaryCV_Sigma_wide,"summaryCV_Sigma_wide.csv")


#Coverage
coverages_calc_wide<-coverages_calc[,-which(colnames(coverages_calc_wide)%in%c("PID","antler","sex","collar","coat"))]
coverages_calc_wide[coverages_calc_wide$param=="N"&coverages_calc_wide$Model=="SC",]
coverages_calc_wide[coverages_calc_wide$param=="sigma"&coverages_calc_wide$Model=="SC",]

write.csv(coverages_calc_wide,"coverages_calc_wide.csv")

coverages_calc_wide_wide<-pivot_wider(coverages_calc_wide,
                                      names_from=param,
                                      values_from=coverage)


##### plot####
colBluefunc <- colorRampPalette(c("lightblue","blue"))
colBluefunc(length(levels(temp$Model))-1)
cbbPalette <- c("white",colBluefunc(length(levels(temp$Model))-1))# "#E69F00", "#56B4E9", "#009E73", "#F0E442")



#comparing SC ad SPIM for the actual estimates
plot_Nmed<-ggplot(temp[temp$param=="N",],
                  aes(x=aggregation, y=Median,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigmaMed<-ggplot(temp[temp$param=="sigma",],
                      aes(x=aggregation, y=Median,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_hline(yintercept=3,lty=2)+
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#rb
temp[temp$param=="N",]%>%
  group_by(p0,cohesion,aggregation)%>%
  summarize(RB_min=min(RB),
            RB_max=max(RB),
            RB_mean=mean(RB),
            RB_SD=sd(RB))
plot_N_rb<-ggplot(temp[temp$param=="N",],
                       aes(x=aggregation, y=RB,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_hline(yintercept=0, linetype="dashed")+
 # ylim(-2,30)+
  labs(#title="Abundance (N)",
    x="Aggregation (Group Size)",y="Relative Bias")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigma_rb<-ggplot(temp[temp$param=="sigma",],
                           aes(x=aggregation, y=RB,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
 # ylim(-2,3.5)+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Relative Bias")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#cv

plot_N_cv<-ggplot(temp[temp$param=="N",],
                       aes(x=aggregation, y=CoV,
                           fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
 # geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sig_cv<-ggplot(temp[temp$param=="sigma",],
                         aes(x=aggregation, y=CoV,
                             fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
 # ylim(0,1)+
   # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#rmse


### coverage ####
# 
cbbPalette2 <-cbbPalette
cbbPalette2[1]<-"black"

coverages_calc$Model<-factor(coverages_calc$Model,levels=levels(temp$Model))


plot_N_coverage<-ggplot() +
  #coverages_boot[coverages_boot$param=="N",],aes(x=aggregation, y=coverage,group=Model)
  # stat_summary(fun.min = function(x) min(x),
  #              fun.max = function(x) max(x),
  #              geom = "linerange",size=1,
  #              aes(color=Model),show.legend = FALSE) +
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=coverages_calc[coverages_calc$param=="N",],
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(x=aggregation, y=coverage,
                 shape=Model,group=Model,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$param=="N",],
            position=position_dodge(width=2),
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
                                                                 cohesion=coh.labs))+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sigma_coverage<-ggplot() +
  #coverages_boot[coverages_boot$param=="sigma",], aes(x=aggregation, y=coverage,group=Model)
  # stat_summary(fun.min = function(x) min(x), 
  #              fun.max = function(x) max(x), 
  #              geom = "linerange",size=1,
  #              aes(color=Model),show.legend = FALSE) +
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=coverages_calc[coverages_calc$param=="sigma",], 
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(x=aggregation, y=coverage,
                 shape=Model,group=Model,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$param=="sigma",],
            position=position_dodge(width=2),
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
                                                                 cohesion=coh.labs))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


######plots #####

plot_Nmed
plot_sigmaMed
plot_N_rb
plot_sigma_rb
plot_N_cv
plot_sig_cv
plot_N_coverage
plot_sigma_coverage
