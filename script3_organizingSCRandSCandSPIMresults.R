#run
# script3_processingSPIMresults.R
# script3_processingSCresults.R
# script3_processingSCResults.R
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

load("workspace_script3_processingSCRresults.RData")
load("workspace_script3_processingSCresults.RData")
load("workspace_script3_processingSPIMresults.RData")

#combine SCR, SPIM, and SC results together
colnames(SCRresults)
colnames(SPIMresults)
colnames(SCresults)

colnames(SCRresults)[which(colnames(SCRresults)%in%c("estimate"))]<-"Mean"
colnames(SCRresults)[which(colnames(SCRresults)%in%c("lcl"))]<-"X2.5."
colnames(SCRresults)[which(colnames(SCRresults)%in%c("ucl"))]<-"X97.5."
colnames(SPIMresults)[which(colnames(SPIMresults)%in%c("X50."))]<-"Median"

SCR_SPIM_and_SC_results<-rbind.fill(cbind(SPIMresults,
                                      Model=SPIMresults$PID),
                                cbind(SCresults,Model="SC"),
                                cbind(SCRresults,Model="SCR"))


## change the labels to use "=" instead of ":"
p0.labs <- gsub(":", " =",p0.labs) 
coh.labs <- gsub(":", " =",coh.labs) 

#combine coverage results for SPIM and SC together

#bootsrapped coverages
# colnames(coverages_boot_SC)
# colnames(coverages_boot_SPIM)
# #get the model name
# coverages_boot<-rbind.fill(cbind(coverages_boot_SPIM,
#                                  Model=paste0("SPIM: ",coverages_boot_SPIM$PID)),
#                                 cbind(coverages_boot_SC,Model="SC"))

colnames(coverages_calc_SC)
#coverages_calc_SC$aggregation<-as.numeric(levels(coverages_calc_SC$aggregation))[coverages_calc_SC$aggregation]
colnames(coverages_calc_SPIM)
colnames(coverages_calc_SCR)
coverages_calc<-rbind.fill(cbind(coverages_calc_SPIM,
                                 Model=paste0("SPIM: ",coverages_calc_SPIM$PID)),
                           cbind(coverages_calc_SC,Model="SC"),
                           cbind(coverages_calc_SCR,Model="SCR"))
coverages_calc$Model<-factor(coverages_calc$Model,
                             levels=c("SCR","SC",unique(coverages_calc$Model)[!startsWith(unique(coverages_calc$Model),"SC")]))

##### assign variables of interest for plots #####
param<-c("N","sigma")
PID<-c("antlerssexcollarcoat", "sexcollarcoat", 
       "sexcoat","sexcollar",
       "sex","collar")

temp<-SCR_SPIM_and_SC_results[SCR_SPIM_and_SC_results$param%in%param,]
temp<-temp[temp$Model=="SC"|temp$Model=="SCR"|temp$PID%in%PID,]
temp$Model[!temp$Model%in%c("SC","SCR")]<-paste0("SPIM: ",temp$Model[!temp$Model%in%c("SC","SCR")])
temp$Model<-factor(temp$Model,
                   levels=c("SCR","SC",unique(temp$Model)[!startsWith(unique(temp$Model),"SC")]))


# for SCR models, copy the means to the median column 
temp$Median[temp$Model=="SCR"]<-temp$Mean[temp$Model=="SCR"]

###########summaries##########
#RB 
#dplyr
summaryRB<-temp%>%
  group_by(Model,param,p0,aggregation,cohesion)%>%
  summarize(RB_min=min(RB),
            RB_mean=mean(RB),
            RB_max=max(RB),
            RB_SD=sd(RB))

summaryRB$RB<-paste0(round(summaryRB$RB_min,2),"-",round(summaryRB$RB_max,2)," ; ",
                     round(summaryRB$RB_mean,2)," (",round(summaryRB$RB_SD,2),")")

summary(summaryRB$RB_mean[summaryRB$param=="N"&summaryRB$p0==0.05])
summary(summaryRB$RB_mean[summaryRB$param=="N"&summaryRB$p0==0.2])

summaryRB_N_wide<-summaryRB[summaryRB$param=="N",][,c(1,3,4,5,10)]
summaryRB_N_wide<-pivot_wider(summaryRB_N_wide,names_from=cohesion,values_from = RB)
write.csv(summaryRB[summaryRB$param=="N",],"summaryRB_N_long.csv")
write.csv(summaryRB_N_wide,"summaryRB_N_wide.csv")

summaryRB_Sigma_wide<-summaryRB[summaryRB$param=="sigma",][,c(1,3,4,5,10)]
summaryRB_Sigma_wide<-pivot_wider(summaryRB_Sigma_wide,names_from=cohesion,values_from = RB)
write.csv(summaryRB[summaryRB$param=="sigma",],"summaryRB_Sigma_long.csv")
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
write.csv(summaryCV[summaryCV$param=="N",],"summaryCV_N_long.csv")
write.csv(summaryCV_N_wide,"summaryCV_N_wide.csv")

summaryCV_Sigma_wide<-summaryCV[summaryCV$param=="sigma",][,c(1,3,4,5,10)]
summaryCV_Sigma_wide<-pivot_wider(summaryCV_Sigma_wide,names_from=cohesion,values_from = CV)
write.csv(summaryCV[summaryCV$param=="sigma",],"summaryCV_Sigma_long.csv")
write.csv(summaryCV_Sigma_wide,"summaryCV_Sigma_wide.csv")


#Coverage
coverages_calc_wide<-coverages_calc[,-which(colnames(coverages_calc)%in%c("PID","antler","sex","collar","coat"))]
coverages_calc_wide[coverages_calc_wide$param=="N"&coverages_calc_wide$Model=="SC",]
coverages_calc_wide[coverages_calc_wide$param=="sigma"&coverages_calc_wide$Model=="SC",]

write.csv(coverages_calc[coverages_calc$param=="N",],"summaryCoverage_N_long.csv")
write.csv(coverages_calc[coverages_calc$param=="sigma",],"summaryCoverage_sigma_long.csv")
write.csv(coverages_calc_wide,"coverages_calc_wide.csv")

coverages_calc_N_wide_wide<-pivot_wider(coverages_calc_wide[coverages_calc_wide$param=="N",],
                                      names_from=cohesion,
                                      values_from=coverage)
coverages_calc_sigma_wide_wide<-pivot_wider(coverages_calc_wide[coverages_calc_wide$param=="sigma",],
                                        names_from=cohesion,
                                        values_from=coverage)

write.csv(coverages_calc_N_wide_wide,"coverages_N_calc_wide.csv")
write.csv(coverages_calc_sigma_wide_wide,"coverages_sigma_calc_wide.csv")


##### plot####
colBluefunc <- colorRampPalette(c("lightblue","blue"))
colBluefunc(length(levels(temp$Model))-2)
cbbPalette <- c("white","pink",colBluefunc(length(levels(temp$Model))-2))# "#E69F00", "#56B4E9", "#009E73", "#F0E442")



#comparing SCR, SC, and SPIM for the actual estimates
plot_Nmean<-ggplot(temp[temp$param=="N",],
                  aes(x=aggregation, y=Median,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Mean/Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free_x", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigmaMean<-ggplot(temp[temp$param=="sigma",],
                      aes(x=aggregation, y=Median,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=3,lty=2)+
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Mean/Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", scales = "free_x", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_means_p05<-ggplot(temp[temp$p0=="0.05",],
       aes(x=aggregation, y=Median,fill=Model)) + 
  geom_hline(aes(yintercept = Truth),lty=2)+
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  labs(title="",x="Aggregation (Group Size)",y="Mean/Median Estimate")+
  facet_grid(rows=vars(param),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(param = c("N",expression(sigma)) ,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_means_p20<-ggplot(temp[temp$p0=="0.2",],
                       aes(x=aggregation, y=Median,fill=Model)) + 
  geom_hline(aes(yintercept = Truth),lty=2)+
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  labs(title="",x="Aggregation (Group Size)",y="Mean/Median Estimate")+
  facet_grid(rows=vars(param),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(param=c("N",expression(sigma)),cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
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
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  # ylim(-2,30)+
  labs(title="Abundance (N)",
    x="Aggregation (Group Size)",y="Relative Bias")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",  scales = "free_x", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigma_rb<-ggplot(temp[temp$param=="sigma",],
                           aes(x=aggregation, y=RB,fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
 # ylim(-2,3.5)+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Relative Bias")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

#bias ax params
ggplot(temp[temp$Model=="SCR",], # or SC
       aes(x=aggregation, y=RB,fill=param)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  # ylim(-2,3.5)+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Relative Bias")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

#cv

plot_N_cv<-ggplot(temp[temp$param=="N",],
                       aes(x=aggregation, y=CoV,
                           fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
 # geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_N_cv_Maxcutoff<-ggplot(temp[temp$param=="N",],
                  aes(x=aggregation, y=CoV,
                      fill=Model)) + 
  ylim(c(0,1.5))+
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  # geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sig_cv<-ggplot(temp[temp$param=="sigma",],
                         aes(x=aggregation, y=CoV,
                             fill=Model)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
 # ylim(0,1)+
   # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sig_cv_Maxcutoff<-ggplot(temp[temp$param=="sigma",],
                    aes(x=aggregation, y=CoV,
                        fill=Model)) + 
  ylim(c(0,0.5))+
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  # ylim(0,1)+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#precision ax params
ggplot(temp[temp$Model=="SC",], # or SC
       aes(x=aggregation, y=CoV,fill=param)) + 
  scale_fill_manual(values=cbbPalette)+
  geom_boxplot(position=position_dodge(1))+
  # ylim(-2,3.5)+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coefficient of Variation")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,



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
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), scales = "free",
             labeller=labeller(p0 = p0.labs,cohesion=coh.labs))+
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
  #ylim(c(0.25,1))+
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=coverages_calc[coverages_calc$param=="sigma",], 
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(x=aggregation, y=coverage,
                 shape=Model,group=Model,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$param=="sigma",],
            position=position_dodge(width=2),
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), scales = "free",
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))



plot_coverage_p05<-ggplot() + 
  #ylim(c(0.25,1))+
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=coverages_calc[coverages_calc$p0=="0.05",], 
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(x=aggregation, y=coverage,
                 shape=Model,group=Model,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$p0=="0.05",],
            position=position_dodge(width=2),
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(param),cols=vars(cohesion),switch="y", scales = "free",
             labeller=labeller(cohesion=coh.labs))+
  labs(#title="Sigma (\u03c3)",
       x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_coverage_p20<-ggplot() + 
  #ylim(c(0.25,1))+
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=coverages_calc[coverages_calc$p0=="0.2",], 
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(x=aggregation, y=coverage,
                 shape=Model,group=Model,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$p0=="0.2",],
            position=position_dodge(width=2),
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(param),cols=vars(cohesion),switch="y", scales = "free",
             labeller=labeller(cohesion=coh.labs))+
  labs(#title="Sigma (\u03c3)",
    x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))



######plots #####

plot_Nmean
plot_sigmaMean
plot_means_p05
plot_means_p20

plot_N_rb
plot_sigma_rb

plot_N_cv
plot_N_cv_Maxcutoff
plot_sig_cv
plot_sig_cv_Maxcutoff

plot_N_coverage
plot_sigma_coverage
plot_coverage_p05
plot_coverage_p20
