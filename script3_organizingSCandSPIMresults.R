#run
# script3_processingSPIMresults.R
# script3_processingSCresults.R


#combine SPIM and SC results together
colnames(SPIMresults)
colnames(SCresults)

colnames(SPIMresults)[which(colnames(SPIMresults)%in%c("X50."))]<-"Median"
SPIM_and_SC_results<-rbind.fill(cbind(SPIMresults,Model=paste0("SPIM: ",rowSums(SPIMresults[,25:28]))),
                                cbind(SCresults,Model="SC"))

#combine coverage results for SPIM and SC together

colnames(coverages_boot_SC)
colnames(coverages_boot_SPIM)
coverages_boot<-rbind.fill(cbind(coverages_boot_SPIM,Model=paste0("SPIM: ",rowSums(coverages_boot_SPIM[,7:10]))),
                                cbind(coverages_boot_SC,Model="SC"))

colnames(coverages_calc_SC)
colnames(coverages_calc_SPIM)
coverages_calc<-rbind.fill(cbind(coverages_calc_SPIM,Model=paste0("SPIM: ",rowSums(coverages_calc_SPIM[,7:10]))),
                           cbind(coverages_calc_SC,Model="SC"))

##### assign variables of interest for plots #####
param<-c("N","sigma")
PID<-c("antlerssexcollarcoat", "sexcollarcoat")

temp<-SPIM_and_SC_results[SPIM_and_SC_results$param%in%param,]
temp<-temp[temp$Model=="SC"|temp$PID%in%PID,]


#comparing SC ad SPIM for the actual estimates
plot_Nmed<-ggplot(temp[temp$param=="N",],
                  aes(x=aggregation, y=Median,color=Model)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigmaMed<-ggplot(temp[temp$param=="sigma",],
                  aes(x=aggregation, y=Median,color=Model)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=3,lty=2)+
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,



#rb
plot_N_rb<-ggplot(temp[temp$param=="N",],
                       aes(x=aggregation, y=RB,color=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="RB")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigma_rb<-ggplot(temp[temp$param=="sigma",],
                           aes(x=aggregation, y=RB,color=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

#cv

plot_N_cv<-ggplot(temp[temp$param=="N",],
                       aes(x=aggregation, y=CoV,
                           color=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="CV")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sig_cv<-ggplot(temp[temp$param=="sigma",],
                         aes(x=aggregation, y=CoV,
                             color=Model)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="CV")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#rmse


### coverage ####
# 
#temporarily manipulate the coverage data so that they dont overlap when plotted
coverages_boot$aggregation<-as.numeric(levels(coverages_boot$aggregation))[coverages_boot$aggregation]
coverages_boot$aggregation[coverages_boot$Model=="SC"]<-coverages_boot$aggregation[coverages_boot$Model=="SC"]-0.8
coverages_boot$aggregation[coverages_boot$Model=="SPIM: 3"]<-coverages_boot$aggregation[coverages_boot$Model=="SPIM: 3"]+0.8

coverages_calc$aggregation<-as.numeric(levels(coverages_calc$aggregation))[coverages_calc$aggregation]
coverages_calc$aggregation[coverages_calc$Model=="SC"]<-coverages_calc$aggregation[coverages_calc$Model=="SC"]-0.8
coverages_calc$aggregation[coverages_calc$Model=="SPIM: 3"]<-coverages_calc$aggregation[coverages_calc$Model=="SPIM: 3"]+0.8


plot_N_coverage<-ggplot(coverages_boot[coverages_boot$param=="N",], 
                             aes(x=aggregation, y=coverage,
                                 group=Model)) +#
  stat_summary(fun.min = function(x) min(x), 
               fun.max = function(x) max(x), 
               geom = "linerange",size=1,
               aes(color=Model),show.legend = FALSE) +
  geom_point(data=coverages_calc[coverages_calc$param=="N",], 
             size=3,
             aes(x=aggregation, y=coverage,
                 shape=cohesion,group=cohesion,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$param=="N",],
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
 
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
                                                                 cohesion=coh.labs))+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))



plot_sigma_coverage<-ggplot(coverages_boot[coverages_boot$param=="sigma",], 
                        aes(x=aggregation, y=coverage,
                            group=Model)) +#
  stat_summary(fun.min = function(x) min(x), 
               fun.max = function(x) max(x), 
               geom = "linerange",size=1,
               aes(color=Model),show.legend = FALSE) +
  geom_point(data=coverages_calc[coverages_calc$param=="sigma",], 
             size=3,
             aes(x=aggregation, y=coverage,
                 shape=cohesion,group=cohesion,color=Model))+
  geom_line(data=coverages_calc[coverages_calc$param=="sigma",],
            aes(x=aggregation, y=coverage, group=Model,color=Model))+ 
  
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
                                                                 cohesion=coh.labs))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#return coverage data so that they dont overlap when plotted
coverages_boot$aggregation[coverages_boot$Model=="SC"]<-coverages_boot$aggregation[coverages_boot$Model=="SC"]+0.8
coverages_boot$aggregation[coverages_boot$Model=="SPIM: 3"]<-coverages_boot$aggregation[coverages_boot$Model=="SPIM: 3"]-0.8
coverages_boot$aggregation<-as.factor(coverages_boot$aggregation)

coverages_calc$aggregation[coverages_calc$Model=="SC"]<-coverages_calc$aggregation[coverages_calc$Model=="SC"]+0.8
coverages_calc$aggregation[coverages_calc$Model=="SPIM: 3"]<-coverages_calc$aggregation[coverages_calc$Model=="SPIM: 3"]-0.8
coverages_calc$aggregation<-as.factor(coverages_calc$aggregation)

