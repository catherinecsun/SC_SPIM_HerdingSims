#packages
library(plyr)
library(ggplot2)
library(ggpubr)
library(Metrics)

####custom bootstrap function####
booties<-function(object,Bsamples=1000){
  boot.samples<-matrix(sample(object, size = Bsamples * length(object), replace = TRUE), Bsamples, length(object))
  return(apply(boot.samples, 1, mean))
}

### read in simulation results####
filePattern<-"SC_simResults_scenario"

results_SC<-list.files("SCresults") #the folder that spim results are in
results_SC<-results_SC[grepl(filePattern,results_SC)]


SCresults<-data.frame()

for(i in 1:length(results_SC)){
  tmp<-read.csv(paste0("SCresults/",results_SC[i]))
  SCresults <-rbind.fill(SCresults,tmp)
}

#remove garbage first column
SCresults<-SCresults[,-1]

#double check things
table(SCresults$scenario)

#remove any sims that didnt converge
which(SCresults$Rhat_pt>1.1)
convergeFail<-unique(SCresults[which(SCresults$Rhat_pt>1.1),colnames(SCresults)%in%c("scenario","sim")])
convergeFail
convergeFail_remove<-c()
if(nrow(convergeFail)>0){
  for(f in 1:nrow(convergeFail)){
    convergeFail_remove<-c(convergeFail_remove,
                           intersect(which(SCresults$scenario==convergeFail$scenario[f]),
                                     which(SCresults$sim==convergeFail$sim[f])))
    
    #btw fwiw i think match_df function in plyr would have worked for this too
  }
  SCresults<-SCresults[-convergeFail_remove,]
}

#double check things
table(SCresults$scenario)

# add columns for the parameter/scenario conditions
#use a different version of parm_combos that focuses 
#only on the 18 combos that are different (and that we evaluated)
parm_combos2<-parm_combos[-c(3:8),]

SCresults$p0<-NA
SCresults$cohesion<-NA
SCresults$aggregation<-NA
for(r in 1:nrow(parm_combos2)){
  SCresults$p0[SCresults$scenario==r]<-parm_combos2$p0[r]
  SCresults$cohesion[SCresults$scenario==r]<-parm_combos2$cohesion[r]
  SCresults$aggregation[SCresults$scenario==r]<-parm_combos2$aggregation[r]
}

SCresults$cohesion<-as.factor(SCresults$cohesion)
SCresults$p0<-as.factor(SCresults$p0)
SCresults$aggregation<-as.factor(SCresults$aggregation)



### relative bias  ####
#(truth-estimate)/estimate
area<-(range(sim_data[[1]][[1]]$mask[,1])[2]-range(sim_data[[1]][[1]]$mask[,1])[1])*
  (range(sim_data[[1]][[1]]$mask[,2])[2]-range(sim_data[[1]][[1]]$mask[,2])[1])
truth<-c("D"=N.inds/area,"N"=N.inds,"psi"=N.inds/M,"sigma"=3) #"lam0"=c(0.05,0.20),

SCresults$Truth[which(SCresults$param%in%names(truth))]<-truth[match(SCresults$param[which(SCresults$param%in%names(truth))],names(truth))]
SCresults$Truth[which(SCresults$param=="lam0"&SCresults$scenario %% 2 == 0)]<-0.2
SCresults$Truth[which(SCresults$param=="lam0"&SCresults$scenario %% 2 == 1)]<-0.05
SCresults$RB<-(SCresults$Truth-SCresults$Mean)/SCresults$Mean
  

### Coefficient of Variation, precision ####
#SD/estimate
SCresults$CoV<-SCresults$St.Dev/SCresults$Mean


##rmse
SCresults$rmse<-NA
for(n in 1:length(unique(SCresults$scenario))){
  SCresults$rmse[SCresults$scenario==n&SCresults$param=="D"]<-
    rmse(truth[names(truth)=="D"], SCresults$Mean[SCresults$scenario==n&SCresults$param=="D"])
  SCresults$rmse[SCresults$scenario==n&SCresults$param=="N"]<-
    rmse(truth[names(truth)=="N"], SCresults$Mean[SCresults$scenario==n&SCresults$param=="N"])
  SCresults$rmse[SCresults$scenario==n&SCresults$param=="sigma"]<-
    rmse(truth[names(truth)=="sigma"], SCresults$Mean[SCresults$scenario==n&SCresults$param=="sigma"])
}

#coverage
SCresults$coverage<-NA
for(co in 1:nrow(SCresults)){
  if(SCresults$param[co]=="D"){
    SCresults$coverage[co]<-ifelse(truth[names(truth)=="D"]>SCresults$X95.CI_low[co]&SCresults$X95.CI_upp[co]>truth[names(truth)=="D"],
                             1,0)
  }
  if(SCresults$param[co]=="N"){
    SCresults$coverage[co]<-ifelse(truth[names(truth)=="N"]>SCresults$X95.CI_low[co]&SCresults$X95.CI_upp[co]>truth[names(truth)=="N"],
                                         1,0)
  }
  if(SCresults$param[co]=="sigma"){
    SCresults$coverage[co]<-ifelse(truth[names(truth)=="sigma"]>SCresults$X95.CI_low[co]&SCresults$X95.CI_upp[co]>truth[names(truth)=="sigma"],
                             1,0)
  }
}

library(dplyr)

# caluclating coverage
coverages_calc_SC<-SCresults[SCresults$param=="N"|SCresults$param=="sigma",]%>%
  group_by(cohesion, aggregation, p0,param)%>%
  summarise(coverage=sum(coverage)/100)

#bootstrapping coverage
coverages_boot_SC<-SCresults[SCresults$param=="N"|SCresults$param=="sigma",]%>%
  group_by(cohesion, aggregation, p0,param)%>%
  summarise(coverage=booties(coverage))


### relative variance ####

###table####
# SC_table<-coverages_SC
# SC_table<-left_join(SC_table,
#                unique(SCresults[SCresults$param%in%unique(SC_table$param),c(match(colnames(coverages_SC)[-ncol(coverages_SC)],colnames(SCresults)),which(colnames(SCresults)=="rmse"))]),
#                by=colnames(SC_table)[-ncol(SC_table)])
# SC_table$rmse<-round(SC_table$rmse,2)

### plot N ####
#   HOW THE HECK DO I ADD A COLOR TO THE ALPHA LEGEND

##plot prep
p0.labs <- c("p0: 0.05", "p0: 0.20")
names(p0.labs) <- c("0.05", "0.2")

coh.labs <- c("Cohesion: 0","Cohesion: 0.3", "Cohesion: 0.67", "Cohesion: 1")
names(coh.labs) <- c("0","0.3", "0.67", "1")
Truth<-"2"


plot_Nmed_SC<-ggplot(data=SCresults[SCresults$param=="N",],
                     aes(x=aggregation, y=Median,fill=aggregation))+# alpha=cohesion,
  geom_boxplot()+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  geom_hline(aes(yintercept=140,linetype="Truth"),color="black")+ #,  color = "N.inds"
  scale_linetype_manual(name="",values=2)+
  guides(color=guide_legend(order=2),
         alpha=guide_legend(order=3),
         linetype=guide_legend(order=1))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Abundance (N)",x="",y="Estimate")+theme_bw()+ #Aggregation (Group Size)
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_Nmed_SC_alt<-ggplot(SCresults[SCresults$param=="N",],
       aes(x=aggregation, y=Median,color=cohesion)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Estimate")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

# ggplot(SCresults[SCresults$param=="N",],
#        aes(x=aggregation, y=Median)) + 
#   geom_boxplot(aes(fill=p0))+
#   # scale_color_brewer(palette="Blues")+
#   geom_hline(yintercept=141,lty=2)+
#   labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Estimate")+
#   facet_grid(vars(cohesion),switch="y")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_Nmean_SC<-ggplot(data=SCresults[SCresults$param=="N",],
                     aes(x=aggregation, y=Mean,fill=aggregation))+# alpha=cohesion,
  geom_boxplot()+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  geom_hline(aes(yintercept=140,linetype="Truth"),color="black")+ #,  color = "N.inds"
  scale_linetype_manual(name="",values=2)+
  guides(color=guide_legend(order=2),
         alpha=guide_legend(order=3),
         linetype=guide_legend(order=1))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Abundance (N)",x="",y="Estimate")+theme_bw()+ #Aggregation (Group Size)
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_Nmean_SC_alt<-ggplot(SCresults[SCresults$param=="N",],
                         aes(x=aggregation, y=Mean,color=cohesion)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Estimate")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_N_rb_SC<-ggplot(data=SCresults[SCresults$param=="N",])+
  geom_boxplot(aes(x=aggregation, y=RB,alpha=cohesion,fill=aggregation))+
  geom_hline(yintercept=0,lty=2)+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
   facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="RB")+theme_bw()+ #Aggregation (Group Size),title="Abundance (N)",title="Abundance (N)",
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #strip.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_N_rb_SC_alt<-ggplot(SCresults[SCresults$param=="N",],
       aes(x=aggregation, y=RB,color=cohesion)) + 
  geom_boxplot()+
  geom_hline(yintercept=0,lty=2)+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size",y="RB")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_N_cv_SC<-ggplot(data=SCresults[SCresults$param=="N",])+
  geom_boxplot(aes(x=aggregation, y=CoV,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  geom_hline(aes(yintercept=0),linetype=2,color="black")+ #,  color = "N.inds"
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="CoV")+theme_bw()+ #title="Abundance (N)",
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
       # strip.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_N_cv_SC_alt<-ggplot(SCresults[SCresults$param=="N",],
       aes(x=aggregation, y=CoV,color=cohesion)) + 
  geom_hline(aes(yintercept=0),linetype=2,color="black")+ #,  color = "N.inds"
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="CoV")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_N_rmse_SC<-ggplot(data=SCresults[SCresults$param=="N",])+
  geom_boxplot(aes(x=aggregation, y=rmse,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  geom_hline(aes(yintercept=0),linetype=2,color="black")+ #,  color = "N.inds"
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="RMSE")+theme_bw()+ #title="Abundance (N)",
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        # strip.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_N_rmse_SC_alt<-ggplot(SCresults[SCresults$param=="N",],
                         aes(x=aggregation, y=rmse,color=cohesion)) + 
  geom_hline(aes(yintercept=0),linetype=2,color="black")+ #,  color = "N.inds"
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="RMSE")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


### plot sigma ####
plot_sigMed_SC<-ggplot(data=SCresults[SCresults$param=="sigma",])+
  geom_boxplot(aes(x=aggregation, y=Median,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  geom_hline(yintercept=3, linetype="dashed", color = "black")+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Estimate")+theme_bw()+ #Aggregation (Group Size)
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


plot_sigMed_SC_alt<-ggplot(SCresults[SCresults$param=="sigma",],
       aes(x=aggregation, y=Median,color=cohesion)) + 
  geom_boxplot()+
  geom_hline(yintercept=3, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma  (\u03c3)", x="Aggregation (Group Size)",y="CoV")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_sigMean_SC<-ggplot(data=SCresults[SCresults$param=="sigma",])+
  geom_boxplot(aes(x=aggregation, y=Mean,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  geom_hline(yintercept=3, linetype="dashed", color = "black")+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Estimate")+theme_bw()+ #Aggregation (Group Size)
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


plot_sigMean_SC_alt<-ggplot(SCresults[SCresults$param=="sigma",],
                           aes(x=aggregation, y=Mean,color=cohesion)) + 
  geom_boxplot()+
  geom_hline(yintercept=3, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma  (\u03c3)", x="Aggregation (Group Size)",y="CoV")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_sig_rb_SC<-ggplot(data=SCresults[SCresults$param=="sigma",])+
  geom_boxplot(aes(x=aggregation, y=RB,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+theme_bw()+ #title="Sigma (\u03c3)", Aggregation (Group Size)
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


plot_sig_rb_SC_alt<-ggplot(SCresults[SCresults$param=="sigma",],
       aes(x=aggregation, y=RB,color=cohesion)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sig_cv_SC<-ggplot(data=SCresults[SCresults$param=="sigma",])+
  geom_boxplot(aes(x=aggregation, y=CoV,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="CoV")+theme_bw()+ #title="Sigma (\u03c3)", 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sig_cv_SC_alt<-ggplot(SCresults[SCresults$param=="sigma",],
       aes(x=aggregation, y=CoV,color=cohesion)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_sig_cv_SC<-ggplot(data=SCresults[SCresults$param=="sigma",])+
  geom_boxplot(aes(x=aggregation, y=rmse,alpha=cohesion,fill=aggregation))+
  scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  scale_alpha_manual(name = "Cohesion",
                     labels=c("0","0.3","0.67","1"),
                     values=c(0.2,0.4,0.6,0.8))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), 
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+ #switch = "y"
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="RMSE")+theme_bw()+ #title="Sigma (\u03c3)", 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sig_rmse_SC_alt<-ggplot(SCresults[SCresults$param=="sigma",],
                           aes(x=aggregation, y=rmse,color=cohesion)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RMSE")+
  facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

#plot_N
#plot_sig

# ggarrange(plot_N + rremove("x.text"), plot_sig + rremove("x.text"),
#           plot_N_rb + rremove("x.text"),plot_sig_rb + rremove("x.text"),
#           plot_N_cov,plot_sig_cov, 
#           heights=c(1.3,1,1),
#           common.legend = TRUE,legend="right",
#           ncol = 2, nrow = 3)


### coverage ####

#temporarily manipulate the coverage data so that they dont overlap when plotted
coverages_boot_SC$aggregation<-as.numeric(levels(coverages_boot_SC$aggregation))[coverages_boot_SC$aggregation]
coverages_calc_SC$aggregation<-as.numeric(levels(coverages_calc_SC$aggregation))[coverages_calc_SC$aggregation]

plot_N_coverage_SC<-ggplot() + # 
  #data=coverages_boot_SC[coverages_boot_SC$param=="N",],
  #aes(x=aggregation, y=coverage,group=cohesion)
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_line(data=coverages_calc_SC[coverages_calc_SC$param=="N",],
            aes(x=aggregation, y=coverage, group=cohesion))+#,color=cohesion))+
  geom_point(data=coverages_calc_SC[coverages_calc_SC$param=="N",], size=3,
             aes(x=aggregation, y=coverage,group=cohesion),pch=21,fill="white")+#shape=cohesion,color=cohesion))+
  # stat_summary(fun.min = function(x) min(x), 
  #              fun.max = function(x) max(x), 
  #              geom = "linerange",size=1,
  #              aes(color=cohesion),show.legend = FALSE) +
  # scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",
       y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


plot_sigma_coverage_SC<-ggplot() + # 
  #data=coverages_boot_SC[coverages_boot_SC$param=="sigma",], aes(x=aggregation, y=coverage,group=cohesion)
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_line(data=coverages_calc_SC[coverages_calc_SC$param=="sigma",],
            aes(x=aggregation, y=coverage, group=cohesion))+#,color=cohesion))+
  geom_point(data=coverages_calc_SC[coverages_calc_SC$param=="sigma",], size=3,
             pch=21,fill="white",
             aes(x=aggregation, y=coverage,group=cohesion))+#, shape=cohesion,color=cohesion))+
  # stat_summary(fun.min = function(x) min(x), 
  #              fun.max = function(x) max(x), 
  #              geom = "linerange",size=1,
  #              aes(color=cohesion),show.legend = FALSE) +
  #scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",
       y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
#wow