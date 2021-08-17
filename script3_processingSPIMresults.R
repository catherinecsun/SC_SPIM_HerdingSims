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
filePattern<-"SPIM_simResults_scenario"

results_SPIM<-list.files("SPIMresults") #the folder that spim results are in
results_SPIM<-results_SPIM[grepl(filePattern,results_SPIM)]

whichPartialIDS<-list(paste(c("_","sex","collar","coat"),collapse=""),
                      paste(c("_","antlers","sex","collar","coat"),collapse=""))

#empty dataframes to put things in
SPIMresults<-data.frame()
coverages_calc_SPIM<-data.frame()
coverages_boot_SPIM<-data.frame()


# some plotty prep
p0.labs <- c("p0: 0.05", "p0: 0.20")
names(p0.labs) <- c("0.05", "0.2")

coh.labs <- c("Cohesion: 0","Cohesion: 0.3", "Cohesion: 0.67", "Cohesion: 1")
names(coh.labs) <- c("0","0.3", "0.67", "1")
Truth<-"2" # for line type


#loop through each batch of partial ids
#bring in the data, calculate stats, and plot

for(pids in 1:length(whichPartialIDS)){
  cat(paste0(pids,":",whichPartialIDS[[pids]]))
  
  results_SPIM_subset<-results_SPIM[grepl(whichPartialIDS[pids],results_SPIM)]
  results_SPIM_subset
  
  SPIMresults_allDF<-data.frame()
  
  for(i in 1:length(results_SPIM_subset)){
    tmp<-read.csv(paste0("SPIMresults/",results_SPIM_subset[i]))
    SPIMresults_allDF <-rbind.fill(SPIMresults_allDF,tmp)
  }
  
  #remove garbage first column
  SPIMresults_allDF<-SPIMresults_allDF[,-1]
  
  #double check things
  table(SPIMresults_allDF$scenario)
  
  #remove any sims that didnt converge
  which(SPIMresults_allDF$Rhat_pt>1.2)
  convergeFail<-unique(SPIMresults_allDF[which(SPIMresults_allDF$Rhat_pt>1.2),colnames(SPIMresults_allDF)%in%c("scenario","sim")])
  convergeFail_remove<-c()
  if(nrow(convergeFail)>0){
    for(f in 1:nrow(convergeFail)){
      convergeFail_remove<-c(convergeFail_remove,
                             intersect(which(SPIMresults_allDF$scenario==convergeFail$scenario[f]),
                                       which(SPIMresults_allDF$sim==convergeFail$sim[f])))
      
      #btw fwiw i think match_df function in plyr would have worked for this too
    } 
    SPIMresults_allDF<-SPIMresults_allDF[-convergeFail_remove,]
  }
 
  #how many sims remain
  table(SPIMresults_allDF$scenario)
  
  # add columns for the parameter/scenario conditions
  #use a different version of parm_combos that focuses 
  #only on the 18 combos that are different (and that we evaluated)
  parm_combos2<-parm_combos[-c(3:8),]
  
  SPIMresults_allDF$p0<-NA
  SPIMresults_allDF$cohesion<-NA
  SPIMresults_allDF$aggregation<-NA
  for(r in 1:nrow(parm_combos2)){
    SPIMresults_allDF$p0[SPIMresults_allDF$scenario==rownames(parm_combos2)[r]]<-parm_combos2$p0[r]
    SPIMresults_allDF$cohesion[SPIMresults_allDF$scenario==rownames(parm_combos2)[r]]<-parm_combos2$cohesion[r]
    SPIMresults_allDF$aggregation[SPIMresults_allDF$scenario==rownames(parm_combos2)[r]]<-parm_combos2$aggregation[r]
  }
  
  SPIMresults_allDF$cohesion<-as.factor(SPIMresults_allDF$cohesion)
  SPIMresults_allDF$p0<-as.factor(SPIMresults_allDF$p0)
  SPIMresults_allDF$aggregation<-as.factor(SPIMresults_allDF$aggregation)
  
  
  ### relative bias  ####
  #(truth-estimate)/estimate
  area<-(range(sim_data[[1]][[1]]$mask[,1])[2]-range(sim_data[[1]][[1]]$mask[,1])[1])*
    (range(sim_data[[1]][[1]]$mask[,2])[2]-range(sim_data[[1]][[1]]$mask[,2])[1])
  truth<-c("D"=N.inds/area,"N"=N.inds,"psi"=N.inds/M,"sigma"=3) #"lam0"=c(0.05,0.20),
  
  SPIMresults_allDF$Truth[which(SPIMresults_allDF$param%in%names(truth))]<-truth[match(SPIMresults_allDF$param[which(SPIMresults_allDF$param%in%names(truth))],names(truth))]
  SPIMresults_allDF$Truth[which(SPIMresults_allDF$param=="lam0"&SPIMresults_allDF$scenario %% 2 == 0)]<-0.2
  SPIMresults_allDF$Truth[which(SPIMresults_allDF$param=="lam0"&SPIMresults_allDF$scenario %% 2 == 1)]<-0.05
  
  SPIMresults_allDF$RB<-(SPIMresults_allDF$Truth-SPIMresults_allDF$Mean)/SPIMresults_allDF$Mean
  
  
  ### Coefficient of Variation, precision ####
  #SD/estimate
  SPIMresults_allDF$CoV<-SPIMresults_allDF$SD/SPIMresults_allDF$Mean
  
  
  ##rmse
  SPIMresults_allDF$rmse<-NA
  for(n in 1:length(unique(SPIMresults_allDF$scenario))){
    SPIMresults_allDF$rmse[SPIMresults_allDF$scenario==n&SPIMresults_allDF$param=="D"]<-
      rmse(truth[names(truth)=="D"], SPIMresults_allDF$Mean[SPIMresults_allDF$scenario==n&SPIMresults_allDF$param=="D"])
    SPIMresults_allDF$rmse[SPIMresults_allDF$scenario==n&SPIMresults_allDF$param=="N"]<-
      rmse(truth[names(truth)=="N"], SPIMresults_allDF$Mean[SPIMresults_allDF$scenario==n&SPIMresults_allDF$param=="N"])
    SPIMresults_allDF$rmse[SPIMresults_allDF$scenario==n&SPIMresults_allDF$param=="sigma"]<-
      rmse(truth[names(truth)=="sigma"], SPIMresults_allDF$Mean[SPIMresults_allDF$scenario==n&SPIMresults_allDF$param=="sigma"])
  }
  
  #coverage
  SPIMresults_allDF$coverage<-NA
  for(co in 1:nrow(SPIMresults_allDF)){
    if(SPIMresults_allDF$param[co]=="D"){
      SPIMresults_allDF$coverage[co]<-ifelse(truth[names(truth)=="D"]>SPIMresults_allDF$X2.5.[co]&SPIMresults_allDF$X97.5.[co]>truth[names(truth)=="D"],
                                             1,0)
    }
    if(SPIMresults_allDF$param[co]=="N"){
      SPIMresults_allDF$coverage[co]<-ifelse(truth[names(truth)=="N"]>SPIMresults_allDF$X2.5.[co]&SPIMresults_allDF$X97.5.[co]>truth[names(truth)=="N"],
                                             1,0)
    }
    if(SPIMresults_allDF$param[co]=="sigma"){
      SPIMresults_allDF$coverage[co]<-ifelse(truth[names(truth)=="sigma"]>SPIMresults_allDF$X2.5.[co]&SPIMresults_allDF$X97.5.[co]>truth[names(truth)=="sigma"],
                                             1,0)
    }
  }
  
  SPIMresults_allDF$PID<-substring(whichPartialIDS[[pids]],2,nchar(whichPartialIDS[[pids]]))
  SPIMresults_allDF$antler<-ifelse(grepl("antler",whichPartialIDS[[pids]]),1,0)
  SPIMresults_allDF$sex<-ifelse(grepl("sex",whichPartialIDS[[pids]]),1,0)
  SPIMresults_allDF$collar<-ifelse(grepl("collar",whichPartialIDS[[pids]]),1,0)
  SPIMresults_allDF$coat<-ifelse(grepl("coat",whichPartialIDS[[pids]]),1,0)
  SPIMresults <-rbind.fill(SPIMresults,SPIMresults_allDF)
  
  library(dplyr)
  
  # caluclating coverage
  coverage_SPIM_calc<-SPIMresults_allDF[SPIMresults_allDF$param=="N"|SPIMresults_allDF$param=="sigma",]%>%group_by(cohesion, aggregation, p0,param)%>%summarise(coverage=sum(coverage)/100)
  
  coverage_SPIM_calc$PID<-substring(whichPartialIDS[[pids]],2,nchar(whichPartialIDS[[pids]]))
  coverage_SPIM_calc$antler<-ifelse(grepl("antler",whichPartialIDS[[pids]]),1,0)
  coverage_SPIM_calc$sex<-ifelse(grepl("sex",whichPartialIDS[[pids]]),1,0)
  coverage_SPIM_calc$collar<-ifelse(grepl("collar",whichPartialIDS[[pids]]),1,0)
  coverage_SPIM_calc$coat<-ifelse(grepl("coat",whichPartialIDS[[pids]]),1,0)
  coverages_calc_SPIM<-rbind(coverages_calc_SPIM,coverage_SPIM_calc)
  
  #bootstrapping coverage
  coverage_SPIM<-SPIMresults_allDF[SPIMresults_allDF$param=="N"|SPIMresults_allDF$param=="sigma",]%>%group_by(cohesion, aggregation, p0,param)%>% summarise(coverage=booties(coverage))
   
  coverage_SPIM$PID<-substring(whichPartialIDS[[pids]],2,nchar(whichPartialIDS[[pids]]))
  coverage_SPIM$antler<-ifelse(grepl("antler",whichPartialIDS[[pids]]),1,0)
  coverage_SPIM$sex<-ifelse(grepl("sex",whichPartialIDS[[pids]]),1,0)
  coverage_SPIM$collar<-ifelse(grepl("collar",whichPartialIDS[[pids]]),1,0)
  coverage_SPIM$coat<-ifelse(grepl("coat",whichPartialIDS[[pids]]),1,0)
  coverages_boot_SPIM <-rbind(coverages_boot_SPIM,coverage_SPIM)
  
  
  ### relative variance ####
  
  
  ### plot N ####
  #   HOW THE HECK DO I ADD A COLOR TO THE ALPHA LEGEND
  
  plot_N<-ggplot(data=SPIMresults_allDF[SPIMresults_allDF$param=="N",])+
    geom_boxplot(aes(x=aggregation, y=X50.,alpha=cohesion,fill=aggregation))+
    scale_fill_discrete(name = "Group Size", labels = c('1', '4','10'))+
    scale_alpha_manual(name = "Cohesion",
                       labels=c("0","0.3","0.67","1"),
                       values=c(0.2,0.4,0.6,0.8))+
    geom_hline(aes(yintercept=140,linetype="Truth"),color="black")+ #,  color = "N.inds"
    scale_linetype_manual(name="",values=2)+
    #geom_hline(yintercept=140/4, linetype="dashed", color = "black")+
    #geom_hline(yintercept=140/10, linetype="dashed", color = "black")+
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
          #strip.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  plot_N_alt<-ggplot(SPIMresults_allDF[SPIMresults_allDF$param=="N",],
                     aes(x=aggregation, y=X50.,color=cohesion)) + 
    geom_boxplot()+
    # scale_color_brewer(palette="Blues")+
    geom_hline(yintercept=141,lty=2)+
    labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Estimate")+
    facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
 
  plot_N_rb<-ggplot(data=SPIMresults_allDF[SPIMresults_allDF$param=="N",])+
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
  
 
  plot_N_rb_alt<-ggplot(SPIMresults_allDF[SPIMresults_allDF$param=="N",],
                        aes(x=aggregation, y=RB,color=cohesion)) + 
    geom_boxplot()+
    geom_hline(yintercept=0,lty=2)+
    # scale_color_brewer(palette="Blues")+
    labs(title="Abundance (N)", x="Aggregation (Group Size",y="RB")+
    facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  
 plot_N_cv<-ggplot(data=SPIMresults_allDF[SPIMresults_allDF$param=="N",])+
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
  
plot_N_cv_alt<-ggplot(SPIMresults_allDF[SPIMresults_allDF$param=="N",],
                        aes(x=aggregation, y=CoV,color=cohesion)) + 
    geom_hline(aes(yintercept=0),linetype=2,color="black")+ #,  color = "N.inds"
    geom_boxplot()+
    # scale_color_brewer(palette="Blues")+
    labs(title="Abundance (N)", x="Aggregation (Group Size)",y="CoV")+
    facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  
  
  ### plot sigma ####
  plot_sig<-ggplot(data=SPIMresults_allDF[SPIMresults_allDF$param=="sigma",])+
    geom_boxplot(aes(x=aggregation, y=X50.,alpha=cohesion,fill=aggregation))+
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
  
  plot_sig_alt<-ggplot(SPIMresults_allDF[SPIMresults_allDF$param=="sigma",],
                       aes(x=aggregation, y=X50.,color=cohesion)) + 
    geom_boxplot()+
    geom_hline(yintercept=3, linetype="dashed", color = "black")+
    # scale_color_brewer(palette="Blues")+
    labs(title="Sigma  (\u03c3)", x="Aggregation (Group Size)",y="CoV")+
    facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  
  
 
  plot_sig_rb<-ggplot(data=SPIMresults_allDF[SPIMresults_allDF$param=="sigma",])+
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
  
  
  plot_sig_rb_alt<-ggplot(SPIMresults_allDF[SPIMresults_allDF$param=="sigma",],
                          aes(x=aggregation, y=RB,color=cohesion)) + 
    geom_boxplot()+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    # scale_color_brewer(palette="Blues")+
    labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+
    facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  
 
  plot_sig_cv<-ggplot(data=SPIMresults_allDF[SPIMresults_allDF$param=="sigma",])+
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
  
 
  plot_sig_cv_alt<-ggplot(SPIMresults_allDF[SPIMresults_allDF$param=="sigma",],
                           aes(x=aggregation, y=CoV,color=cohesion)) + 
    geom_boxplot()+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    # scale_color_brewer(palette="Blues")+
    labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+
    facet_grid(vars(p0),switch="y",labeller = labeller(p0 = p0.labs))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  

 
  ### coverage ####
  
  #temporarily manipulate the coverage data so that they dont overlap when plotted
  coverage_SPIM$aggregation<-as.numeric(levels(coverage_SPIM$aggregation))[coverage_SPIM$aggregation]
  coverage_SPIM$aggregation[coverage_SPIM$cohesion==0]<-coverage_SPIM$aggregation[coverage_SPIM$cohesion==0]-0.6
  coverage_SPIM$aggregation[coverage_SPIM$cohesion==0.3]<-coverage_SPIM$aggregation[coverage_SPIM$cohesion==0.3]-0.3
  coverage_SPIM$aggregation[coverage_SPIM$cohesion==1]<-coverage_SPIM$aggregation[coverage_SPIM$cohesion==1]+0.3
  
  coverage_SPIM_calc$aggregation<-as.numeric(levels(coverage_SPIM_calc$aggregation))[coverage_SPIM_calc$aggregation]
  coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0]<-coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0]-0.6
  coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0.3]<-coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0.3]-0.3
  coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==1]<-coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==1]+0.3
  
  plot_coverage<- ggplot(data=coverage_SPIM, aes(x=aggregation, y=coverage,group=cohesion)) + # 
    geom_point(data=coverage_SPIM_calc, size=3,
               aes(x=aggregation, y=coverage,
                   shape=cohesion,group=cohesion,color=cohesion))+
    stat_summary(fun.min = function(x) min(x), 
                 fun.max = function(x) max(x), 
                 geom = "linerange",size=1,
                 aes(color=cohesion),show.legend = FALSE) +
    # stat_summary(fun = mean,
    #              geom = "line",aes(color=cohesion)) +
    scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
    facet_grid(rows=vars(p0),cols=vars(param),labeller = labeller(p0 = p0.labs))+
    labs(title="Abundance (N)       and Sigma (\u03c3)", x="Aggregation (Group Size)",
         y="Coverage")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
 
  #return coverage data so aggregation values are correct. 
  
  coverage_SPIM$aggregation[coverage_SPIM$cohesion==0]<-coverage_SPIM$aggregation[coverage_SPIM$cohesion==0]+0.6
  coverage_SPIM$aggregation[coverage_SPIM$cohesion==0.3]<-coverage_SPIM$aggregation[coverage_SPIM$cohesion==0.3]+0.3
  coverage_SPIM$aggregation[coverage_SPIM$cohesion==1]<-coverage_SPIM$aggregation[coverage_SPIM$cohesion==1]-0.3
  coverage_SPIM$aggregation<-as.factor(coverage_SPIM$aggregation)
  
  coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0]<-coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0]+0.6
  coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0.3]<-coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==0.3]+0.3
  coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==1]<-coverage_SPIM_calc$aggregation[coverage_SPIM_calc$cohesion==1]-0.3
  coverage_SPIM_calc$aggregation<-as.factor(coverage_SPIM_calc$aggregation)
  
    
  ### the plots, giving them more specific names now
 assign(paste0("plot_Nmed",whichPartialIDS[[pids]]),plot_N) 
 assign(paste0("plot_Nmed",whichPartialIDS[[pids]],"_alt"),plot_N_alt)
 assign(paste0("plot_N_rb",whichPartialIDS[[pids]]),plot_N_rb) 
 assign(paste0("plot_N_rb",whichPartialIDS[[pids]],"_alt"),plot_N_rb_alt)
 assign(paste0("plot_N_cv",whichPartialIDS[[pids]]),plot_N_cv) 
 assign(paste0("plot_N_cv",whichPartialIDS[[pids]],"_alt"),plot_N_cv_alt)
 
 assign(paste0("plot_sigMed",whichPartialIDS[[pids]]),plot_sig) 
 assign(paste0("plot_sigMed",whichPartialIDS[[pids]],"_alt"),plot_sig_alt)
 assign(paste0("plot_sig_rb",whichPartialIDS[[pids]]),plot_sig_rb) 
 assign(paste0("plot_sig_rb",whichPartialIDS[[pids]],"_alt"),plot_sig_rb_alt)
 assign(paste0("plot_sig_cv",whichPartialIDS[[pids]]),plot_N_cv) 
 assign(paste0("plot_sig_cv",whichPartialIDS[[pids]],"_alt"),plot_sig_cv_alt)
 
 assign(paste0("plot_coverage",whichPartialIDS[[pids]]),plot_coverage) 
  
}
rm(plot_N,plot_N_alt,plot_N_rb,plot_N_rb_alt,plot_N_cv,plot_N_cv_alt,
   plot_sig,plot_sig_alt,plot_sig_rb,plot_sig_rb_alt,plot_sig_cv,plot_sig_cv_alt,
   plot_coverage)



#### all scenarios all together now

#median and means
plot_Nmed_SPIM<-ggplot(SPIMresults[SPIMresults$param=="N",],
       aes(x=aggregation, y=X50.,color=PID)) + 
      geom_boxplot()+
      # scale_color_brewer(palette="Blues")+
     geom_hline(yintercept=141,lty=2)+
      labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Median Estimate")+
      facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
                 labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
     theme_bw()+
     theme(plot.title = element_text(hjust = 0.5),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
 
plot_Nmean_SPIM<-ggplot(SPIMresults[SPIMresults$param=="N",],
                       aes(x=aggregation, y=Mean,color=PID)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=141,lty=2)+
  labs(title="Abundance (N)",x="Aggregation (Group Size)",y="Median Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


plot_sigmaMed_SPIM<-ggplot(SPIMresults[SPIMresults$param=="sigma",],
                    aes(x=aggregation, y=X50.,color=PID)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=3,lty=2)+
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Mean Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigmaMean_SPIM<-ggplot(SPIMresults[SPIMresults$param=="sigma",],
                           aes(x=aggregation, y=Mean,color=PID)) + 
  geom_boxplot()+
  # scale_color_brewer(palette="Blues")+
  geom_hline(yintercept=3,lty=2)+
  labs(title="Sigma (\u03c3)",x="Aggregation (Group Size)",y="Mean Estimate")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

#rb
plot_N_rb_SPIM<-ggplot(SPIMresults[SPIMresults$param=="N",],
       aes(x=aggregation, y=RB,color=PID)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="RB")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sigma_rb_SPIM<-ggplot(SPIMresults[SPIMresults$param=="sigma",],
                  aes(x=aggregation, y=RB,color=PID)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="RB")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

#cv

plot_N_cv_SPIM<-ggplot(SPIMresults[SPIMresults$param=="N",],
                        aes(x=aggregation, y=CoV,color=PID)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="CV")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,

plot_sig_cv_SPIM<-ggplot(SPIMresults[SPIMresults$param=="sigma",],
                        aes(x=aggregation, y=CoV,color=PID)) + 
  geom_boxplot()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="CV")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#rmse


### coverage ####
# 
# plot_N_coverage_SPIM<-ggplot(coverages_calc_SPIM[coverages_calc_SPIM$param=="N",], 
#                       aes(x=aggregation, y=coverage, group=PID)) +#
#   geom_line(aes(color=PID))+ 
#   geom_point(aes(color=PID))+
#   facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
#                                                                  cohesion=coh.labs))+
#   labs(title="Abundance (N) ", x="Aggregation (Group Size)",y="Coverage")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))
# 
# plot_sig_coverage_SPIM<-ggplot(coverages_SPIM[coverages_SPIM$param=="sigma",], 
#                              aes(x=aggregation, y=coverage, group=PID)) +#
#   geom_line(aes(color=PID))+ 
#   geom_point(aes(color=PID))+
#   facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
#                                                                    cohesion=coh.labs))+
#   labs(title="Sigma (\u03c3) ", x="Aggregation (Group Size)",y="Coverage")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))

#temporarily manipulate the coverage data so that they dont overlap when plotted
coverages_boot_SPIM$aggregation<-as.numeric(levels(coverages_boot_SPIM$aggregation))[coverages_boot_SPIM$aggregation]
coverages_boot_SPIM$aggregation[coverages_boot_SPIM$PID=="sexcollarcoat"]<-coverages_boot_SPIM$aggregation[coverages_boot_SPIM$PID=="sexcollarcoat"]-0.8

coverages_calc_SPIM$aggregation<-as.numeric(levels(coverages_calc_SPIM$aggregation))[coverages_calc_SPIM$aggregation]
coverages_calc_SPIM$aggregation[coverages_calc_SPIM$PID=="sexcollarcoat"]<-coverages_calc_SPIM$aggregation[coverages_calc_SPIM$PID=="sexcollarcoat"]-0.8


plot_N_coverage_SPIM<-ggplot(coverages_boot_SPIM[coverages_boot_SPIM$param=="N",], 
                             aes(x=aggregation, y=coverage, group=PID)) +#
  geom_point(data=coverages_calc_SPIM[coverages_calc_SPIM$param=="N",], 
             size=3,
             aes(x=aggregation, y=coverage,
                 shape=cohesion,group=cohesion,color=PID))+
  geom_line(data=coverages_calc_SPIM[coverages_calc_SPIM$param=="N",],
            aes(x=aggregation, y=coverage, group=PID,color=PID))+ 
  stat_summary(fun.min = function(x) min(x), 
               fun.max = function(x) max(x), 
               geom = "linerange",size=1,
               aes(color=PID),show.legend = FALSE) +
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),labeller=labeller(p0 = p0.labs,
                                                                 cohesion=coh.labs))+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sigma_coverage_SPIM<-ggplot(coverages_boot_SPIM[coverages_boot_SPIM$param=="sigma",], 
                             aes(x=aggregation, y=coverage, group=PID)) +#
  geom_point(data=coverages_calc_SPIM[coverages_calc_SPIM$param=="sigma",], 
             size=3,
             aes(x=aggregation, y=coverage,
                 shape=cohesion,group=cohesion,color=PID))+
  geom_line(data=coverages_calc_SPIM[coverages_calc_SPIM$param=="sigma",],
            aes(x=aggregation, y=coverage, group=PID,color=PID))+ 
  stat_summary(fun.min = function(x) min(x), 
               fun.max = function(x) max(x), 
               geom = "linerange",size=1,
               aes(color=PID),show.legend = FALSE) +
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
coverages_boot_SPIM$aggregation[coverages_boot_SPIM$PID=="sexcollarcoat"]<-coverages_boot_SPIM$aggregation[coverages_boot_SPIM$PID=="sexcollarcoat"]+0.6
coverages_boot_SPIM$aggregation<-as.factor(coverages_boot_SPIM$aggregation)

coverages_calc_SPIM$aggregation[coverages_calc_SPIM$PID=="sexcollarcoat"]<-coverages_calc_SPIM$aggregation[coverages_calc_SPIM$PID=="sexcollarcoat"]+0.6
coverages_calc_SPIM$aggregation<-as.factor(coverages_calc_SPIM$aggregation)
