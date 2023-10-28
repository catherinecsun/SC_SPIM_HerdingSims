nits <-c(10,25,50,75,80,85,90,95,100)

scenarios<-unique(SCRresults$scenario)


###SCR ####
load("workspace_script3_processingSCRresults.RData")
colnames(SCRresults)


#first one
whichIts<-sample(1:100,nits[1])
tempy<-SCRresults[SCRresults$param=="N"&
                    SCRresults$sim%in%whichIts,]
tempy$nits<-nits[1]

for(n in 2:length(nits)){
  print(nits[n])
  #get the next random batch of simulations
  whichIts<-sample(1:100,nits[n])

  tempy_2<-SCRresults[SCRresults$param=="N"&
                      SCRresults$sim%in%whichIts,]
  tempy_2$nits<-nits[n]
  
  tempy<-rbind(tempy,tempy_2)
}


# plot pattern in median estimates (hopefully see a tapering)

colGrayfunc <- colorRampPalette(c("white","darkgray"))
colGrayfunc(length(nits))
palette <-colGrayfunc(length(nits))

plot_SCR_Nmed_axNits<-ggplot(tempy,
                   aes(x=aggregation, y=Median,fill=as.factor(nits))) + 
  geom_boxplot(position=position_dodge(1))+
  scale_fill_manual(values=palette)+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=141,lty=2)+
  labs(title="SCR Median Abundance Estimates, with Increasing No. of Iterations",
       x="Aggregation (Group Size)",y="N",
       fill="Iterations")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             scales = "free_x", space = "free"#,
             #labeller = labeller(p0 = p0.labs,cohesion=coh.labs)
             )+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#bias
plot_SCR_Nrb_axNits<-ggplot(tempy,
                  aes(x=aggregation, y=RB,fill=as.factor(nits))) + 
  scale_fill_manual(values=palette)+
  geom_boxplot(position=position_dodge(1))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  # ylim(-2,30)+
  labs(title="Relative Bias of SCR Estimates, with Increasing No. of Iterations",
       x="Aggregation (Group Size)",y="RB",
       fill="Iterations")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",  scales = "free_x", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,



plot_SCR_N_CV_axNits<-ggplot(tempy,
                  aes(x=aggregation, y=CoV,
                      fill=as.factor(nits))) + 
  scale_fill_manual(values=palette)+
  geom_boxplot(position=position_dodge(1))+
  # geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Coefficient Variation of SCR Estimates, with Increasing No. of Iterations", 
       x="Aggregation (Group Size)",y="CV",
       fill="Iterations")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


###SC ####
load("workspace_script3_processingSCresults.RData")
colnames(SCresults)


#first one
whichIts<-sample(1:100,nits[1])
tempy<-SCresults[SCresults$param=="N"&
                   SCresults$sim%in%whichIts,]
tempy$nits<-nits[1]

for(n in 2:length(nits)){
  print(nits[n])
  #get the next random batch of simulations
  whichIts<-sample(1:100,nits[n])
  
  tempy_2<-SCresults[SCresults$param=="N"&
                       SCresults$sim%in%whichIts,]
  tempy_2$nits<-nits[n]
  
  tempy<-rbind(tempy,tempy_2)
}


# plot pattern in median estimates (hopefully see a tapering)

colGrayfunc <- colorRampPalette(c("white","darkgray"))
colGrayfunc(length(nits))
palette <-colGrayfunc(length(nits))

plot_SC_Nmed_axNits<-ggplot(tempy,
                             aes(x=aggregation, y=Median,fill=as.factor(nits))) + 
  geom_boxplot(position=position_dodge(1))+
  scale_fill_manual(values=palette)+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=141,lty=2)+
  labs(title="SC Median Abundance Estimates, with Increasing No. of Iterations",
       x="Aggregation (Group Size)",y="N",
       fill="Iterations")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
             scales = "free_x", space = "free"#,
             #labeller = labeller(p0 = p0.labs,cohesion=coh.labs)
  )+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,


#bias
plot_SC_Nrb_axNits<-ggplot(tempy,
                            aes(x=aggregation, y=RB,fill=as.factor(nits))) + 
  scale_fill_manual(values=palette)+
  geom_boxplot(position=position_dodge(1))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  # ylim(-2,30)+
  labs(title="Relative Bias of SC Estimates, with Increasing No. of Iterations",
       x="Aggregation (Group Size)",y="RB",
       fill="Iterations")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",  scales = "free_x", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,



plot_SC_N_CV_axNits<-ggplot(tempy,
                             aes(x=aggregation, y=CoV,
                                 fill=as.factor(nits))) + 
  scale_fill_manual(values=palette)+
  geom_boxplot(position=position_dodge(1))+
  # geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(col="gray",xintercept = c(1.5,2.5))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
  # scale_color_brewer(palette="Blues")+
  labs(title="Coefficient Variation of SC Estimates, with Increasing No. of Iterations", 
       x="Aggregation (Group Size)",y="CV",
       fill="Iterations")+
  facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free", space = "free",
             labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,









###Select SPIMs ####
load("workspace_script3_processingSPIMresults.RData")
colnames(SPIMresults)


#first one
whichIts<-sample(1:100,nits[1])
tempy<-SPIMresults[SPIMresults$param=="N"&
                   SPIMresults$sim%in%whichIts,]
tempy$nits<-nits[1]

for(n in 2:length(nits)){
  print(nits[n])
  #get the next random batch of simulations
  whichIts<-sample(1:100,nits[n])
  
  tempy_2<-SPIMresults[SPIMresults$param=="N"&
                         SPIMresults$sim%in%whichIts,]
  tempy_2$nits<-nits[n]
  
  tempy<-rbind(tempy,tempy_2)
}

table(tempy$PID)

# plot pattern in median estimates (hopefully see a tapering)

colGrayfunc <- colorRampPalette(c("white","darkgray"))
colGrayfunc(length(nits))
palette <-colGrayfunc(length(nits))

for(p in 1:length(unique(tempy$PID))){
  print(as.character(unique(tempy$PID)[p]))
  plot_SPIM_Nmed_axNits<-ggplot(tempy[tempy$PID==unique(tempy$PID)[p],],
                              aes(x=aggregation, y=Median,fill=as.factor(nits))) + 
    geom_boxplot(position=position_dodge(1))+
    scale_fill_manual(values=palette)+
    geom_vline(col="gray",xintercept = c(1.5,2.5))+
    geom_hline(yintercept=141,lty=2)+
    labs(title=paste0("SPIM ",unique(tempy$PID)[p], " Median Abundance Estimates, with Increasing No. of Iterations"),
         x="Aggregation (Group Size)",y="N",
         fill="Iterations")+
    facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",
               scales = "free_x", space = "free"#,
               #labeller = labeller(p0 = p0.labs,cohesion=coh.labs)
    )+
    theme_bw()+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  assign(paste0("plot_SPIM_",as.character(unique(tempy$PID)[p]),"_Nmed_axNits"),
         plot_SPIM_Nmed_axNits)
  
  
  #bias
  plot_SPIM_Nrb_axNits<-ggplot(tempy[tempy$PID==unique(tempy$PID)[p],],
                             aes(x=aggregation, y=RB,fill=as.factor(nits))) + 
    scale_fill_manual(values=palette)+
    geom_boxplot(position=position_dodge(1))+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_vline(col="gray",xintercept = c(1.5,2.5))+
    # ylim(-2,30)+
    labs(title=paste0("Relative Bias of SPIM ",unique(tempy$PID)[p], " Estimates, with Increasing No. of Iterations"),
         x="Aggregation (Group Size)",y="RB",
         fill="Iterations")+
    facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",  scales = "free_x", space = "free",
               labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
    theme_bw()+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  
  assign(paste0("plot_SPIM_",as.character(unique(tempy$PID)[p]),"_Nrb_axNits"),
         plot_SPIM_Nrb_axNits)
  
  
  
  #cv
  plot_SPIM_N_CV_axNits<-ggplot(tempy[tempy$PID==unique(tempy$PID)[p],],
                              aes(x=aggregation, y=CoV,
                                  fill=as.factor(nits))) + 
    scale_fill_manual(values=palette)+
    geom_boxplot(position=position_dodge(1))+
    # geom_hline(yintercept=0, linetype="dashed", color = "black")+
    geom_vline(col="gray",xintercept = c(1.5,2.5))+
    geom_hline(yintercept=0.2, linetype="dashed", color = "black")+
    # scale_color_brewer(palette="Blues")+
    labs(title=paste0("Coefficient Variation of SPIM ",unique(tempy$PID)[p]," Estimates, with Increasing No. of Iterations"), 
         x="Aggregation (Group Size)",y="CV",
         fill="Iterations")+
    facet_grid(rows=vars(p0),cols=vars(cohesion),switch="y",scales = "free", space = "free",
               labeller = labeller(p0 = p0.labs,cohesion=coh.labs))+
    theme_bw()+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())#,
  
  assign(paste0("plot_SPIM_",as.character(unique(tempy$PID)[p]),"_NCV_axNits"),
         plot_SPIM_N_CV_axNits)
  
}

#### all plots ####
plot_SCR_Nmed_axNits
plot_SCR_Nrb_axNits
plot_SCR_N_CV_axNits

plot_SC_Nmed_axNits
plot_SC_Nrb_axNits
plot_SC_N_CV_axNits

plot_SPIM_collar_Nmed_axNits
plot_SPIM_collar_Nrb_axNits
plot_SPIM_collar_NCV_axNits

plot_SPIM_sex_Nmed_axNits
plot_SPIM_sex_Nrb_axNits
plot_SPIM_sex_NCV_axNits

plot_SPIM_sexcollar_Nmed_axNits
plot_SPIM_sexcollar_Nrb_axNits
plot_SPIM_sexcollar_NCV_axNits

plot_SPIM_sexcoat_Nmed_axNits
plot_SPIM_sexcoat_Nrb_axNits
plot_SPIM_sexcoat_NCV_axNits

plot_SPIM_sexcollarcoat_Nmed_axNits
plot_SPIM_sexcollarcoat_Nrb_axNits
plot_SPIM_sexcollarcoat_NCV_axNits

plot_SPIM_antlerssexcollarcoat_Nmed_axNits
plot_SPIM_antlerssexcollarcoat_Nrb_axNits
plot_SPIM_antlerssexcollarcoat_NCV_axNits
