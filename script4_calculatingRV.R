### relative variance ####
# the ratio of the empirical variance amount simulations for any give scenario
# and the variance associated with the independent scenario
# and correcting coverage  using chat

library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

load("workspace_script1b_calculatingOverdispersion.RData")
load("workspace_script3_processingSCRresults.RData")
load("workspace_script3_processingSCresults.RData")
load("workspace_script3_processingSPIMresults.RData")

#create a dataframe to hold chat and rv per sim of each scenario
chat_and_var<-chat_df_long
chat_and_var<-chat_and_var[-which(chat_and_var$scenario%in% c(3:8)),]
chat_and_var$sim<-rep(1:100,length(unique(chat_and_var$scenario)))
colnames(chat_and_var)

### for SCR ####
#first do for N
tmp<-SCRresults[SCRresults$param=="N",]
colnames(tmp)

chat_and_var_tmpa<-left_join(chat_and_var,
                        tmp,
                       # tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","SE.estimate"))],
                        by=c("sim","p0","cohesion","aggregation"))
colnames(chat_and_var_tmpa)[which(colnames(chat_and_var_tmpa)=="SE.estimate")]<-"N_var"

chat_and_var_tmpa$N_var<-chat_and_var_tmpa$N_var^2

#chat corrected variance etc.
chat_and_var_tmpa$N_var_chat <-chat_and_var_tmpa$N_var * chat_and_var_tmpa$chat
chat_and_var_tmpa$lcl_chat <-chat_and_var_tmpa$estimate-1.96*sqrt(chat_and_var_tmpa$N_var_chat)
chat_and_var_tmpa$ucl_chat <-chat_and_var_tmpa$estimate+1.96*sqrt(chat_and_var_tmpa$N_var_chat)
for(v in 1:nrow(chat_and_var_tmpa)){
  chat_and_var_tmpa$coverage_chat[v] <-ifelse(chat_and_var_tmpa$Truth[v]>chat_and_var_tmpa$lcl_chat[v]&
                                                chat_and_var_tmpa$ucl_chat[v]>chat_and_var_tmpa$Truth[v]
                                         ,1,0)  
}
coverages_chat_SCR_N<-chat_and_var_tmpa%>%
  group_by(cohesion, aggregation, p0)%>%
  summarise(N_coverage_chat=sum(coverage_chat,na.rm = TRUE)/100)

chat_and_var_tmpa<-chat_and_var_tmpa[,which(colnames(chat_and_var_tmpa)%in%
                                    c("chat","p0","cohesion","aggregation","sim","N_var" ))]

#then do for sigma
tmp<-SCRresults[SCRresults$param=="sigma",]
colnames(tmp)

chat_and_var_tmpb<-left_join(chat_and_var,
                        tmp,
                        #tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","SE.estimate"))],
                        by=c("sim","p0","cohesion","aggregation"))
colnames(chat_and_var_tmpb)[which(colnames(chat_and_var_tmpb)=="SE.estimate")]<-"sigma_var"
chat_and_var_tmpb$sigma_var<-chat_and_var_tmpb$sigma_var^2

#chat corrected variance etc.
chat_and_var_tmpb$sigma_var_chat <-chat_and_var_tmpb$sigma_var * chat_and_var_tmpb$chat
chat_and_var_tmpb$lcl_chat <-chat_and_var_tmpb$estimate-1.96*sqrt(chat_and_var_tmpb$sigma_var_chat)
chat_and_var_tmpb$ucl_chat <-chat_and_var_tmpb$estimate+1.96*sqrt(chat_and_var_tmpb$sigma_var_chat)
for(v in 1:nrow(chat_and_var_tmpb)){
  chat_and_var_tmpb$coverage_chat[v] <-ifelse(chat_and_var_tmpb$Truth[v]>chat_and_var_tmpb$lcl_chat[v]&
                                               chat_and_var_tmpb$ucl_chat[v]>chat_and_var_tmpb$Truth[v]
                                             ,1,0)  
}
coverages_chat_SCR_sigma<-chat_and_var_tmpb%>%
  group_by(cohesion, aggregation, p0)%>%
  summarise(sigma_coverage_chat=sum(coverage_chat,na.rm = TRUE)/100)

chat_and_var_tmpb<-chat_and_var_tmpb[,which(colnames(chat_and_var_tmpb)%in%
                                         c("chat","p0","cohesion","aggregation","sim","sigma_var" ))]

chat_and_var_SCR<-full_join(chat_and_var_tmpa,chat_and_var_tmpb)

#get mean across sims 
chat_and_var_SCR_means<-chat_and_var_SCR%>%
  group_by(p0, cohesion, aggregation)%>% #scenario, 
  summarize(chat=mean(chat),
            NVar=mean(N_var,na.rm=TRUE),SigmaVar=mean(sigma_var,na.rm=TRUE))

# add new coverage
chat_and_var_SCR_means<-left_join(chat_and_var_SCR_means,
                                  coverages_chat_SCR_N)
chat_and_var_SCR_means<-left_join(chat_and_var_SCR_means,
                                  coverages_chat_SCR_sigma)

#calculate RV for the means ax sims
chat_and_var_SCR_means$RV_N<-NA
chat_and_var_SCR_means$RV_Sigma<-NA
chat_and_var_SCR_means$RV_N[1:2]<-1
chat_and_var_SCR_means$RV_Sigma[1:2]<-1

for(i in 3:nrow(chat_and_var_SCR_means)){
  if(chat_and_var_SCR_means$p0[i]=="0.05"){
    chat_and_var_SCR_means$RV_N[i]<-chat_and_var_SCR_means$NVar[i]/chat_and_var_SCR_means$NVar[1]
    chat_and_var_SCR_means$RV_Sigma[i]<-chat_and_var_SCR_means$SigmaVar[i]/chat_and_var_SCR_means$SigmaVar[1]
  }
  else{#for p0=0.20
    chat_and_var_SCR_means$RV_N[i]<-chat_and_var_SCR_means$NVar[i]/chat_and_var_SCR_means$NVar[2]
    chat_and_var_SCR_means$RV_Sigma[i]<-chat_and_var_SCR_means$SigmaVar[i]/chat_and_var_SCR_means$SigmaVar[2]
  }
}


### for SC ####
#first do for N
tmp<-SCresults[SCresults$param=="N",]
colnames(tmp)

chat_and_var_tmpa<-left_join(chat_and_var,
                        tmp,
                        #tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","St.Dev."))],
                        by=c("sim","p0","cohesion","aggregation"))
colnames(chat_and_var_tmpa)[which(colnames(chat_and_var_tmpa)=="St.Dev.")]<-"N_var"
chat_and_var_tmpa$N_var<-chat_and_var_tmpa$N_var^2

#chat corrected variance etc.
chat_and_var_tmpa$N_var_chat <-chat_and_var_tmpa$N_var * chat_and_var_tmpa$chat
chat_and_var_tmpa$lcl_chat <-chat_and_var_tmpa$Mean-1.96*sqrt(chat_and_var_tmpa$N_var_chat)
chat_and_var_tmpa$ucl_chat <-chat_and_var_tmpa$Mean+1.96*sqrt(chat_and_var_tmpa$N_var_chat)
for(v in 1:nrow(chat_and_var_tmpa)){
  chat_and_var_tmpa$coverage_chat[v] <-ifelse(chat_and_var_tmpa$Truth[v]>chat_and_var_tmpa$lcl_chat[v]&
                                               chat_and_var_tmpa$ucl_chat[v]>chat_and_var_tmpa$Truth[v]
                                             ,1,0)  
}

coverages_chat_SC_N<-chat_and_var_tmpa%>%
  group_by(cohesion, aggregation, p0)%>%
  summarise(N_coverage_chat=sum(coverage_chat,na.rm = TRUE)/100)

chat_and_var_tmpa<-chat_and_var_tmpa[,which(colnames(chat_and_var_tmpa)%in%
                                             c("chat","p0","cohesion","aggregation","sim","N_var" ))]

#then do for sigma
tmp<-SCresults[SCresults$param=="sigma",]
colnames(tmp)

chat_and_var_tmpb<-left_join(chat_and_var,
                        tmp,
                        #tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","St.Dev."))],
                        by=c("sim","p0","cohesion","aggregation"))
colnames(chat_and_var_tmpb)[which(colnames(chat_and_var_tmpb)=="St.Dev.")]<-"sigma_var"
chat_and_var_tmpb$sigma_var<-chat_and_var_tmpb$sigma_var^2

#chat corrected variance etc.
chat_and_var_tmpb$sigma_var_chat <-chat_and_var_tmpb$sigma_var * chat_and_var_tmpb$chat
chat_and_var_tmpb$lcl_chat <-chat_and_var_tmpb$Mean-1.96*sqrt(chat_and_var_tmpb$sigma_var_chat)
chat_and_var_tmpb$ucl_chat <-chat_and_var_tmpb$Mean+1.96*sqrt(chat_and_var_tmpb$sigma_var_chat)
for(v in 1:nrow(chat_and_var_tmpb)){
  chat_and_var_tmpb$coverage_chat[v] <-ifelse(chat_and_var_tmpb$Truth[v]>chat_and_var_tmpb$lcl_chat[v]&
                                                chat_and_var_tmpb$ucl_chat[v]>chat_and_var_tmpb$Truth[v]
                                              ,1,0)  
}

coverages_chat_SC_sigma<-chat_and_var_tmpb%>%
  group_by(cohesion, aggregation, p0)%>%
  summarise(sigma_coverage_chat=sum(coverage_chat,na.rm = TRUE)/100)

chat_and_var_tmpb<-chat_and_var_tmpb[,which(colnames(chat_and_var_tmpb)%in%
                                              c("chat","p0","cohesion","aggregation","sim","sigma_var" ))]

chat_and_var_SC<-full_join(chat_and_var_tmpa,chat_and_var_tmpb)

#get mean across sims 
chat_and_var_SC_means<-chat_and_var_SC%>%
  group_by( p0, cohesion, aggregation)%>%
  summarize(chat=mean(chat),
            NVar=mean(N_var,na.rm=TRUE),SigmaVar=mean(sigma_var,na.rm=TRUE))

# add new coverage
chat_and_var_SC_means<-left_join(chat_and_var_SC_means,
                                  coverages_chat_SC_N)
chat_and_var_SC_means<-left_join(chat_and_var_SC_means,
                                  coverages_chat_SC_sigma)

#calculate RV
chat_and_var_SC_means$RV_N<-NA
chat_and_var_SC_means$RV_Sigma<-NA
chat_and_var_SC_means$RV_N[1:2]<-1
chat_and_var_SC_means$RV_Sigma[1:2]<-1

for(i in 3:nrow(chat_and_var_SC_means)){
  if(chat_and_var_SC_means$p0[i]=="0.05"){
    chat_and_var_SC_means$RV_N[i]<-chat_and_var_SC_means$NVar[i]/chat_and_var_SC_means$NVar[1]
    chat_and_var_SC_means$RV_Sigma[i]<-chat_and_var_SC_means$SigmaVar[i]/chat_and_var_SC_means$SigmaVar[1]
  }
  else{#for p0=0.20
    chat_and_var_SC_means$RV_N[i]<-chat_and_var_SC_means$NVar[i]/chat_and_var_SC_means$NVar[2]
    chat_and_var_SC_means$RV_Sigma[i]<-chat_and_var_SC_means$SigmaVar[i]/chat_and_var_SC_means$SigmaVar[2]
  }
}



### for SPIM #### 
PID <- unique(SPIMresults$PID)#collar, sex, sexcollar, etc. 

chat_and_var_SPIM_means<-data.frame()
for(p in 1:length(PID)){ #for each PID
 
  #first do for N
  tmp<-SPIMresults[SPIMresults$param=="N"&SPIMresults$PID==PID[p],]
  colnames(tmp)
  
  chat_and_var_tmpa<-left_join(chat_and_var,
                          tmp,
                          #tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","SD"))],
                          by=c("sim","p0","cohesion","aggregation"))
  colnames(chat_and_var_tmpa)[which(colnames(chat_and_var_tmpa)=="SD")]<-"N_var"
  chat_and_var_tmpa$N_var<-chat_and_var_tmpa$N_var^2
  
  #chat corrected variance etc.
  chat_and_var_tmpa$N_var_chat <-chat_and_var_tmpa$N_var * chat_and_var_tmpa$chat
  chat_and_var_tmpa$lcl_chat <-chat_and_var_tmpa$Mean-1.96*sqrt(chat_and_var_tmpa$N_var_chat)
  chat_and_var_tmpa$ucl_chat <-chat_and_var_tmpa$Mean+1.96*sqrt(chat_and_var_tmpa$N_var_chat)
  for(v in 1:nrow(chat_and_var_tmpa)){
    chat_and_var_tmpa$coverage_chat[v] <-ifelse(chat_and_var_tmpa$Truth[v]>chat_and_var_tmpa$lcl_chat[v]&
                                                  chat_and_var_tmpa$ucl_chat[v]>chat_and_var_tmpa$Truth[v]
                                                ,1,0)  
  }
  coverages_chat_SPIM_N_tmp<-chat_and_var_tmpa%>%
    group_by(cohesion, aggregation, p0)%>%
    summarise(N_coverage_chat=sum(coverage_chat,na.rm = TRUE)/100)
  
  chat_and_var_tmpa<-chat_and_var_tmpa[,which(colnames(chat_and_var_tmpa)%in%
                                                c("chat","p0","cohesion","aggregation","sim","N_var" ))]
  
  
  #then do for sigma
  tmp<-SPIMresults[SPIMresults$param=="sigma"&SPIMresults$PID==PID[p],]
  colnames(tmp)
  
  chat_and_var_tmpb<-left_join(chat_and_var,
                          tmp,
                          #tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","SD"))],
                          by=c("sim","p0","cohesion","aggregation"))
  colnames(chat_and_var_tmpb)[which(colnames(chat_and_var_tmpb)=="SD")]<-"sigma_var"
  chat_and_var_tmpb$sigma_var<-chat_and_var_tmpb$sigma_var^2
  
  #chat corrected variance etc.
  chat_and_var_tmpb$sigma_var_chat<-chat_and_var_tmpb$sigma_var * chat_and_var_tmpb$chat
  chat_and_var_tmpb$lcl_chat <-chat_and_var_tmpb$Mean-1.96*sqrt(chat_and_var_tmpb$sigma_var_chat)
  chat_and_var_tmpb$ucl_chat <-chat_and_var_tmpb$Mean+1.96*sqrt(chat_and_var_tmpb$sigma_var_chat)
  for(v in 1:nrow(chat_and_var_tmpb)){
    chat_and_var_tmpb$coverage_chat[v] <-ifelse(chat_and_var_tmpb$Truth[v]>chat_and_var_tmpb$lcl_chat[v]&
                                                  chat_and_var_tmpb$ucl_chat[v]>chat_and_var_tmpb$Truth[v]
                                                ,1,0)  
  }
  coverages_chat_SPIM_sigma_tmp<-chat_and_var_tmpb%>%
    group_by(cohesion, aggregation, p0)%>%
    summarise(sigma_coverage_chat=sum(coverage_chat,na.rm = TRUE)/100)
  
  chat_and_var_tmpb<-chat_and_var_tmpb[,which(colnames(chat_and_var_tmpb)%in%
                                                c("chat","p0","cohesion","aggregation","sim","sigma_var" ))]
  
  
  chat_and_var_SPIM_part<-full_join(chat_and_var_tmpa,chat_and_var_tmpb)
  
  #get mean across sims 
  chat_and_var_SPIM_part_means<-chat_and_var_SPIM_part%>%
    group_by(p0, cohesion, aggregation)%>% #scenario,
    summarize(chat=mean(chat),
              NVar=mean(N_var,na.rm=TRUE),SigmaVar=mean(sigma_var,na.rm=TRUE))
 
   # add new coverage
  chat_and_var_SPIM_part_means<-left_join(chat_and_var_SPIM_part_means,
                                          coverages_chat_SPIM_N_tmp )
  chat_and_var_SPIM_part_means<-left_join(chat_and_var_SPIM_part_means,
                                          coverages_chat_SPIM_sigma_tmp )
  #calculate RV
  chat_and_var_SPIM_part_means$RV_N<-NA
  chat_and_var_SPIM_part_means$RV_Sigma<-NA
  chat_and_var_SPIM_part_means$RV_N[1:2]<-1
  chat_and_var_SPIM_part_means$RV_Sigma[1:2]<-1
  
  for(i in 3:nrow(chat_and_var_SPIM_part_means)){
    if(chat_and_var_SPIM_part_means$p0[i]=="0.05"){
      chat_and_var_SPIM_part_means$RV_N[i]<-chat_and_var_SPIM_part_means$NVar[i]/chat_and_var_SPIM_part_means$NVar[1]
      chat_and_var_SPIM_part_means$RV_Sigma[i]<-chat_and_var_SPIM_part_means$SigmaVar[i]/chat_and_var_SPIM_part_means$SigmaVar[1]
    }
    else{#for p0=0.20
      chat_and_var_SPIM_part_means$RV_N[i]<-chat_and_var_SPIM_part_means$NVar[i]/chat_and_var_SPIM_part_means$NVar[2]
      chat_and_var_SPIM_part_means$RV_Sigma[i]<-chat_and_var_SPIM_part_means$SigmaVar[i]/chat_and_var_SPIM_part_means$SigmaVar[2]
    }
  }
  
  #add PID
  chat_and_var_SPIM_part_means$PID<-PID[p]
  
  chat_and_var_SPIM_means<-rbind(chat_and_var_SPIM_means,chat_and_var_SPIM_part_means)
}
  
table(chat_and_var_SPIM_means$PID)
chat_and_var_SPIM_means$Model<-paste0("SPIM: ",chat_and_var_SPIM_means$PID)


#### SCR, SC, and SPIM togehetr ####
colnames(chat_and_var_SCR_means)
colnames(chat_and_var_SC_means)
colnames(chat_and_var_SPIM_means)

chat_and_var_SCR_means$Model<-"SCR"
chat_and_var_SC_means$Model<-"SC"

chat_and_var_means<-rbind.fill(chat_and_var_SCR_means,
                         chat_and_var_SC_means,
                         chat_and_var_SPIM_means)
chat_and_var_means$Model<-factor(chat_and_var_means$Model,levels=c("SCR","SC",
                                                       "SPIM: collar",
                                                       "SPIM: sex",
                                                       "SPIM: sexcollar",
                                                       "SPIM: sexcoat",
                                                       "SPIM: sexcollarcoat",
                                                       "SPIM: antlerssexcollarcoat"))
colBluefunc <- colorRampPalette(c("lightblue","blue"))
colBluefunc(length(levels(chat_and_var_means$Model))-2)
cbbPalette <- c("white","pink",colBluefunc(length(levels(chat_and_var_means$Model))-2))# "#E69F00", "#56B4E9", "#009E73", "#F0E442")
cbbPalette2 <-cbbPalette
cbbPalette2[1]<-"black"


plot_chatRV_N_mean<-ggplot()+
  geom_abline(lty=2)+
  geom_point(data=chat_and_var_means,
             size=3,position=position_dodge(width=2),
             aes(x=chat,y=RV_N,color=Model,group=Model,shape=Model))+
  geom_line(data=chat_and_var_means,
            aes(x=chat, y=RV_N, group=Model,color=Model),
            position=position_dodge(width=2))+
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  #  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(cols=vars(cohesion),rows=vars(p0))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_chatRV_Sigma_mean<-ggplot()+
  geom_abline(lty=2)+
  geom_point(data=chat_and_var_means,
             size=3,position=position_dodge(width=1),
             aes(x=chat,y=RV_Sigma,color=Model,group=Model,shape=Model))+
  geom_line(data=chat_and_var_means,
            aes(x=chat, y=RV_Sigma, group=Model,color=Model),
            position=position_dodge(width=1))+
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  #  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(cols=vars(cohesion),rows=vars(p0))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))




### write CHat, RV, and corrected Coverages to file ####
colnames(chat_and_var_means)
write.csv(chat_and_var_means,"summaryChat_RV_coverage_long.csv")
summaryChat_wide<-chat_and_var_means[,colnames(chat_and_var_means)%in%c("p0","cohesion","aggregation","chat","Model")]
summaryChat_wide<-pivot_wider(summaryChat_wide,names_from=cohesion,values_from = chat)
write.csv(summaryChat_wide,"summaryChat_wide.csv")

summaryRV_N_wide<-chat_and_var_means[,colnames(chat_and_var_means)%in%c("p0","cohesion","aggregation","RV_N","Model")]
summaryRV_N_wide<-pivot_wider(summaryRV_N_wide,names_from=cohesion,values_from = RV_N)
write.csv(summaryRV_N_wide,"summaryRV_N_wide.csv")

summaryRV_Sigma_wide<-chat_and_var_means[,colnames(chat_and_var_means)%in%c("p0","cohesion","aggregation","RV_Sigma","Model")]
summaryRV_Sigma_wide<-pivot_wider(summaryRV_Sigma_wide,names_from=cohesion,values_from = RV_Sigma)
write.csv(summaryRV_N_wide,"summaryRV_Sigma_wide.csv")


#### new coverage plots corrected with chat #####

plot_N_coverage_chat_SCR<-ggplot() + # 
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_line(data=coverages_chat_SCR_N,
            aes(x= as.numeric(as.character(coverages_chat_SCR_N$aggregation)), 
                y=N_coverage_chat, group=cohesion))+#,color=cohesion))+
  geom_point(data=coverages_chat_SCR_N, size=3,
             aes(x= as.numeric(as.character(coverages_chat_SCR_N$aggregation)),
                 y=N_coverage_chat,group=cohesion), pch=21,fill="white")+#shape=cohesion,color=cohesion))+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  ylim(c(0,1))+
  labs(title="Abundance (N),chat corrected", x="Aggregation (Group Size)",
       y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sigma_coverage_chat_SCR<-ggplot() + # 
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_line(data=coverages_chat_SCR_sigma,
            aes(x= as.numeric(as.character(coverages_chat_SCR_sigma$aggregation)), 
                y=sigma_coverage_chat,group=cohesion))+#,color=cohesion))+
  geom_point(data=coverages_chat_SCR_sigma, size=3,
             aes(x= as.numeric(as.character(coverages_chat_SCR_sigma$aggregation)), 
                 y=sigma_coverage_chat,group=cohesion),pch=21,fill="white")+#shape=cohesion,color=cohesion))+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  ylim(c(0,1))+
  labs(title="Sigma (\u03c3),chat corrected", x="Aggregation (Group Size)",
       y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_N_coverage_chat_SC<-ggplot() + # 
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_line(data=coverages_chat_SC_N,
            aes(x= as.numeric(as.character(coverages_chat_SC_N$aggregation)), 
                y=N_coverage_chat, group=cohesion))+#,color=cohesion))+
  geom_point(data=coverages_chat_SC_N, size=3,
             aes(x= as.numeric(as.character(coverages_chat_SC_N$aggregation)),
                 y=N_coverage_chat, group=cohesion),pch=21,fill="white")+#shape=cohesion,color=cohesion))+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  ylim(c(0,1))+
  labs(title="Abundance (N),chat corrected", x="Aggregation (Group Size)",
       y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sigma_coverage_chat_SC<-ggplot() + # 
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_line(data=coverages_chat_SC_sigma,
            aes(x= as.numeric(as.character(coverages_chat_SC_sigma$aggregation)), 
                y=sigma_coverage_chat, group=cohesion))+#,color=cohesion))+
  geom_point(data=coverages_chat_SC_sigma, size=3,
             aes(x= as.numeric(as.character(coverages_chat_SC_sigma$aggregation)),
                 y=sigma_coverage_chat,group=cohesion),pch=21,fill="white")+#shape=cohesion,color=cohesion))+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  ylim(c(0,1))+
  labs(title="Sigma (\u03c3),chat c0rrected", x="Aggregation (Group Size)",
       y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#create a temp df of the corrected spim coverages to exclude 
#p0.20 sexcoat at 4 ag and 0 cohesion, which failed
tmp<-chat_and_var_SPIM_means[-which(is.na(chat_and_var_SPIM_means$NVar)),]
plot_N_coverage_chat_SPIM<-ggplot()+
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=tmp, 
             size=3,position=position_dodge(width=2),
             aes(as.numeric(as.character(tmp$aggregation)), 
                 y=N_coverage_chat, shape=PID,group=PID,color=PID))+ #group=cohesion
  geom_line(data=tmp, 
            aes(as.numeric(as.character(tmp$aggregation)), 
                y=N_coverage_chat, group=PID,color=PID),
            position=position_dodge(width=2))+
  scale_color_brewer(palette="Blues")+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs,cohesion=coh.labs),
             scales="free")+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sigma_coverage_chat_SPIM<-ggplot() +
  geom_point(data=tmp,
             size=3,position=position_dodge(width=2),
             aes(as.numeric(as.character(tmp$aggregation)), 
                 y=sigma_coverage_chat,
                 shape=PID,group=PID,color=PID))+
  geom_line(data=tmp,
            position=position_dodge(width=2),
            aes(as.numeric(as.character(tmp$aggregation)), 
                y=sigma_coverage_chat, group=PID,color=PID))+ 
   scale_color_brewer(palette="Blues")+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion),
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs),
             scales="free")+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


# all together now
tmp<-chat_and_var_means[-which(is.na(chat_and_var_means$NVar)),]

plot_N_coverage_chat<-ggplot() +
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_hline(yintercept = 0.90,lty=2,)+
  geom_point(data=tmp,
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(as.numeric(as.character(tmp$aggregation)), 
                 y=N_coverage_chat,shape=Model,group=Model,color=Model))+
  geom_line(data=tmp,
            position=position_dodge(width=2),
            aes(as.numeric(as.character(tmp$aggregation)), 
                y=N_coverage_chat, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), scales = "free",
             labeller=labeller(p0 = p0.labs,cohesion=coh.labs))+
  labs(title="Abundance (N)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_sigma_coverage_chat_SC<-ggplot() +
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=tmp,
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(as.numeric(as.character(tmp$aggregation)), 
                 y=sigma_coverage_chat,shape=Model,group=Model,color=Model))+
  geom_line(data=tmp,
            position=position_dodge(width=2),
            aes(as.numeric(as.character(tmp$aggregation)), 
                y=sigma_coverage_chat, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(p0),cols=vars(cohesion), scales = "free",
             labeller=labeller(p0 = p0.labs, cohesion=coh.labs))+
  labs(title="Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#rearrange from wide to long
colnames(tmp)
chat_and_var_means_long_a<-tmp[,c(1:4,11,12)]
chat_and_var_means_long_a$param<-"N"
chat_and_var_means_long_a<-cbind(chat_and_var_means_long_a,tmp[,which(grepl("N",colnames(tmp)))])

chat_and_var_means_long_b<-tmp[,c(1:4,11,12)]
chat_and_var_means_long_b$param<-"Sigma"
chat_and_var_means_long_b<-cbind(chat_and_var_means_long_b,tmp[,which(grepl("igma",colnames(tmp)))])

colnames(chat_and_var_means_long_b)<-c("p0","cohesion","aggregation",
                                       "chat","Model","PID","param",
                                       "Var","Coverage_chat","RV")
colnames(chat_and_var_means_long_a)<-colnames(chat_and_var_means_long_b)
chat_and_var_means_long<-rbind(chat_and_var_means_long_a,
                               chat_and_var_means_long_b)
plot_coverage_chat_p05<-ggplot() + 
  geom_hline(yintercept = 0.95,lty=2,)+
  geom_point(data=chat_and_var_means_long[chat_and_var_means_long$p0=="0.05",], #0.2
             size=3,position=position_dodge(width=2),alpha=0.5,
             aes(as.numeric(as.character(aggregation)), 
                 y=Coverage_chat, shape=Model,group=Model,color=Model))+
  geom_line(data=chat_and_var_means_long[chat_and_var_means_long$p0=="0.05",],
            position=position_dodge(width=2),
            aes(as.numeric(as.character(aggregation)), 
                y=Coverage_chat, group=Model,color=Model))+ 
  scale_shape_manual(values=c(1,17,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(rows=vars(param),cols=vars(cohesion),switch="y", scales = "free",
             labeller=labeller(cohesion=coh.labs))+
  labs(x="Aggregation (Group Size)",y="Coverage")+
  ylim(c(0,1))+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

plot_coverage_chat_p05
