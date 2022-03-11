
### relative variance ####
# the ratio of the empirical variance amount simulations for any give scenario
# and the variance associated with the independent scenario


### for SC ####
#create a dataframe to hold chat and rv per sim of each scenario
chat_and_var<-chat_df_long
chat_and_var$sim<-rep(1:100,length(unique(chat_and_var$scenario)))
colnames(chat_and_var)

  #first do for N
tmp<-SCresults[SCresults$param=="N",]
colnames(tmp)

chat_and_var<-left_join(chat_and_var,
                       tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","St.Dev."))],
                       by=c("sim","p0","cohesion","aggregation"))
colnames(chat_and_var)[which(colnames(chat_and_var)=="St.Dev.")]<-"N_var"
chat_and_var$N_var<-chat_and_var$N_var^2

  #then do for sigma
tmp<-SCresults[SCresults$param=="sigma",]
colnames(tmp)

chat_and_var<-left_join(chat_and_var,
                       tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","St.Dev."))],
                       by=c("sim","p0","cohesion","aggregation"))
colnames(chat_and_var)[which(colnames(chat_and_var)=="St.Dev.")]<-"sigma_var"
chat_and_var$sigma_var<-chat_and_var$sigma_var^2

chat_and_var_SC<-chat_and_var

#get mean across sims 
chat_and_var_SC<-chat_and_var_SC%>%
  group_by(scenario, p0, cohesion, aggregation)%>%
  summarize(chat=mean(chat),
            NVar=mean(N_var,na.rm=TRUE),SigmaVar=mean(sigma_var,na.rm=TRUE))

#remove scenarios 3:8 that are equivalent to scenario 1
chat_and_var_SC<-chat_and_var_SC[-which(is.na(chat_and_var_SC$NVar)),]

#calculate RV
chat_and_var_SC$RV_N<-NA
chat_and_var_SC$RV_Sigma<-NA
chat_and_var_SC$RV_N[1:2]<-1
chat_and_var_SC$RV_Sigma[1:2]<-1

for(i in 3:nrow(chat_and_var_SC)){
  if(chat_and_var_SC$p0[i]=="0.05"){
    chat_and_var_SC$RV_N[i]<-chat_and_var_SC$NVar[i]/chat_and_var_SC$NVar[1]
    chat_and_var_SC$RV_Sigma[i]<-chat_and_var_SC$SigmaVar[i]/chat_and_var_SC$SigmaVar[1]
  }
  else{#for p0=0.20
    chat_and_var_SC$RV_N[i]<-chat_and_var_SC$NVar[i]/chat_and_var_SC$NVar[2]
    chat_and_var_SC$RV_Sigma[i]<-chat_and_var_SC$SigmaVar[i]/chat_and_var_SC$SigmaVar[2]
  }
}

plot_chat_SC_mean<-ggplot()+
  geom_point(data=chat_and_var_SC,aes(x=aggregation,y=chat,color=cohesion))+
  geom_line(data=chat_and_var_SC,aes(x=aggregation,y=chat,color=cohesion,group=cohesion))+
  facet_wrap(vars(p0))

plot_chatRV_SC_N<-ggplot()+
  geom_point(data=chat_and_var_SC,aes(x=chat,y=RV_N,color=cohesion))+
  geom_line(data=chat_and_var_SC,aes(x=chat,y=RV_N,color=cohesion,group=cohesion))+
  facet_wrap(vars(p0))

plot_chatRV_SC_Sigma<-ggplot()+
  geom_point(data=chat_and_var_SC,aes(x=chat,y=RV_Sigma,color=cohesion))+
  geom_line(data=chat_and_var_SC,aes(x=chat,y=RV_Sigma,color=cohesion,group=cohesion))+
  facet_wrap(vars(p0))


### for SPIM #### 
PID #collar, sex, sexcollar, etc. 
chat_and_var_SPIM<-data.frame()
for(p in 1:length(PID)){ #for each PID
  #create a dataframe to hold chat and rv per sim of each scenario
  chat_and_var<-chat_df_long
  chat_and_var$sim<-rep(1:100,length(unique(chat_and_var$scenario)))
  colnames(chat_and_var)
  
  #first do for N
  tmp<-SPIMresults[SPIMresults$param=="N"&SPIMresults$PID==PID[p],]
  colnames(tmp)
  
  chat_and_var<-left_join(chat_and_var,
                          tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","SD"))],
                          by=c("sim","p0","cohesion","aggregation"))
  colnames(chat_and_var)[which(colnames(chat_and_var)=="SD")]<-"N_var"
  chat_and_var$N_var<-chat_and_var$N_var^2
  
  #then do for sigma
  tmp<-SPIMresults[SPIMresults$param=="sigma"&SPIMresults$PID=="antlerssexcollarcoat",]
  colnames(tmp)
  
  chat_and_var<-left_join(chat_and_var,
                          tmp[,which(colnames(tmp)%in%c("sim","p0","cohesion","aggregation","SD"))],
                          by=c("sim","p0","cohesion","aggregation"))
  colnames(chat_and_var)[which(colnames(chat_and_var)=="SD")]<-"sigma_var"
  chat_and_var$sigma_var<-chat_and_var$sigma_var^2
  
  chat_and_var_SPIM_part<-chat_and_var
  
  #get mean across sims 
  chat_and_var_SPIM_part<-chat_and_var_SPIM_part%>%
    group_by(scenario, p0, cohesion, aggregation)%>%
    summarize(chat=mean(chat),
              NVar=mean(N_var,na.rm=TRUE),SigmaVar=mean(sigma_var,na.rm=TRUE))
  
  #remove scenarios 3:8 that are equivalent to scenario 1
  chat_and_var_SPIM_part<-chat_and_var_SPIM_part[-which(is.na(chat_and_var_SPIM_part$NVar)),]
  
  #calculate RV
  chat_and_var_SPIM_part$RV_N<-NA
  chat_and_var_SPIM_part$RV_Sigma<-NA
  chat_and_var_SPIM_part$RV_N[1:2]<-1
  chat_and_var_SPIM_part$RV_Sigma[1:2]<-1
  
  for(i in 3:nrow(chat_and_var_SPIM_part)){
    if(chat_and_var_SPIM_part$p0[i]=="0.05"){
      chat_and_var_SPIM_part$RV_N[i]<-chat_and_var_SPIM_part$NVar[i]/chat_and_var_SPIM_part$NVar[1]
      chat_and_var_SPIM_part$RV_Sigma[i]<-chat_and_var_SPIM_part$SigmaVar[i]/chat_and_var_SPIM_part$SigmaVar[1]
    }
    else{#for p0=0.20
      chat_and_var_SPIM_part$RV_N[i]<-chat_and_var_SPIM_part$NVar[i]/chat_and_var_SPIM_part$NVar[2]
      chat_and_var_SPIM_part$RV_Sigma[i]<-chat_and_var_SPIM_part$SigmaVar[i]/chat_and_var_SPIM_part$SigmaVar[2]
    }
  }
  chat_and_var_SPIM_part$PID<-PID[p]
  chat_and_var_SPIM<-rbind(chat_and_var_SPIM,chat_and_var_SPIM_part)

}
  
table(chat_and_var_SPIM$PID)

chat_and_var_SPIM$Model<-paste0("SPIM: ",chat_and_var_SPIM$PID)

plot_chat_SPIM_mean<-ggplot()+
  geom_point(data=chat_and_var_SPIM,
             size=3,position=position_dodge(width=2),
             aes(x=aggregation,y=chat,color=PID,group=PID,shape=PID))+
  geom_line(data=chat_and_var_SPIM,
            aes(x=aggregation, y=chat, group=PID,color=PID),
            position=position_dodge(width=2))+
  scale_color_brewer(palette="Blues")+
  facet_wrap(vars(p0,cohesion))

plot_chatRV_SC_N<-ggplot()+
   geom_point(data=chat_and_var_SPIM,
             size=3,position=position_dodge(width=2),
             aes(x=chat,y=RV_N,color=PID,group=PID,shape=PID))+
  geom_line(data=chat_and_var_SPIM,
            aes(x=chat,y=RV_N, group=PID,color=PID),
            position=position_dodge(width=2))+
  scale_color_brewer(palette="Blues")+
  facet_wrap(vars(p0,cohesion))

plot_chatRV_SC_Sigma<-ggplot()+
  geom_point(data=chat_and_var_SPIM,
           size=3,position=position_dodge(width=2),
           aes(x=chat,y=RV_Sigma,color=PID,group=PID,shape=PID))+
  geom_line(data=chat_and_var_SPIM,
            aes(x=chat,y=RV_Sigma, group=PID,color=PID),
            position=position_dodge(width=2))+
  scale_color_brewer(palette="Blues")+
  facet_wrap(vars(p0,cohesion))

#### SC ad SPIM togehetr ####
colnames(chat_and_var_SC)
colnames(chat_and_var_SPIM)

chat_and_var_SC$Model<-"SC"
chat_and_var_SPIM

chat_and_var<-rbind.fill(chat_and_var_SC,chat_and_var_SPIM)
chat_and_var$Model<-factor(chat_and_var$Model,levels=c("SC",
                                                       "SPIM: collar",
                                                       "SPIM: sex",
                                                       "SPIM: sexcollar",
                                                       "SPIM: sexcoat",
                                                       "SPIM: sexcollarcoat",
                                                       "SPIM: antlerssexcollarcoat"))
colBluefunc <- colorRampPalette(c("lightblue","blue"))
colBluefunc(length(levels(chat_and_var$Model))-1)
cbbPalette <- c("white",colBluefunc(length(levels(chat_and_var$Model))-1))# "#E69F00", "#56B4E9", "#009E73", "#F0E442")
cbbPalette2 <-cbbPalette
cbbPalette2[1]<-"black"

plot_chat_mean<-ggplot()+
  geom_point(data=chat_and_var,
             size=3,position=position_dodge(width=1),
             aes(x=aggregation,y=chat,color=Model,group=Model,shape=Model))+
  geom_line(data=chat_and_var,
            aes(x=aggregation, y=chat, group=Model,color=Model),
            position=position_dodge(width=1))+
  scale_shape_manual(values=c(1,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
#  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(cols=vars(cohesion),rows=vars(p0))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))



plot_chatRV_N_mean<-ggplot()+
  geom_point(data=chat_and_var,
             size=3,position=position_dodge(width=1),
             aes(x=chat,y=RV_N,color=Model,group=Model,shape=Model))+
  geom_line(data=chat_and_var,
            aes(x=chat, y=RV_N, group=Model,color=Model),
            position=position_dodge(width=1))+
  scale_shape_manual(values=c(1,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  #  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(cols=vars(cohesion),rows=vars(p0))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))



plot_chatRV_Sigma_mean<-ggplot()+
  geom_point(data=chat_and_var,
             size=3,position=position_dodge(width=1),
             aes(x=chat,y=RV_Sigma,color=Model,group=Model,shape=Model))+
  geom_line(data=chat_and_var,
            aes(x=chat, y=RV_Sigma, group=Model,color=Model),
            position=position_dodge(width=1))+
  scale_shape_manual(values=c(1,17,17,17,17,17,17))+#,17,15,3,7,8,1))+
  scale_color_manual(values=cbbPalette2)+
  #  scale_x_continuous(breaks=c(1,4,10),labels=c(1,4,10))+
  facet_grid(cols=vars(cohesion),rows=vars(p0))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))