#Fletcher et al. 2012 overdispersion
#overdispersion factor >1 signifies overdispersed data, 
#i.e. when the variance in the data exceeds that predicted by the statistical model

#this is of the individual level capture histories
#before the data are simplified into SC or SPIM datasets for model/estimation.
#ie pattern of overdispersion in the actual population


#packages
#cant seem to most load packages so everything will just be base
library(ggplot2)
library(tidyr)
library(dplyr)

#remove duplicated scenarios
#parm_combos<-parm_combos[-c(3:8),]
#sim_data<-sim_data[-c(3:8)]


### chat  ####
#create an empty dataframe to put the calculated chats 
chat_df<-data.frame(row.names=c(1:nsims))

#for each scenario, and each sim in that scenario
for(i in 1:length(sim_data)){
  print(i)
  #create an empty vector to put the calculated chats for this scenario
  chat_vec<-c()
  
  for(j in 1:nsims){
    #get the simulated 
    tmp_data<-sim_data[[i]][[j]]
    
    #expected number of detections
    area<-(range(tmp_data$mask[,1])[2]-range(tmp_data$mask[,1])[1])*(range(tmp_data$mask[,2])[2]-range(tmp_data$mask[,2])[1])
    D <- tmp_data$parms$N.inds/area #density
    fn <- function (r) r * tmp_data$detectpar$g0 * exp(-r^2/2/(tmp_data$detectpar$sigma/100)^2)
    expected.nk <- 2 * pi *D * integrate(fn, 0, 5*tmp_data$detectpar$sigma/100)$value
    
    #realized number of detections 
    #check if tmp_data$capthist is truly ch2!
    nk <- apply(apply(tmp_data$capthist,c(1,3),sum)>0,2,sum)  
    
    X2 <- sum((nk - expected.nk)^2 / expected.nk)
    si <- mean((nk - expected.nk) / expected.nk) 
    
    # number of traps - number of parameters in model
    # number of parameters in model is 1 if we assume homogeneous density
    nu <- dim(tmp_data$capthist)[3] - 1
    
    chat <- X2/nu / (1 + si) 
    chat_vec<-c(chat_vec,chat)
  } # for each sim
  
  chat_df<-cbind(chat_df,chat_vec)
  colnames(chat_df)[i]<-i
}#for each scenario


#turn chat_df from wide to long
chat_df_long<-c()
for(c in 1:ncol(chat_df)){
  chat_df_long<-c(chat_df_long,chat_df[,c])
}
chat_df_long<-as.data.frame(chat_df_long)
colnames(chat_df_long)<-"chat"

# add columns for the parameter/scenario conditions
chat_df_long$scenario<-rep(1:length(sim_data),each=nsims)
chat_df_long$p0<-NA
chat_df_long$cohesion<-NA
chat_df_long$aggregation<-NA
for(r in 1:nrow(parm_combos)){
  chat_df_long$p0[chat_df_long$scenario==r]<-parm_combos$p0[r]
  chat_df_long$cohesion[chat_df_long$scenario==r]<-parm_combos$cohesion[r]
  chat_df_long$aggregation[chat_df_long$scenario==r]<-parm_combos$aggregation[r]
}

chat_df_long$cohesion<-as.factor(chat_df_long$cohesion)
chat_df_long$p0<-as.factor(chat_df_long$p0)
chat_df_long$aggregation<-as.factor(chat_df_long$aggregation)

chat_df_long%>%group_by(scenario)%>%summarize(mean=mean(chat))
str(chat_df_long)



chat_StatsSum<-chat_df_long%>%group_by(scenario,p0,aggregation,cohesion)%>%
  summarize(#p0,aggregation, cohesion,
    min(chat),
    mean(chat),
    max(chat),
    sd(chat))

#make chat_StatsSum into pretty tables
chat_StatsSum_Wide<-chat_StatsSum[,which(colnames(chat_StatsSum)%in% c("mean(chat)","sd(chat)","p0","cohesion","aggregation"))]
chat_StatsSum_Wide$chat<-paste0(round(chat_StatsSum_Wide$`mean(chat)`,1)," (",round(chat_StatsSum_Wide$`sd(chat)`,1),")")
chat_StatsSum_Wide<-chat_StatsSum_Wide[,-c(4,5)]
chat_StatsSum_Wide<-pivot_wider(chat_StatsSum_Wide,names_from = cohesion,values_from = chat)
chat_StatsSum_Wide<-chat_StatsSum_Wide[order(chat_StatsSum_Wide$p0),]



#plot
ggplot(data=chat_df_long,aes(x=cohesion, y=chat))+
  geom_boxplot()+
  ylim(c(0,20))

ggplot(data=chat_df_long,aes(x=aggregation, y=chat))+
  geom_boxplot()+
  ylim(c(0,20))

p0.labs <- c("p0: 0.05", "p0: 0.20")
names(p0.labs) <- c("0.05", "0.2")

coh.labs <- c("Cohesion: 0","Cohesion: 0.3", "Cohesion: 0.67", "Cohesion: 1")
names(coh.labs) <- c("0","0.3", "0.67", "1")

plot_chat<-ggplot(data=chat_df_long,aes(x=aggregation, y=chat))+
  geom_hline(yintercept=1,lty=2)+
  geom_boxplot(aes(fill=as.factor(p0)))+ #notch=TRUE
  facet_grid(cols=vars(cohesion), #rows=vars(p0),
             labeller = labeller(cohesion=coh.labs))+ #switch = "y" p0 = p0.labs,
  labs(title="Overdispersion ",
       x="Aggregation (Group Size)",y="Fletcher's C-hat",fill="Detection \nProbability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
plot_chat   



### relative variance ####
# the ratio of the empirical variance amount simulations for any give scenario
# and the variance associated with the independent scenario
