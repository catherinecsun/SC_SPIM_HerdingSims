#Run SCR model on data

#libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(secr)

##### The data ####
#scenario 1:8 all have aggregation of 1 so the cohesion is irrelevant
#what changes is the detection probability
# so 1, 3, 5, and 7 are identical (p0 of 0.05 but coohesoon of 0,0.3, 0.67 and 1)
# and2 4, 6, and 8 are the same
#so can remove 3:8

sim_data<-readRDS("../Data/simulatedData_from24Scenarios.rds")

View(sim_data)
sim_data<-sim_data[-c(3:8)]
parm_combos<-parm_combos[-c(3:8),]

##### how many batches of sims to do/cores to use? ####
batches<-100
#the buffer is 3*sigma (9)
# and the mask area is  0.1729688 ha 
# as is reported in summary of scr.fit
# and the reported is  Fitted (ie the Real D)  per ha
# so the N is the reported/0.1729688

for(w in 1:length(sim_data)){ #for every parmcomb
  summary_SCR<-data.frame()
  for(b in 1:batches){ 
    tmp_capthist<-sim_data[[scen]][[b]]$capthist
    tmp_fit<-secr.fit(capthist=tmp_capthist,
               buffer=3*sigma,model=list(D~1,g0~1,sigma~1)) #since sigma is 3km or 3000m
    tmp_pred<-predict(tmp_fit)
    tmp_pred[4,]<-c("",tmp_pred[which(rownames(tmp_pred)=="D"),2:5]*0.1729688)
    rownames(tmp_pred)[4]<-"N"                
    tmp_pred<-cbind(scenario=scen,param=rownames(tmp_pred),sim=b,tmp_pred)
    summary_SCR<-rbind(summary_SCR,tmp_pred)
  }
  write.csv(summary_SCR,paste0("../SCRresults/SCR_simResults_scenario",w,
                              "_N",parm_combos$N.inds[w],
                              "_p0",parm_combos$p0[w],
                              "_coh",parm_combos$cohesion[w],
                              "_agg",parm_combos$aggregation[w],
                              ".csv"))
  
}
esa.plot(tmp_fit)
