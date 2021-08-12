#Run SPIM model on data

#libraries
library(doParallel)
library(devtools)

Sys.setenv(TZ='UTC')
#install_github("benaug/SPIM")
library(SPIM)
library(coda)
library(tidyr)
library(ggplot2)


##### The data ####
#scenario 1:8 all have aggregation of 1 so the cohesion is irrelevant
#what changes is the detection probability
# so 1, 3, 5, and 7 are identical (p0 of 0.05 but coohesoon of 0,0.3, 0.67 and 1)
# and2 4, 6, and 8 are the same
#so can remove 3:8

View(sim_data)
#make sure this is the most recent version 
sim_data[[1]][[1]]$parms_IDcovs$sex_bin # should be [0.5, 0.5]
sim_data[[1]][[1]]$pop_ids$antlers # shouldnt be 999
sum(sim_data[[1]][[1]]$parms_IDcovs$antlers_cont) # should be 1
#sim_data<-sim_data[-c(3:8)]
#parm_combos<-parm_combos[-c(3:8),]

##### which scenarios? ####
whichScenario<-c(1)#c(14:18)

whichIDS<-c("antlers","sex", "collar" ,"coat") #out ofantlers sex collar coat


##### how many batches of sims to do/cores to use? ####
batches<-10


##### how many sims per batches to do? ####
#actually i have hardwired this below; not ideal but it gets the job done for now
#sims<-2
#batches*sims #total number of sims to do

##### how many iterations to run and burn per simulation? ####
niters<-  16666*3 #30000 650#
burn<- 4998#5000 50#  

#data<-sim_data[[9]][[1]]

##### the function(s) to do the modeling ####
#the model needs to compile inside the function because it uses the data
#not ideal because having to compile for every dataset takes a while

myFun2<-function(data,whichIDS,niters,burn,j){ #
  
  #assemble the data 
  
  #the static things
  M<-600 
  K<-data$parms$noccasions
  X<-data.frame(x=data$traps$x,y=data$traps$y)
  buff<-data$parms$habitat.buffer
  obstype<-'bernoulli'
  
  #the particular data
  data$capthist<-aperm(data$capthist, c(1,3,2)) # original but now: n x j x k
  
  #but detections across traps and occasions cant be associated with the same ind because we dont know ind id
  #so need to expand y.obs in the n dimension so that a [1,1] capture history becomes two histories [1,0] and [0,1]
  dets<- which(data$capthist>0,arr.ind = T)
  
  y.obs<-array(0,c(nrow(dets),dim(data$capthist)[2],dim(data$capthist)[3]))#NA,c(nrow(data$capthist)*dim(y.obs)[3],dim(data$capthist)[2],dim(data$capthist)[3]))
  
  for(d in 1:nrow(dets)){
    tmp_whichTrap<-dets[d,2]
    tmp_whichOcc<-dets[d,3]
    y.obs[d,tmp_whichTrap,tmp_whichOcc]<-1
  }
  
  #fix rownames of dets so that that the number after the (.) doesnt have a leading 0
  rownames(dets)<-paste0( sapply(strsplit(rownames(dets),"[.]"),"[[",1),".",
                         as.numeric(sapply(strsplit(rownames(dets),"[.]"),"[[",2)))
  
  data$pop_ids$Group_ID<-sapply(strsplit(data$pop_ids$Group_ID,"[.]"),"[[",1)
  data$pop_ids$ID<-(1:data$parms$group.size)
  data$pop_ids$ID<-paste0(data$pop_ids$Group_ID,".",data$pop_ids$ID)
  
  #partial IDs
  G.obs<-data$pop_ids[match(rownames(dets),data$pop_ids$ID),which(colnames(data$pop_ids)%in%whichIDS)] #(of the detected inds)
  
  #possible values and their proportions (evenly distributed across all values)
  IDlist<-list()
  IDlist$ncat<-ncol(G.obs)
  IDlist$IDcovs<-list()
  gammaVals<-list()
  
  for(c in 1:ncol(G.obs)){
    IDlist$IDcovs[[c]]<- 1:max(G.obs[,c])
    gammaVals[[c]]<-rep(1/max(G.obs[,c]),max(G.obs[,c])) #as.vector(table(G.obs[,c])/nrow(G.obs))
  }
  names(IDlist$IDcovs)<-colnames(G.obs)
  names(gammaVals)<-colnames(G.obs)
  
  #if collar is a category but there's only 1 value (ie actually no collared inds)
  # remove it from the list of partial covs
  # ie as if we dont know anything about collars aside from what we get on the cameras
  if("collar" %in% whichIDS){
    if(length(unique(G.obs$collar))==1){
      G.obs<-G.obs[,-which(colnames(G.obs)=="collar")]
      IDlist$ncat<- IDlist$ncat-1
      IDlist$IDcovs<-IDlist$IDcovs[-which(names(IDlist$IDcovs)=="collar")]
      gammaVals<-gammaVals[-which(names(gammaVals)=="collar")]
    }
  }

    
  dataSPIM<-list(y.obs=y.obs,G.obs=as.matrix(G.obs),IDlist=IDlist,X=X,K=K,buff=buff,obstype=obstype)
 
 
  nswap=nrow(y.obs)/2 
  inits <- list(psi=0.5,lam0=runif(1),sigma=rnorm(1,3),gamma=gammaVals)
  priors=list(sigma=c(24,8))
  
  #MCMC proposal parameters. Tuning these is the most difficult part.
  #shoot for rejection rates between 0.2 and 0.4; acceptance between 0.6-0.8 
  #If a parameter is not moving, the proposal parameter is too large.
  proppars=list(lam0=0.005,sigma=0.09,sx=1,sy=1) #for 1:16 with all partial ID covs, 
  #0.002 resulted in too high of an acceptance for lam0; 0.02 a touch too low.
  proppars=list(lam0=0.008,sigma=0.06,sx=1,sy=1) # for 17:24 with all partial ID covs
  IDup="MH" #Gibbs or metropolis-hastings latent ID update. Must use MH for binomial model.
  #Both about equally as fast for Poisson
  
  keepACs=TRUE #Do we store the activity center posteriors and other latent structure including the ID vector?
  keepGamma=TRUE #Do we store the gamma posterior?
  
  time_a=Sys.time()
  out<-mcmc.CatSPIM(dataSPIM,
                    niter=niters,nburn=0,nthin=1, nswap=nswap,
                    M = M, inits=inits,proppars=proppars,obstype=obstype,priors=priors,
                    IDup=IDup,keepACs=keepACs,keepGamma=keepGamma)
  time_b=Sys.time()

  
  #This will get you  acceptance probabilies. Can't change N, n, or psi.
 # 1-rejectionRate(mcmc(out$out))
  out$acceptRate<-1-rejectionRate(mcmc(out$out))
  
  #burn 
  out$out<-out$out[-c(1:burn),]
  
  
  #and split into 3 chains to check for convergence and effective size
  chain1<-as.mcmc(out$out[1:((niters-burn)/3),])
  chain2<-as.mcmc(out$out[(1+(niters-burn)/3):((2*(niters-burn)/3)),])
  chain3<-as.mcmc(out$out[(1+(2*(niters-burn)/3)):nrow(out$out),])
  chains<-as.mcmc.list(list(chain1,chain2,chain3))
  
  out$gel<-gelman.diag(chains,autoburnin = FALSE,multivariate=FALSE)
  
  #effective sample size;should shoot for at least 400 for N. 
  out$eff<-effectiveSize(chains)
  
  #clean up the summaries
  out$summary<-as.data.frame(cbind(summary(chains)$statistics,summary(chains)$quantiles))
  out$summary$sim<-j
  
  gammsProps<-matrix(ncol=8)
  for(g in 1:length(out$gammaOut)){ #for the partial id proportions
    for(v in 1:ncol(out$gammaOut[[g]])){
      gammsProps<-rbind(gammsProps,
                        cbind(paste0(colnames(G.obs)[g],v),
                              round(t(summary(out$gammaOut[[g]][,v])),3),
                              round(sd(summary(out$gammaOut[[g]][,v])),3)))
    }
  }
  gammsProps<-gammsProps[-1,]
  colnames(gammsProps)[ncol(gammsProps)]<-"SD"
  gammsProps<-as.data.frame(gammsProps)
  gammsProps$sim<-j
  out$gammaProps<-gammsProps
  
  
  #add the elapsed time, and the iterations
  out$time<-difftime(time_b,time_a,units="mins")
  out$itsTot<-niters-burn
  
  
  PID<-matrix(0,nrow=ncol(out$IDout),ncol=1+ncol(out$IDout))
  PID[,1]<-seq(1,ncol(out$IDout),1)
  PID<-as.data.frame(PID)
  colnames(PID)<-c("ID",seq(1,ncol(out$IDout),1))
  rownames(PID)<-colnames(PID)[-1]
  for(z in 1:ncol(out$IDout)){
    #print(z)
    check=z #Which sample to check
    storematch=matrix(0,nrow=ncol(out$IDout),ncol=niters)
    for(i in 1:niters){
      storematch[,i]=out$IDout[i,check]==out$IDout[i,]
    }
    PID[z,2:ncol(PID)]<-rowSums(storematch)/niters
  }
  PID<-PID%>%
    pivot_longer(-ID,names_to="ComparedTo_Ind",values_to="Probability")
  PID$ComparedTo_Ind<-as.numeric(PID$ComparedTo_Ind)
  
  out$PID<-PID
  
  #keep only the gel, summary, elapsed times, and its
  out2<-out[7:13]
  
  return(out2)
}

##### the parallelizing of sims across cores ####
for(w in whichScenario){# which scenario of the 24
  
  ##### prepping the cores ####
  doParallel::registerDoParallel(cores = batches)
  
  paste0("SPIM_simResults_scenario",w,
         "_N",parm_combos$N.inds[w],
         "_p0",parm_combos$p0[w],
         "_coh",parm_combos$cohesion[w],
         "_agg",parm_combos$aggregation[w],
         "_",paste(whichIDS,collapse=""))
  
  #subset the data to the scenario
  datasets<-sim_data[[w]]
  
  #run each dataset for that scenario in parallel
  #.export = c("simsFromThisCore")
  out = foreach(i = 1:batches ,.packages = c("coda","SPIM","dplyr","tidyr")) %dopar% {
    
    set.seed(i)#i
    
    # we dont want the core to just do 1 sim
    #we want it to do multiple;
    #if we do 20 batches/cores, with  3 sims per batch
    # one after the other
    #ie the 1st ,21th , 31, 41 41st on the first core
    
    test1<-myFun2(data=datasets[[i]],whichIDS=whichIDS,niters=niters,burn=burn,j=i) 
    test2<-myFun2(data=datasets[[i+batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+batches)
    test3<-myFun2(data=datasets[[i+2*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+2*batches)
    test4<-myFun2(data=datasets[[i+3*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+3*batches)
    test5<-myFun2(data=datasets[[i+4*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+4*batches)
    
    test6<-myFun2(data=datasets[[i+5*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+5*batches)
    test7<-myFun2(data=datasets[[i+6*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+6*batches)
    test8<-myFun2(data=datasets[[i+7*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+7*batches)
    test9<-myFun2(data=datasets[[i+8*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+8*batches)
    test10<-myFun2(data=datasets[[i+9*batches]],whichIDS=whichIDS,niters=niters,burn=burn,j=i+9*batches)

    test<-list(test1,test2,test3,test4,test5,
               test6,test7,test8,test9,test10)
  }
  
  #put everything in a large dataframe instead of a nested list
   out <- unlist(out, recursive = FALSE)
   out<-unlist(out,recursive=FALSE)
   # acceptRate, gel, summary, gammaProps, time, itsTot, PID

   summary_df <- do.call("rbind", out[which(names(out)=="summary")])
   summary_df$param<-rep(rownames(out$summary),max(summary_df$sim))
  
   # 
   gel_df<-unlist(out[which(names(out)=="gel")],recursive = FALSE)
   gel_df<-as.data.frame(do.call("rbind",gel_df[which(names(gel_df)=="gel.psrf")]))
   colnames(gel_df)<-c("Rhat_pt","Rhat_upper")
   summary_df<-cbind(summary_df,gel_df)
  # 
   
   accept_df<-unlist(out[which(names(out)=="acceptRate")],recursive = FALSE)
   neff_df<-unlist(out[which(names(out)=="eff")],recursive = FALSE)
   summary_df$accept<-accept_df
   summary_df$effSize<-neff_df
   
   #
   summary_df$its<-rep(do.call("rbind", out[which(names(out)=="itsTot")]),each=5)
   #
   time<-do.call("rbind", out[which(names(out)=="time")])
   time<-time[rep(seq_len(nrow(time)), each = 5), ]
   
   summary_df<-cbind(summary_df,time)
   summary_df$scenario<-w
   summary_df$M<-M
   
   # colnames(summary_df)
   # [1] "Mean"           "SD"             "Naive SE"       "Time-series SE"
   # [5] "2.5%"           "25%"            "50%"            "75%"           
   # [9] "97.5%"          "sim"            "param"          "Rhat_pt"       
   # [13] "Rhat_upper"     "accept"         "effSize"        "its"           
   # [17] "time"           "scenario"       "M"             
   #  
   summary_df<-summary_df[,c(18,11,10, 1,7,2,5,9,16,19,17,12:15)]
  
  write.csv(summary_df,paste0("SPIM_simResults_scenario",w,
                              "_N",parm_combos$N.inds[w],
                              "_p0",parm_combos$p0[w],
                              "_coh",parm_combos$cohesion[w],
                              "_agg",parm_combos$aggregation[w],
                              "_",paste(whichIDS,collapse=""),
                              ".csv"))
  #which(summary_df$Rhat_pt>1.2)
   
   gamma_df <- do.call("rbind", out[which(names(out)=="gammaProps")])
   gamma_df$scenario<-w
   write.csv(gamma_df,paste0("SPIM_gammaResults_scenario",w,
                               "_N",parm_combos$N.inds[w],
                               "_p0",parm_combos$p0[w],
                               "_coh",parm_combos$cohesion[w],
                               "_agg",parm_combos$aggregation[w],
                               "_",paste(whichIDS,collapse=""),
                               ".csv"))
   
   PID_list<-out[which(names(out)=="PID")]
   saveRDS(PID_list,paste0("SPIM_PID_scenario",w,
                           "_N",parm_combos$N.inds[w],
                           "_p0",parm_combos$p0[w],
                           "_coh",parm_combos$cohesion[w],
                           "_agg",parm_combos$aggregation[w],
                           "_",paste(whichIDS,collapse=""),
                           ".rds"))
   
  stopImplicitCluster()
}

