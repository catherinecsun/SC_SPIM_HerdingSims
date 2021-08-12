#Run SC model on data

#libraries
library(doParallel)
library(coda)
library(nimble)

##### The data ####
#scenario 1:8 all have aggregation of 1 so the cohesion is irrelevant
#what changes is the detection probability
# so 1, 3, 5, and 7 are identical (p0 of 0.05 but coohesoon of 0,0.3, 0.67 and 1)
# and2 4, 6, and 8 are the same
#so can remove 3:8

View(sim_data)
sim_data<-sim_data[-c(3:8)]
parm_combos<-parm_combos[-c(3:8),]

##### how many batches of sims to do/cores to use? ####
batches<-50

##### which scenarios? ####
whichScenario<-c(14:18)

##### how many sims per batches to do? ####
#actually i have hardwired this below; not ideal but it gets the job done for now
#sims<-2
#batches*sims #total number of sims to do

##### how many iterations to run and burn per simulation? ####
niters<-16666*3 #30000
burn<-4998#5000


##### the function(s) to do the modeling ####
#the model needs to compile inside the function because it uses the data
#not ideal because having to compile for every dataset takes a while
myFun<-function(data, niters,burn){ #,j
  
  #nimble prep
  
  #the static things
  M<-400
  X<-data$traps
  xlim <- range(data$mask[,1]) 
  ylim <- range(data$mask[,2])
  
  #the particular data
  simSC<-t(apply(data$capthist,c(2,3),sum))
  
  #nimble
  nimbleconstants <- list(M = M,
                          X=X,J = nrow(simSC),
                          K=ncol(simSC),
                          xlim=xlim, ylim=ylim, 
                          area=((xlim[2]-xlim[1])*(ylim[2]-ylim[1])))
  inits <- list(sigma=rnorm(1,3),lam0=runif(1), z=rep(1,M),psi=0.5)
  nimbledata <- list( n = simSC)
  params <- c("N", "D", "lam0", "sigma","psi")#,"z","s")
  
  
  nimbleOptions(showCompilerOutput = TRUE)
  time_a<-Sys.time()
  model <- nimbleModel(code = sc_nimble, name = "sc_nimble",
                       constants = nimbleconstants, data = nimbledata,
                       inits = inits)
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  time_b<-Sys.time()
  out1 <- runMCMC(CmodelMCMC, niter = niters,samples=TRUE, summary=TRUE)
  time_c<-Sys.time()
  
  #burn and split into 3 chains to check for convergence
  out1$samples<-out1$samples[-c(1:burn),]
  
  chain1<-as.mcmc(out1$samples[1:((niters-burn)/3),colnames(out1$samples)%in%params])
  chain2<-as.mcmc(out1$samples[(1+(niters-burn)/3):((2*(niters-burn)/3)),colnames(out1$samples)%in%params])
  chain3<-as.mcmc(out1$samples[(1+(2*(niters-burn)/3)):nrow(out1$samples),colnames(out1$samples)%in%params])
  chains<-as.mcmc.list(list(chain1,chain2,chain3))
  
  
  gel<-gelman.diag(chains,autoburnin = FALSE,multivariate=FALSE)
  out1$gel<-gel
  
  #clean up the summary
  out1$summary<-out1$summary[rownames(out1$summary)%in%params,]
  out1$summary<-as.data.frame(out1$summary)
  out1$summary$param<-rownames(out1$summary)
  out1$summary$sim<-i#j
  
  #add the elapsed times, total iterations
  out1$time<-c(time_b-time_a,time_c-time_b)
  out1$itsTotal<-niters-burn
  
  #keep only the gel, summary,  elapsed time, and total iterations
  out2<-out1[2:5]
  
  return(out2)
  
  
}

myFun2<-function(data, niters,burn,j){ #
  
  #nimble prep
  
  #the static things
  M<-400
  X<-data$traps
  xlim <- range(data$mask[,1]) 
  ylim <- range(data$mask[,2])
  
  #the particular data
  simSC<-t(apply(data$capthist,c(2,3),sum))
  
  #nimble
  nimbleconstants <- list(M = M,
                          X=X,J = nrow(simSC),
                          K=ncol(simSC),
                          xlim=xlim, ylim=ylim, 
                          area=((xlim[2]-xlim[1])*(ylim[2]-ylim[1])))
  inits <- list(sigma=rnorm(1,3),lam0=runif(1), z=rep(1,M),psi=0.5)
  nimbledata <- list( n = simSC)
  params <- c("N", "D", "lam0", "sigma","psi")#,"z","s")
  
  
  nimbleOptions(showCompilerOutput = TRUE)
  time_a<-Sys.time()
  model <- nimbleModel(code = sc_nimble, name = "sc_nimble",
                       constants = nimbleconstants, data = nimbledata,
                       inits = inits)
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  time_b<-Sys.time()
  out1 <- runMCMC(CmodelMCMC, niter = niters,samples=TRUE, summary=TRUE)
  time_c<-Sys.time()
  
  #burn and split into 3 chains to check for convergence
  out1$samples<-out1$samples[-c(1:burn),]
  
  chain1<-as.mcmc(out1$samples[1:((niters-burn)/3),colnames(out1$samples)%in%params])
  chain2<-as.mcmc(out1$samples[(1+(niters-burn)/3):((2*(niters-burn)/3)),colnames(out1$samples)%in%params])
  chain3<-as.mcmc(out1$samples[(1+(2*(niters-burn)/3)):nrow(out1$samples),colnames(out1$samples)%in%params])
  chains<-as.mcmc.list(list(chain1,chain2,chain3))
  
  
  gel<-gelman.diag(chains,autoburnin = FALSE,multivariate=FALSE)
  out1$gel<-gel
  
  #clean up the summary
  out1$summary<-out1$summary[rownames(out1$summary)%in%params,]
  out1$summary<-as.data.frame(out1$summary)
  out1$summary$param<-rownames(out1$summary)
  out1$summary$sim<-j
  
  #add the elapsed time, and the iterations
  out1$time<-c(time_b-time_a,time_c-time_b)
  out1$itsTot<-niters-burn
  
  #keep only the gel, summary, elapsed times, and its
  out2<-out1[2:5]
  
  return(out2)
  
  
}


#double check that params 
params <- c("N", "D", "lam0", "sigma","psi")

##### the parallelizing of sims across cores ####
for(w in whichScenario){
  
  ##### prepping the cores ####
  doParallel::registerDoParallel(cores = batches)
  
  paste0("SC_simResults_scenario",w,
         "_N",parm_combos$N.inds[w],
         "_p0",parm_combos$p0[w],
         "_coh",parm_combos$cohesion[w],
         "_agg",parm_combos$aggregation[w])
  
  #subset the data to the scenario
  datasets<-sim_data[[w]]
  
  #run each dataset for that scenario in parallel
  #.export = c("simsFromThisCore")
  out = foreach(i = 1:batches ,.packages = c("coda","nimble")) %dopar% {
    
    sc_nimble<-nimbleCode({
      sigma ~ dgamma(24,8) # weakly informative prior - 38-619 km2
      lam0 ~ dunif(0,10) #baseline detection rate
      psi ~dbeta(1,1) # data augmentation
      
      for(i in 1:M) {
        z[i] ~ dbern(psi)
        s[i,1] ~ dunif(xlim[1], xlim[2])
        s[i,2] ~ dunif(ylim[1], ylim[2])
        
        for(j in 1:J) {
          distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
          lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i] 
        }
      } # End of 1:M
      for(j in 1:J) {
        for(k in 1:K) {
          bigLambda[j, k] <- sum(lam[1:M, j])
          n[j, k] ~ dpois(bigLambda[j, k])
          
        }
      }
      N <- sum(z[1:M])
      D <- N/area
    })
    
    set.seed(i)#i
    
    #if we were just doing 1 sim per core, then this would be enough
    #myFun(data=datasets[[i]],niters=niters,burn=burn)
    
    #but we dont want the core to just do 1 sim
    #we want it to do multiple;
    #if we do 20 batches/cores, with  3 sims per batch
    # one after the other
    #ie the 1st ,21th , and 41st on the first core
    
    test1<-myFun2(data=datasets[[i]],niters=niters,burn=burn,j=i) 
    test2<-myFun2(data=datasets[[i+batches]],niters=niters,burn=burn,j=i+batches)
    # test3<-myFun2(data=datasets[[i+2*batches]],niters=niters,burn=burn,j=i+2*batches)
    test<-list(test1,test2)#,test3)
    
    #THIS DIDNDT WORK
    #test <- vector("list", sims)
    # 
    # for(j in seq(from=i,by=batches,length.out=sims)){
    #   tmp<-myFun2(data=datasets[[j]],niters=niters,burn=burn,j=j)
    #   test[[length(test)+1]]<-tmp
    # }
    #use a new counter,j, with a starting value of i
    #and continue to do sims until batches*sims
    # j<-i
    #create a list in which to put the output of the sims
    # simsFromThisCore <- vector("list", sims)
    # while(j<(batches*sims)+1){ #(length(datasets)+1)
    #   myFun(data=datasets[[j]],niters=niters,burn=burn,j=j)
    
    #   simsFromThisCore[[length(simsFromThisCore)+1]]<- myFun(data=datasets[[j]],niters=niters,burn=burn,j=j)
    #   j<-j+batches
    # }
    # 
    
  }
  
  #put everything in a large dataframe instead of a nested list
  out <- unlist(out, recursive = FALSE)
  out<-unlist(out,recursive=FALSE)
  # summary , gel, time tot, its 
  summary_df <- do.call("rbind", out[c(TRUE,FALSE,FALSE,FALSE)])
  
  gel_df<-unlist(out[c(FALSE,TRUE,FALSE,FALSE)],recursive = FALSE)
  gel_df<-as.data.frame(do.call("rbind",gel_df[c(TRUE,FALSE)]))
  colnames(gel_df)<-c("Rhat_pt","Rhat_upper")
  summary_df<-cbind(summary_df,gel_df)
  
  summary_df$its<-rep(do.call("rbind", out[c(FALSE,FALSE,FALSE,TRUE)]),each=length(params))
  time<-do.call("rbind", out[c(FALSE,FALSE,TRUE,FALSE)])
  time<-time[rep(seq_len(nrow(time)), each = length(params)), ]
  colnames(time)<-c("time_buildCompile","time_run")
  summary_df<-cbind(summary_df,time)
  summary_df$scenario<-w
  
  colnames(summary_df)
  # [1] "Mean"              "Median"            "St.Dev."           "95%CI_low"        
  # [5] "95%CI_upp"         "param"             "sim"               "Rhat_pt"          
  # [9] "Rhat_upper"        "its"               "time_buildCompile" "time_run"         
  # [13] "scenario"         
  summary_df<-summary_df[,c(13,6,7,1:5,8:12)]
  
  write.csv(summary_df,paste0("SC_simResults_scenario",w,
                              "_N",parm_combos$N.inds[w],
                              "_p0",parm_combos$p0[w],
                              "_coh",parm_combos$cohesion[w],
                              "_agg",parm_combos$aggregation[w],
                              ".csv"))
  #which(summary_df$Rhat_pt>1.2)
  stopImplicitCluster()
}

