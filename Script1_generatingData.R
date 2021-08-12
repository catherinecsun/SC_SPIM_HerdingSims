# annotated R code based on code from Bischof et al. 2020 Group paper
# - added option of multiple sampling occasions,
# which makes it so that the cohesion for an individual can change per occasion
# ie fission-fusion, with strengh determined by the cohesion parameter. 
# - adding partial ids to inds

# we assume that there is no spatial pattern in adhesion and cohesion of partial ids

# bischof simulated nonindpeendent data and then used the independent model
# we can do the same for SC and SPIM  but it would also be cool to 
# then use a nonindpendent model, written in jags, 
# to estimate at least cohesion and maybe adhesion (Group size?)

library(secr)
library(raster)
library(sp)
library(extraDistr)
library(dplyr)
library(ggplot2)

#bischof pop simulation parameters
#sigma=1.5
#p0=0.1
#n=128 (ie 0.32 per unit)
#cohesion<-c(0,0.25, 0.5, 0.75, 1)
#aggregation<-c(1,2,4,8,16,32,64,128)
#resulting in 40 scenarios, each simulated and run 1000x.

# pop parms for us
grid.size.x<-5 #traps in the x direction
grid.size.y<-15 #traps in the y direction
cell.size<- 3# unit spacing between traps
habitat.buffer<-9
nocc<-4
nsims<-100

sigma<-3
N.inds<-140
p0<-c(0.05,0.2) # low, good,
cohesion<-c(0,0.3, 0.67, 1)
aggregation<-c(1,4,10) # max group size detected in algar
#partial ids: sex + antlers;  + collar;  + collar + coat
#resulting in 2*4*3*3=72 scenarios


parm_combos<-expand.grid(grid.size.x,grid.size.y,cell.size,
                         habitat.buffer,nocc,
                         sigma,N.inds,
                         p0,cohesion, aggregation,
                         nsims)
colnames(parm_combos)<-c("grid.size.x","grid.size.y","cell.size",
                      "habitat.buffer","nocc",
                      "sigma","N.inds",
                      "p0","cohesion", "aggregation",
                      "nsims")


#ID covariates per indivdiual 
parms_IDcovs<-list(antlers_cont=c(0.015,0.180,0.036,0.015,0.064,0.039,0.049,0.052,0.090,0.078,0.053,0.058,0.055,0.015,0.015,0.036,0.037,0.004),
                   sex_bin= c(0.27,0.73), #  F/M
                   collar_bin=c(0.96,0.04), #no collar/ collar
                   coat_bin=rep(1/5,5))# black, white, cinnamon, blonde, piebald

####functions #####
#to add all-zero capture histories back in to allow alignment between capture histories of all group members,
expandCH <- function(CH) {# Function from M. Efford (pers. comm., 5/13/2019).
  pop <- attr(CH, 'pop')
  fullCH <- array(0, dim = c(nrow(pop), dim(CH)[2:3]))
  rownames(fullCH) <- rownames(pop)
  fullCH[rownames(CH),,] <- CH
  traps(fullCH) <- traps(CH)
  class(fullCH) <- 'capthist'
  fullCH
}

#create one data replicate
createOne <- function (parms,parms_IDcovs) {
  # Function modified/expanded from version shared by M. Efford (pers. comm., 5/13/2019)
  out <- list()
  detectpar <- list(g0 = parms$p0, sigma = parms$sigma)
  
  # detector grid
  tr <- make.grid(nx = parms$grid.size.x,ny = parms$grid.size.y,
                  spacing = parms$cell.size, detector = 'proximity')
  
  # available habitat
  msk <- make.mask(traps = tr, buffer = parms$habitat.buffer)
  
  # expected n etc.
  p <- sum(pdot(X = msk, traps = tr,
                detectfn = 'HN', detectpar = detectpar, noccasions = 1) *
             attr(msk, 'area')) / maskarea(msk)
  
  # clone at level of capthist for maximal cohesion (identical AC location and identical detection pattern)
  #maximal cohesion via the number of inds to sim, Nbuffer
  #rows correspond to individuals/groups
  pop <- sim.popn(Nbuffer = parms$N.inds/parms$group.size, core = tr,
                  buffer = parms$habitat.buffer,
                  Ndist = "fixed")
  
  # id covariates per individual 
  pop_ids<-data.frame(Ind=seq(1,parms$N.inds,1),
                      Group_ID=paste0(sort(rep(1:(parms$N.inds/parms$group.size),parms$group.size)),
                                      ".",rep(1:(parms$N.inds/parms$group.size),parms$group.size)))
  for(idcov in 1:length(parms_IDcovs)){
    pop_ids[,idcov+2]<-rcat(parms$N.inds,parms_IDcovs[[idcov]])
    colnames(pop_ids)[idcov+2]<- strsplit(names(parms_IDcovs)[idcov],"_")[[1]][1]
  }
  
  #change the largest category for "continuous" covs to 999, to represent 0 
  #because 0 will later be used to represent missing covariate values.
  whichCovs_cont<-which(sapply(strsplit(names(parms_IDcovs),"_"),"[[",2)=="cont")
  for(w in 1:length(whichCovs_cont)){
    #are there any with the largest value?
    actuallyMissing<-which(pop_ids[,w+2]==length(parms_IDcovs[[w]]))
    #if so,
    if(length(actuallyMissing)>0){
      pop_ids[,w+2][]<-999 
    }
  }
  
  #capt hist, and retain the number(which group, which ind) of the pop from which they were drawn
  ch <- sim.capthist(tr, popn = pop, detectfn = 'HN', renumber = FALSE,
                     detectpar = detectpar, noccasions = parms$noccasions, savepopn = TRUE)
  
  out$capthist_aggCoh<-ch
  # now clone  popn for minimal cohesion (just identical AC location)
  #since  above was for all groups, clone the pop to explicitly get all inds 
  clonedpop <- secr::clone(pop, 'constant', parms$group.size)
  
  #detected inds in groups
  # so that it's n x K x J
  chg <- sim.capthist(tr, popn = clonedpop, detectfn = 'HN', renumber = FALSE,
                      detectpar = detectpar, noccasions = parms$noccasions, savepopn = TRUE)
  
  #add back in the non-detected individuals into the capture history
  #so that it's N x K X J
  chg <- expandCH(chg)
  
  # mix individual and group detection patterns
  # first, create the detections to draw from if individual separated from group
  chi <- secr::clone(expandCH(ch), 'constant', parms$group.size)
  ch <- chi
  
  #binary 0/1 according to cohesion parameter to determine if indivdiual
  #follows the group or follows its own activity pattern
  # with only 1 occasion, an ind is only ever entirely cohesive or not at all.
  # with mutliple occasions, creates fission-fusion dynamics, whereby
  #the stronger the cohesion, the more consistent across sampling occasions the mixture will be
  #the weaker the cohesion, ind movements will be more often differ from the group's.
  #so intermediate values will generate strong fission-fusion
  p <- runif(prod(dim(ch))) < parms$cohesion # Bernoulli by animal and detector
  ch[] <- p * chi + (1-p) * chg # mixture of group and individual detection patterns
  ch2 <- subset(ch) # drops null histories
  
  #capture histories
  
  out$parms<-parms
  out$parms_IDcovs<-parms_IDcovs
  out$detectpar<-detectpar
  out$clonedpop<-clonedpop
  out$traps<-tr
  out$mask<-msk
  out$pop_locs<-pop
  out$pop_ids<-pop_ids
  out$expectation <- c(p = p, exp.n = p * parms$N.inds,
                       sd.n = parms$N.inds* (p * (1-p) / parms$N.inds)^0.5)
  out$capthist<-ch2
  return(out)
  
}
# in this function,  out$capthist and out$capthist_aggCoh
# may differ, with the former having detections of more inds
# that is because some indivdiuals, compared to out$capthist_aggCoh,
# now stray from the group and can be detected elsewhere by themselves
# ie., the former (out$capthist) contains all the detections
# in the latter(out$capthist_aggCoh) but with indivdiual ids 
# and then also more detections of others.
#the number of more detections depends on the 
# cohesion parameter of that scenario


#### now use all the functions to simulate data ##### 
sim_data<-list()# list to put the sims per scenario

for(r in 1:nrow(parm_combos)){ #
  print(paste0("scenario ", r))
  #gather all the parameters and conditions
  tmp_parms <- list(grid.size.x = parm_combos$grid.size.x[r],
                grid.size.y = parm_combos$grid.size.y[r],
                cell.size = parm_combos$cell.size[r],
                habitat.buffer = parm_combos$habitat.buffer[r],
                group.size = parm_combos$aggregation[r],
                N.inds = parm_combos$N.inds[r],
                p0 = parm_combos$p0[r],
                noccasions=parm_combos$nocc[r], # this is new compared to Bischof
                sigma = parm_combos$sigma[r],
                cohesion = parm_combos$cohesion[r])
  
  tmp_list<-list()
  #create the s simulated datasets for each r scenario
  for(s in 1:parm_combos$nsims[r]){
    print(paste0("sim ",s))
    tmp_sim<-createOne(parms = tmp_parms,parms_IDcovs=parms_IDcovs)
    
    #put it into a growing list of simulated datasets for this scenario
    tmp_list[[s]]<-tmp_sim
    names(tmp_list)[[s]]<-paste0("sim_",s)
    
  }
  
  #put all the simulated datasets for that scenario into a growing list
  sim_data[[r]]<-tmp_list
  names(sim_data)[[r]]<-paste0("scenario",r,"_",tmp_parms$N.inds,"N_",tmp_parms$p0,"det_",
         tmp_parms$cohesion,"coh_",tmp_parms$group.size,"agg")
  
}




# PLOTTING simulated data
whichSim<-1
for(w in 1:nrow(parm_combos)){
  print(w)
  data_tmp<-sim_data[[w]][[whichSim]]
  
  tr.sp <- data.frame(data_tmp$traps)
  coordinates(tr.sp) <- tr.sp
  # habitat
  mask.sp <- data.frame(data_tmp$mask)
  coordinates(mask.sp) <- mask.sp
  e <- extent(mask.sp)
  e.sp <- as(e, 'SpatialPolygons')
  # ACs without jitter
  AC.orig.sp <- data.frame(data_tmp$clonedpop)
  coordinates(AC.orig.sp) <- AC.orig.sp
  # ACs with jitter
  AC.sp <- data.frame(data_tmp$clonedpop)
  AC.sp$ID <- dimnames(data_tmp$clonedpop)[[1]]
  jitter.f <- data_tmp$parms$cell.size/10*sqrt(data_tmp$parms$group.size)
  coordinates(AC.sp) <- apply(AC.sp[,1:2],c(1,2),function(x)
    runif(1,x-jitter.f,x+jitter.f))
  
  jpeg(paste0("plot_simData_","N",data_tmp$parms$N.inds,
              "_p0",data_tmp$parms$p0,"_sig",data_tmp$parms$sigma,
              "_agg", data_tmp$parms$group.size,
              "_coh",data_tmp$parms$cohesion,
              "_nocc",data_tmp$parms$noccasions,".jpg"))
  
  # PLOT HABITAT EXTENT
  plot(e.sp,pch = 19,col = NA,border = "black")
  # PLOT DETECTIONS
  ids.detected <- dimnames(data_tmp$capthist)[[1]]
  AC.colors <- sample(colors(),length(AC.sp),replace = TRUE)
  names(AC.colors) <- AC.sp$ID
  for(i in 1:length(ids.detected)){
    x <- data_tmp$capthist[i,,]
    locs <- data.frame(data_tmp$traps[which(colSums(x)> 0),])
    coordinates(locs) <- locs
    xx <- coordinates(locs)
    yy <- coordinates(AC.sp[AC.sp$ID == ids.detected[i],])
    segments(yy[,1],yy[,2],xx[,1],xx[,2],col = AC.colors[i],lwd = 0.8)
  }
  # PLOT ACs
  plot(AC.sp,pch = 21,bg = AC.colors,add = TRUE,col = "black",lwd = 0.2,cex = 0.7)
  plot(tr.sp,col = grey(0.4),cex = 0.7,add = TRUE)
  title(main=paste0("N = ",data_tmp$parms$N.inds,
                    "; p0: ",data_tmp$parms$p0,"; sigma: ",data_tmp$parms$sigma,
                    "\n group size (aggregation): ", data_tmp$parms$group.size,
                    "; cohesion: ",data_tmp$parms$cohesion,
                    "\n sampling occasions: ",data_tmp$parms$noccasions),
        sub=paste0("(",length(ids.detected), " inds detected at ",
                   length(which(colSums(apply(data_tmp$capthist,c(1,3),sum))>0)),
                   " of ", nrow(tr.sp)," sites)"),line=0.1)
  dev.off()
}


saveRDS(sim_data, file = "simulatedData_from24Scenarios.rds")
# Restore the object
readRDS(file = "my_data.rds")

# }



## on July 12 I changed the sex ratio from the skewed c(0.27,0.73) to an even c(0.5,0.5)
# and then replaced this partial covariate ID in the sim data
# and also fixed the odd phenomenon of some sims having all inds with 999 for antler counts
# the antler issue must be where all antlers in a sim were changed to 999 
# if there was even 1 with the last/18th antler category 
# let's not worry abtou considering missing partial id covs in these sims, and just call it a new category
# but on July 30, I noticed the antler proportions sum to <1.
# so i fixed this, using the mean posterior values from 'outputSPIM_caribou_Algar_2019_all_wCollars_postProcessed_Gammas.RData' 
#  as.vector(round(colMeans(antlers_chainsInOne),3))
#  and added 0.001 to the first value (0.021) to make it 0.022 so that the sum of the proportions would be 1
sum(sim_data[[1]][[1]]$parms_IDcovs$antlers_cont)

View(sim_data)
sex_bin_new=c(0.5,0.5) #  F/M c(0.27,0.73)
antlers_cont_new<-c(0.022,0.185,0.043,0.021,0.070,0.045,0.055,0.059,0.096,0.085,0.060,0.064,0.061,0.021,0.021,0.043,0.043,0.006)

for(r in 1: length(sim_data)){
  for(s in 1:parm_combos$nsims[1]){
  sim_data[[r]][[s]]$parms_IDcovs$sex_bin<-sex_bin_new
  sim_data[[r]][[s]]$pop_ids$sex<-rcat(sim_data[[r]][[s]]$parms$N.inds,sex_bin_new)
  
  sim_data[[r]][[s]]$parms_IDcovs$antlers_cont<-antlers_cont_new
  sim_data[[r]][[s]]$pop_ids$antlers <-rcat(sim_data[[r]][[s]]$parms$N.inds,antlers_cont_new)
  }
}

#check to make sure that there arent any funky sims with 999 for antler counts anymore
errors<-c()
for(r in 1: length(sim_data)){
  for(s in 1:parm_combos$nsims[1]){
   
    if(length(which(sim_data[[r]][[s]]$pop_ids$antlers==999))==140){
      errors<-c(errors,paste0("scenario ",r," sim ", s))
    }
  }
}



saveRDS(sim_data, file = "simulatedData_from24Scenarios.rds")

## summary statistics  ####
# about the individual-detection data, which will be latent
# to be able to say things like 
# we detected, on average, about 40% (range = 15–60%) of the true population size N under any given scenario

sim_data_Stats<-data.frame(scenario=NA,sim=NA,
                           p0=NA,aggregation=NA, cohesion=NA,
                           ninds=NA, ntraps=NA)
for(p in 1:length(sim_data)){
  print(p)
  for (s in 1:length(sim_data[[p]])){
    tmp_data<-sim_data[[p]][[s]]
    
    # number of detected individuals
    tmp_ninds<-dim(tmp_data$capthist)[1]
    # number of traps that detected individuals
    tmp_ntraps<-length(which(colSums(apply(tmp_data$capthist,c(2,3),sum))>0))
    
    sim_data_Stats<-rbind(sim_data_Stats,
                          c(p,s,
                            tmp_data$parms$p0,tmp_data$parms$group.size,tmp_data$parms$cohesion,
                            tmp_ninds,tmp_ntraps))
  }
}
sim_data_Stats<-sim_data_Stats[-1,]

sim_data_StatsSum<-sim_data_Stats%>%group_by(scenario)%>%
  summarize(#p0,aggregation, cohesion,
    min(ninds),
    mean(ninds),
    max(ninds),
    sd(ninds),
    min(ntraps),
    mean(ntraps),
    max(ntraps),
    sd(ntraps))
sim_data_StatsSum$p0<-parm_combos$p0
sim_data_StatsSum$cohesion<-parm_combos$cohesion
sim_data_StatsSum$aggregation<-parm_combos$aggregation


###plots 

p0.labs <- c("p0: 0.05", "p0: 0.20")
names(p0.labs) <- c("0.05", "0.2")

coh.labs <- c("Cohesion: 0","Cohesion: 0.3", "Cohesion: 0.67", "Cohesion: 1")
names(coh.labs) <- c("0","0.3", "0.67", "1")

plot_ninds<-ggplot(data=sim_data_Stats,aes(x=as.factor(aggregation),y=ninds))+
  geom_boxplot(aes(fill=as.factor(p0)))+
  ylim(c(0,N.inds))+
  facet_grid(cols=vars(cohesion), #rows=vars(p0),
             labeller = labeller(cohesion=coh.labs))+ #switch = "y" p0 = p0.labs,
  labs(title="Number of Detected Individuals",
       x="Aggregation (Group Size)",y="n",fill="Detection \nProbability")+theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5), 
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
plot_ninds

plot_ntraps<-ggplot(data=sim_data_Stats,aes(x=as.factor(aggregation),y=ntraps))+
  geom_boxplot(aes(fill=as.factor(p0)))+
  ylim(c(0,grid.size.x*grid.size.y))+
  facet_grid(cols=vars(cohesion), #rows=vars(p0),
             labeller = labeller(cohesion=coh.labs))+ #switch = "y" p0 = p0.labs,
  labs(title="Number of Traps with Detections",
       x="Aggregation (Group Size)",y="n",fill="Detection \nProbability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
plot_ntraps


#this is in the original bischof code, to run the scr model and then calc stats.
# and included these arguments in the function: plot.check = TRUE,do.fit = TRUE
#
# MODEL FITTING to normal Scr model
# if(do.fit){
#   out$fit <- out$error <- NULL
#   try(out$fit <- summary(secr.fit(ch2, mask = msk, trace = F,
#                                   details = list(distribution = "binomial")))
#       ,silent = TRUE)
# }
