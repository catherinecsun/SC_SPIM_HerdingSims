# Calculating probablity of identity
# for partial identities

library(gtools) #for combinations()
library(plyr) #for rbind.fill()
library(ggplot2) #for ggplot()

#the (corrected) partial ids
parms_IDcovs<-sim_data[[1]][[1]]$parms_IDcovs

#4 covariates (all 4); 1 combo
# antlers, sex, collar, coat
combos_4<-as.data.frame(combinations(4,4,names(parms_IDcovs)))
combos_4$probs_unique<-NA

for(com in 1:nrow(combos_4)){
  probs<-c() #temporary housing
  
  #for each combination of 4 covariates
  #step 1a: find all unique Identity possibilities 
  
  #names(parms_IDcovs) which(names(parms_IDcovs)%in%combos_4[com,])
  partial_IDsprobs<-expand.grid(parms_IDcovs[[1]],
                                parms_IDcovs[[2]],
                                parms_IDcovs[[3]],
                                parms_IDcovs[[4]])
  #there are this may identity possibilities
  nrow(partial_IDsprobs)
  
  #step 1b: get probability of each identity 
  partial_IDsprobs$probability<-apply(partial_IDsprobs, 1, prod)
  
  #step 2: what's the probability that 2 independent draws 
  # will yield the same identity?
  
  # multinomial: [n!/(x_1! ... x_k!)] * p_1^x1 ... p_k^xk
  # the first term in [] is equal to 2!/(2!*however many 0!'s)= 1
  #the second term reduces to the square of the probabiltiy of the focal identity
  
  for(identity in 1:nrow(partial_IDsprobs)){
    probs<-c(probs,partial_IDsprobs$probability[identity]^2)
  }
  
  sum(probs) #total probability of drawing 2 of the same identities, whatever that identity is.
  
  1-sum(probs) #probability that 2 draws will yield 2 different identities
  combos_4$probs_unique[com]<-(1-sum(probs))
}

combos_4<-cbind(nIDs=4,combos_4)
colnames(combos_4)<-c("nIDs","Antlers","Coat","Collar","Sex","Probs_unique")
combos_4[,2:5]<-1



####3 covariates; 4 combos####
combos_3<-as.data.frame(combinations(4,3,names(parms_IDcovs)))
combos_3$probs_unique<-NA

for(com in 1:nrow(combos_3)){
  probs<-c() #temporary housing
  
  #for each combination of 4 covariates
  #step 1a: find all unique Identity possibilities 
  
  whichPartialIDS<-which(names(parms_IDcovs)%in%combos_3[com,])
  partial_IDsprobs<-expand.grid(parms_IDcovs[[whichPartialIDS[1]]],
                                parms_IDcovs[[whichPartialIDS[2]]],
                                parms_IDcovs[[whichPartialIDS[3]]])
  #there are this may identity possibilities
  nrow(partial_IDsprobs)
  
  #step 1b: get probability of each identity 
  partial_IDsprobs$probability<-apply(partial_IDsprobs, 1, prod)
  
  #step 2: what's the probability that 2 independent draws 
  # will yield the same identity?
  
  # multinomial: [n!/(x_1! ... x_k!)] * p_1^x1 ... p_k^xk
  # the first term in [] is equal to 2!/(2!*however many 0!'s)= 1
  #the second term reduces to the square of the probabiltiy of the focal identity
  
  for(identity in 1:nrow(partial_IDsprobs)){
    probs<-c(probs,partial_IDsprobs$probability[identity]^2)
  }
  
  
  sum(probs) #total probaility of drawing 2 of the same identities
  
  1-sum(probs) #probaility that 2 draws will yield 2 different identities
  combos_3$probs_unique[com]<-(1-sum(probs))
}

combos_3<-cbind(nIDs=3,combos_3)

combos_3final<-matrix(NA,nrow=nrow(combos_3),ncol=6)
colnames(combos_3final)<-c("nIDs","Antlers","Coat","Collar","Sex","Probs_unique")
combos_3final[,1]<-3
combos_3final[,6]<-combos_3$probs_unique
combos_3final<-as.data.frame(combos_3final)
for(com in 1:nrow(combos_3)){
  combos_3final$Antlers[com]<-ifelse("antlers_cont"%in%combos_3[com,2:4],1,0)
  combos_3final$Sex[com]<-ifelse("sex_bin"%in%combos_3[com,2:4],1,0)
  combos_3final$Collar[com]<-ifelse("collar_bin"%in%combos_3[com,2:4],1,0)
  combos_3final$Coat[com]<-ifelse("coat_bin"%in%combos_3[com,2:4],1,0)
  
  
}
combos_3<-combos_3final

####2 covariates; 6 combos#####
combos_2<-as.data.frame(combinations(4,2,names(parms_IDcovs)))
combos_2$probs_unique<-NA
for(com in 1:nrow(combos_2)){
  
  probs<-c() #temporary housing
  
  #for each combination of 4 covariates
  #step 1a: find all unique Identity possibilities 
  
  whichPartialIDS<-which(names(parms_IDcovs)%in%combos_2[com,])
  partial_IDsprobs<-expand.grid(parms_IDcovs[[whichPartialIDS[1]]],
                                parms_IDcovs[[whichPartialIDS[2]]])
  #there are this may identity possibilities
  nrow(partial_IDsprobs)
  
  #step 1b: get probability of each identity 
  partial_IDsprobs$probability<-apply(partial_IDsprobs, 1, prod)
  
  #step 2: what's the probability that 2 independent draws 
  # will yield the same identity?
  
  # multinomial: [n!/(x_1! ... x_k!)] * p_1^x1 ... p_k^xk
  # the first term in [] is equal to 2!/(2!*however many 0!'s)= 1
  #the second term reduces to the square of the probabiltiy of the focal identity
  
  for(identity in 1:nrow(partial_IDsprobs)){
    probs<-c(probs,partial_IDsprobs$probability[identity]^2)
  }
  
  sum(probs) #total probaility of drawing 2 of the same identities
  
  1-sum(probs) #probaility that 2 draws will yield 2 different identities
  combos_2$probs_unique[com]<-(1-sum(probs))
}

combos_2<-cbind(nIDs=2,combos_2)

combos_2final<-matrix(NA,nrow=nrow(combos_2),ncol=6)
colnames(combos_2final)<-c("nIDs","Antlers","Coat","Collar","Sex","Probs_unique")
combos_2final[,1]<-2
combos_2final[,6]<-combos_2$probs_unique
combos_2final<-as.data.frame(combos_2final)
for(com in 1:nrow(combos_2)){
  combos_2final$Antlers[com]<-ifelse("antlers_cont"%in%combos_2[com,2:3],1,0)
  combos_2final$Sex[com]<-ifelse("sex_bin"%in%combos_2[com,2:3],1,0)
  combos_2final$Collar[com]<-ifelse("collar_bin"%in%combos_2[com,2:3],1,0)
  combos_2final$Coat[com]<-ifelse("coat_bin"%in%combos_2[com,2:3],1,0)
  }
combos_2<-combos_2final

####1 covariate; 4 combos####
combos_1<-as.data.frame(combinations(4,1,names(parms_IDcovs)))
combos_1$probs_unique<-NA

for(com in 1:nrow(combos_1)){
  probs<-c() #temporary housing
  
  #for each combination of 4 covariates
  #step 1a: find all unique Identity possibilities 
  
  whichPartialIDS<-which(names(parms_IDcovs)%in%combos_1[com,])
  partial_IDsprobs<-expand.grid(parms_IDcovs[[whichPartialIDS[1]]])
  #there are this may identity possibilities
  nrow(partial_IDsprobs)
  
  #step 1b: get probability of each identity 
  partial_IDsprobs$probability<-apply(partial_IDsprobs, 1, prod)
  
  #step 2: what's the probability that 2 independent draws 
  # will yield the same identity?
  
  # multinomial: [n!/(x_1! ... x_k!)] * p_1^x1 ... p_k^xk
  # the first term in [] is equal to 2!/(1!1!0!0!)= 2
  #the second term reduces to the square of the probabiltiy of the focal identity
  
  for(identity in 1:nrow(partial_IDsprobs)){
    probs<-c(probs,(partial_IDsprobs$probability[identity]^2))
  }
  
  sum(probs) #total probaility of drawing 2 of the same identities
  
  1-sum(probs) #probaility that 2 draws will yield 2 different identities
  combos_1$probs_unique[com]<-(1-sum(probs))
}

combos_1<-cbind(nIDs=1,combos_1)

combos_1final<-matrix(NA,nrow=nrow(combos_1),ncol=6)
colnames(combos_1final)<-c("nIDs","Antlers","Coat","Collar","Sex","Probs_unique")
combos_1final[,1]<-1
combos_1final[,6]<-combos_1$probs_unique
combos_1final<-as.data.frame(combos_1final)

for(com in 1:nrow(combos_1)){
  combos_1final$Antlers[com]<-ifelse("antlers_cont"%in%combos_1[com,2],1,0)
  combos_1final$Sex[com]<-ifelse("sex_bin"%in%combos_1[com,2],1,0)
  combos_1final$Collar[com]<-ifelse("collar_bin"%in%combos_1[com,2],1,0)
  combos_1final$Coat[com]<-ifelse("coat_bin"%in%combos_1[com,2],1,0)
}
combos_1<-combos_1final

###all together ####
combos<-rbind.fill(combos_1,combos_2,combos_3,combos_4)


## we considered
# 4: antler, sex, coat, collar
# 3: sex, collar, coat
# 2: sex, coat
combos[c(15,14,9),]

colors <- c("Antlers" = "red", "Coat" = "blue",
            "Sex"="green","Collar" = "gray")
shapes<-c("Combinations"=19,"Evaluated"=0)

ggplot()+
  scale_color_manual(values = colors)+
  scale_shape_manual(values=shapes)+
  stat_summary(data=combos[combos$Antlers==1,],
               aes(x=nIDs,y=Probs_unique,color="Antlers"),
              fun = mean, geom = "line") +
  stat_summary(data=combos[combos$Coat==1,],
               aes(x=nIDs,y=Probs_unique,color="Coat"),
               fun = mean, geom = "line") +
  stat_summary(data=combos[combos$Sex==1,],
               aes(x=nIDs,y=Probs_unique,color="Sex"),
               fun = mean, geom = "line") +
  stat_summary(data=combos[combos$Collar==1,],
               aes(x=nIDs,y=Probs_unique,color="Collar"),
               fun = mean, geom = "line")+
  geom_point(data=combos[c(15,14,9),],aes(x=nIDs, y=Probs_unique,shape="Evaluated"),size=4)+
  geom_point(data=combos,aes(x=nIDs, y=Probs_unique,shape="Combinations"),size=3)+
  ylim(c(0,1))+
  labs(title="Theoretical Expected Probability of Identity ",
       x="# of Partial Identity Covariates",y="Probability",
       color = "Mean",
       shape ="")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


