# Calculating probablity of identity
# for partial identities

#the partial ids
parms_IDcovs

#possible identity combinations
prod(sapply(parms_IDcovs,length))
partial_IDs<-expand.grid(1:length(parms_IDcovs[[1]]),
                  1:length(parms_IDcovs[[2]]),
                  1:length(parms_IDcovs[[3]]),
                  1:length(parms_IDcovs[[4]]))
partial_IDsprobs<-expand.grid(parms_IDcovs[[1]],parms_IDcovs[[2]],parms_IDcovs[[3]],parms_IDcovs[[4]])

names(partial_IDs)<-names(parms_IDcovs)
names(partial_IDsprobs)<-names(parms_IDcovs)

#but they are not all equally weighted 
for(i in 1:nrow(partial_IDsprobs)){
  partial_IDsprobs$probability[i]<-partial_IDsprobs$antlers_cont[i]*partial_IDsprobs$sex_bin[i]*partial_IDsprobs$collar_bin[i]*partial_IDsprobs$coat_bin[i]
}
