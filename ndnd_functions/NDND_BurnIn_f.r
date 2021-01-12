###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Burn-In function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# This function deletes the burn-in from the output data
NDND.BI<-function(AutoCor,Out.file){
  
  # extract BI value from the AutoCor file
  BI<-AutoCor$Burn_In
  Biomass<-Out.file$Biomass
  Flows<-Out.file$Flows
  Fish<-Out.file$Fish[,1:9]
  Tmax<-Out.file$Data$Tmax
  sim_l<-as.numeric(length(table(Biomass[,ncol(Biomass)])))
  
  Biomass_BI<-NULL
  Flows_BI<-NULL
  Fish_BI<-NULL
  
  for (i in 1:sim_l){
    a<-Biomass[as.matrix(Biomass$Simulation)==as.matrix(1:sim_l)[i],]
    b<-as.matrix(a[(BI+1):Tmax,])
    Biomass_BI<-rbind(Biomass_BI,b)
    
    c<-Flows[as.matrix(Flows$Simulation)==as.matrix(1:sim_l)[i],]
    d<-as.matrix(c[(BI+1):(Tmax-1),])
    Flows_BI<-rbind(Flows_BI,d)
    
    e<-Fish[as.matrix(Fish[,9])==as.matrix(1:sim_l)[i],]
    f<-as.matrix(e[(BI+1):(Tmax-1),])
    Fish_BI<-rbind(Fish_BI,f)
  }
  
  # Gather the time series of biomass, flows and fisheries into one list
  g<-list(Biomass=Biomass_BI,Flows=Flows_BI,Fish=Fish_BI) 
  
  # gather the two list in one
  out<-list(Sim.tag=Output$Sim.tag,
            Data=Output$Data,
            Burn.In=BI,
            Output=g)
  
  return(out)
}

# # Delete the Burn-in in the output with the incomplete time series
# a<-Output$with_NaN$Biomass[BI:Output$Data$Tmax,,]
# b<-Output$with_NaN$Flows[BI:(Output$Data$Tmax-1),,]
# c<-Output$with_NaN$Fisheries[BI:(Output$Data$Tmax-1),,]
# 
# # Gather the time series of biomass, flows and fisheries into one list
# d<-list(Biomass=a,Flows=b,Fisheries=c)
# 
# # Delete the Burn-in in the output with only complete time series
# e<-Output$without_NaN$Biomass[BI:Output$Data$Tmax,,]
# f<-Output$without_NaN$Flows[BI:(Output$Data$Tmax-1),,]
# g<-Output$without_NaN$Fisheries[BI:(Output$Data$Tmax-1),,]