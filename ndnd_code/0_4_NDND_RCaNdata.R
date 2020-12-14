###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### RCaN Data
#### Version v1.0
#### 14.12.20
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

# Load RCaN simulation outputs and transform in the NDND format to use them for further analysis 

load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_code/RCaN_Elliot_1e6_samples.Rdata") # Load the RCaN files

# Two objects are contained in the RCaN outputs: myCaNmod and myCaNmodFit
# Vizualize the structure of each objects
str(myCaNmodFit)
str(myCaNmod)

obs<-do.call(rbind, myCaNmodFit) # Rbind the two objects

# The obs object has 806 columns, from 625, there are the biomass data, before it's the flow data 
sub_obs<-obs[,625:ncol(obs)]/1600 # Separate the biomass data
f_obs<-obs[,1:624]/1600 # Separate the flow data

# Format the data to have one column for one species
observation<-NULL
for (i in 1:nrow(sub_obs)){
  print(i)
  ss<-sub_obs[i,]
  ss<-matrix(ss,ncol = 7)
  observation<-rbind(observation,ss)
}

# Format the data to have one column per flow
obs_flows<-NULL
for (i in 1:nrow(f_obs)){
  print(i)
  ss<-f_obs[i,]
  ss<-matrix(ss,ncol = 24,byrow = T)
  obs_flows<-rbind(obs_flows,ss)
}

# The RCaN consider fishing landings as a flow, as such we need to identify the fishing flows and remove them
fish<-obs_flows[,c(18,22)]
obs_flows<-obs_flows[,as.numeric(rownames(myCaNmod$fluxes_def[myCaNmod$fluxes_def$Trophic==TRUE,]))] # Consider only the trophic flows
obs_flows<-obs_flows[,c(1:5,10,8,9,6,7,11,12,16,14,15,13,17,18)] # Re-arrange the flows in the order of the NDND 

# Save the biomass and the flow data
save(observation,file=paste(NDNDData$directories$data_dir,"/RCaN_Obs.RData",sep = ""))
save(obs_flows,file=paste(NDNDData$directories$data_dir,"/RCaN_flows.RData",sep = ""))