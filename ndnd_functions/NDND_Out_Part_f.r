###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Adding cluster function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# The function is recreating the vector of cluster with NaN and 0
# it also merge the vector of cluster with the initial data
NDND.out.partition<-function(Data,Output,Cluster){
  
  # Extract the simulation biomass array 
  Biomass<-Output_BI$Output$Biomass
  
  a<-c(ta(Biomass))                                         # Vectorizing the biomass array
  a1<-matrix(a,ncol=dim(Biomass)[[2]],byrow=T)              # Shaping the Biomass vector in a matrix
  colnames(a1)<-Data$Species                                # Define columns names
  rownames(a1)<-c(1:nrow(a1))                               # Define rows names
  
  # Create the vector of cluster with NaN and 0
  b<-which(Biomass[,1]==0 | Biomass[,1]=="NaN")                                               # Identify the lines for which we have NaN or 0
  c<-numeric(nrow(Biomass))            # Create an empty vector of length Tmax*s_length
  
  # Filling the empty vector c
  c[b]<-NA                                                            # Fill NA where the biomass values are NaN and 0
  c[!is.na(c)]<-Output_MDS4$Partition$Cluster             # Fill the remaining of the vector with the cluster vector 
  a<-cbind(Biomass,c)                                                     # Merging the biomass matrix with cluster vector
  
  output<-list(Sim.tag=Output$Sim.tag,Biomass.cluster=a1)
  return(output)
}