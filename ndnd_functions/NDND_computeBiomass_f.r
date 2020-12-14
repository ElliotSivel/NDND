###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### ComputeBiomass function
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel
###############################################################################################

##### ComputeBiomass function
# The models samples flows randomly to compute a system dynamic under constraint
# Constraints expressed in the form of inequalities : Ax =< b
# ComputeA computes the matrix A (constraints on flows)
# Computeb computes the matrix b (constraints on biomasses)
# After sample of flows (in the main code) --> application of the model's master equation
# Need Biomass input
# Need a vector of flows
# Need Import and Export data for all species
# Need parameters values for all species for Gama (assimilation efficiency rate), Mu (metabolic losses) and Kapa (food digestability)
# Need thee number of species ns

ComputeBiomass<-function(NDNDData,CurrentBiomass,Fluxes,t){
  
  # Load parameters
  Import=NDNDData$Importall[t,]
  Export=NDNDData$Exportall[t,]
  Gama=NDNDData$Gama
  Mu=NDNDData$Mu
  # Phi=NDNDData$Phi
  Kapa=NDNDData$Kapa
  ns=NDNDData$ns
  
  # As in Computeb, Computation of constants C and D
  C=as.matrix((1-exp(-(Mu)))/(Mu))                    # Computes C
  D=as.matrix(exp(-(Mu)))                                 # Computes D
  
  # The sampling of flows is going to be done for a vector of length ns^2
  # To apply compute the new biomasses -- Need elements of length ns
  
  dim(Fluxes)<-c(ns,ns)                                      # Resahpe it in a matrix of dimension ns*ns
  FM<-t(Fluxes)                                              # Transpose F
  NewBiomass<-matrix(0,ns)                                   # Define an empty vector which length is ns
  for (i in 1:ns){                                           # For each species we compute the new biomass by applying the master equation
    NewBiomass[i]=D[i]*CurrentBiomass[i]+C[i]*(Gama[i]*sum(FM[,i]*Kapa)-sum(FM[i,])+Import[i]-Export[i])
  }
  
  return(NewBiomass)
}