###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### TEOMR function (Temperature Effect On Metabolic Rates)
#### Version v1.0
#### 09.04.2019
#### Author : Planque, Lindstrøm, Sivel
###############################################################################################

##### TEOMR function
# Based on the metabolic theory of ecology
# 1K results in a change of MR by 10%

TEOMR<-function(NDNDData,TDS,deltaT=0){
  # NDNDData is the main NDND data object
  # TDS is the list of temperature dependent species (a vector of zeroes and ones)
  # deltaT is the change in temperature in Kelvin (+ or -)
  which.species=which(TDS==1)
  NDNDData$Mu[which.species]=NDNDData$Mu[which.species]*(1.1^deltaT)
  NDNDData$Rho[which.species]=NDNDData$Rho[which.species]*(1.1^deltaT)
  NDNDData$Sgma[which.species]=NDNDData$Sgma[which.species]*(1.1^deltaT)
  return(NDNDData)
}

  
  