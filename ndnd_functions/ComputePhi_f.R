###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Compute Phi function
#### Version v1.0
#### 24.06.19
#### Author : Benjamin Planque
#### License ???
###############################################################################################

##### ComputePhi function
# Implementation of a simple harvest control rule to the fishing mortality value

# 1.6 million km2 for the BS
# Bmp for demersals = 1 million tonnes = 0.625t/km2
# Fmp for demersals = 0.4y-1
# Bmp for pelagics = 1 million tonnes
# Fmp for pelagics = 0.4y-1

ComputePhi <- function(NDNDData,CurrentBiomass){
  Fmp=NDNDData$Fmp                       # maximum fishing mortality
  Bmp=NDNDData$Bmp                       # Biomass trigger (from Management Plan)
  slopes=Fmp/Bmp                         # slope of the biomass:F relationship
  Fslope=slopes*CurrentBiomass           # fishing mortality along the slope
  Phi=apply(rbind(Fmp,Fslope),2,min)     # fishing mortality following management plan
  Phi[is.na(Phi)]=0                      # To avoid the simulation to crash because of NA, replace NA by 0
  return(Phi)
}
  