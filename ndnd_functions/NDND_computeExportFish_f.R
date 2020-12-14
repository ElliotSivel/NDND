###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### ExportFish function
#### Version v1.0
#### 14.12.20
#### Author : Elliot Sivel
###############################################################################################

ExportFish <- function(NDNDData,Biomass){
  
  Fmp=NDNDData$Fmp                       # Maximum fishing mortality
  Bmp=NDNDData$Bmp                       # Biomass trigger (from Management Plan
  Blim=NDNDData$Blim                     # Biomass limit for escapement
  M=NDNDData$M                           # Natural mortality
  
  CatchBiomass<-(Fmp/(Fmp+M))*(1-exp(-(Fmp+M)))         # Compute C/B after the Baranov equation
  
  Fslope<-CatchBiomass/(Bmp-Blim)                       # Compute catch on the slope, when the fishing decrease because the stock is smaller
  
  # We draw two line matching the catch depending on the size of the stock
  FBiomass<-apply(rbind((CatchBiomass*Biomass),(Fslope*(Biomass-Blim))),2,min)
  FBiomass[(Biomass<Blim)==T]<-0 # Check if the biomass is below Blim, when it is the case then set fishing to 0
  FBiomass[is.na(FBiomass)]<-0 # Check if there are NA, when it is the case, replace NA by 0
  
  return(FBiomass)
}