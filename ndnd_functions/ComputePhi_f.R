# 1.6 million km2 for the BS
# Bmp for demersals = 1 million tonnes = 0.625t/km2
# Fmp for demersals = 0.4y-1
# Bmp for pelagics = 1 million tonnes
# Fmp for pelagics = 0.4y-1

ComputePhi <- function(NDNDData,CurrentBiomass){
  Fmp=NDNDData$Fmp  # maximum fishing mortality
  Bmp=NDNDData$Bmp # Biomass trigger (from Management Plan)
  slopes=Fmp/Bmp # slope of the biomass:F relationship
  Fslope=slopes*CurrentBiomass # fishing mortality along the slope
  Phi=apply(rbind(Fmp,Fslope),2,min) # fishing mortality following management plan
  Phi[is.na(Phi)]=0
  return(Phi)
}
  