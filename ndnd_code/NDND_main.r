###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model simulation
#### Version v1.1.5
#### 08.04.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# 1. Initialization ----------------------------------------------------------

# Initialize the simulation : Clear your work environment

graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

# Load specific libraries

library(rstudioapi)                                    # Loading the library for directory and file interactive choice
library(LIM)                                           # Loading the library LIM -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)
library(ggplot2)                                       # Loading the ggplot package -- graphics package 


# 2a. prepare data file ---------------------------------------------------

#source('./ndnd_code/Dataprep.r')

# 2b. load data file --------------------------------------------------------------

load(file="NDNDData.Rdata")

# 3. Source the functions -------------------------------------------------
## source functions

NDNDfunctions<-list.files(NDNDData$directories$functions_dir)
for (i in 1:length(NDNDfunctions)){
  function2source<-paste(NDNDData$directories$functions_dir,"/",NDNDfunctions[i],sep="")
  source(function2source)
}

# 4. run the simulation ------------------------------------------------------

NDNDSimulation=SimNDND(NDNDData)

# 5. save simulation outputs  -------------------------------------------

save(NDNDSimulation,file=paste(NDNDData$directories$outputs_dir,"/NDNDSim_",format(NDNDSimulation$Simulation.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))

# 6. Plots ---------------------------------------------------------------

#Fig<-plot.NDND(BiomassSeries,NDNDData$Tmax,NDNDData$Species,NDNDData$ns,NDNDData$Plotting)
#Fig

