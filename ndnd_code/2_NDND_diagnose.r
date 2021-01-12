###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Output Diagnose
#### Version v1.0
#### 03.12.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################


# 0. Libraries ------------------------------------------------------------

# Needed Libraries
###
require(tidyverse)
require(ggplot2)
###

# Source the initialization file, functions, directories, etc,...
source(file = "C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_code/0_1_NDND_initialization.R") # Used to load all functions in the work environment

# 1. Load outputs ---------------------------------------------------------

load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/NDNDsim_out.RData") # Load the latest simulation
# load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/NDNDSim_2020_11_18_19_11_20.RData") # Load the chose simulation after a specific time tag

s_length<-length(unique(Output$Biomass$Simulation))

pdf(file = paste(Output$Data$directories$outputs_dir,"/Diagnostic/NDNDdiag_",format(Output$Sim.tag,"%Y_%m_%d_%H_%M_%S"),".pdf",sep = ""),onefile = T)
ndnddiagnose(Output)
dev.off()



