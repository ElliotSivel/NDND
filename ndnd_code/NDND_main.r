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
library(tidyverse)
library(pracma)
library(linprog)
library(RcppEigen)
library(Rcpp)
library(inline)
library(devtools)
library(memisc)

# 2. Load functions -------------------------------------------------------

setwd(selectDirectory(caption = "Select your work directory : "))
wd<-getwd()

NDNDfunctions<-list.files(paste(wd,"/ndnd_functions",sep=""))
for (i in 1:length(NDNDfunctions)){
  function2source<-paste(wd,"/ndnd_functions","/",NDNDfunctions[i],sep="")
  source(function2source)
}

# 3. Load directories -------------------------------------------------
## source functions

dir_tab<-directories(wd=wd)

# 4a. prepare data file ---------------------------------------------------

source('./ndnd_code/Dataprep.r')

# 4b. load data file --------------------------------------------------------------

load(file="NDNDData.Rdata")

# 5. run the simulation ------------------------------------------------------

T1<-Sys.time()
NDNDSimulation=SimNDND(NDNDData)
T2<-Sys.time()

Tdif<-difftime(T2,T1)
Tdif

# 6. save simulation outputs  -------------------------------------------

save(NDNDSimulation,file=paste(NDNDData$directories$outputs_dir,"/NDNDSim_",format(NDNDSimulation$Simulation.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))

# 7. Plots ---------------------------------------------------------------

Fig<-plot.NDND(NDNDSimulation = NDNDSimulation)
Fig[[2]]

