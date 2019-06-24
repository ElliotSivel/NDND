###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model simulation
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

# 0. Warning --------------------------------------------------------------

# The code needs to run parts in C++ to increase the speed of simulation
# Therefore, the some functions are coded in c++ and need to be compiled
# If you are working with Linux or Mac, the compiler is included, you don't need to include any other package
# If you are working on Windows, Rtools needs to be downloaded and installed on the computer
# https://cran.r-project.org/bin/windows/Rtools/
# Once installed you need to add the the path to Rtools to the environment path with the function sys.setenv()

#####
# Needs to be generalized
Sys.setenv(PATH = "C:/Users/a22073/Documents/Rtools/bin;
           C:/Users/a2207/Documents/R-3.6.0/bin/x64;
           C:/Program Files(x86)/Common Files/Oracle/Java/javapath;
           C:/WINDOWS/system32;
           C:/WINDOWS; 
           C:/WINDOWS/System32/Wbem;
           C:/WINDOWS/System32/WindowsPowerShell/v1.0/;
           C:/Program Files/Intel/WiFi/bin/;
           C:/Program Files/Common Files/Intel/WirelessCommon/;
           C:/Program Files/Git/cmd;
           C:/Users/Administrator/AppData/Local/Microsoft/WindowsApps;
           C:/Users/Administrator/AppData/Local/Box/Box Edit/")

# 1. Initialization ----------------------------------------------------------

# Initialize the simulation : Clear your work environment

graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

# Load specific libraries

library(abind)                              # Loading the abind package -- Bindind arrays
library(devtools)                           # Loading the devtools package -- Tools for package development
library(ggplot2)                            # Loading the ggplot package -- Graphics package
library(inline)                             # Loading the inline package -- Read C++ functions called in R
library(LIM)                                # Loading the LIM package -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)
library(linprog)                            # Loading the linprog package -- Linear programming and optimization
library(memisc)                             # Loading the memisc package -- 
library(plyr)                               # Loading the plyr package -- Dataframe manipulation
library(pracma)                             # Loading the pracma package -- Mathematical tools and functions
library(Rcpp)                               # Loading the Rcpp package -- Integration of R and C++
library(RcppEigen)                          # Loading the RcppEigen package -- Linear algebra on C++ integrated in R
library(rstudioapi)                         # Loading the rstudioapi package -- Directory and file interactive choice
library(tidyverse)                          # Loading the tidyverse package -- Set of packages for data manipulations

# 2. Load functions -------------------------------------------------------

setwd(selectDirectory(caption = "Select your work directory : "))                  # Chose your initial work directory
wd<-getwd()                                                                        # Save initial work directory

NDNDfunctions<-list.files(paste(wd,"/ndnd_functions",sep=""))                      # List the functions in the folder ndnd_functions
# Source the functions in the folder ndnd_functions
for (i in 1:length(NDNDfunctions)){
  function2source<-paste(wd,"/ndnd_functions","/",NDNDfunctions[i],sep="")
  source(function2source)
}

# 3. Load directories -------------------------------------------------

dir_tab<-directories(wd=wd)                    # Save all directories available from the initial directory in a vector

# 4. Load data ---------------------------------------------------

# The Data file is not varying from unless the initial data has been changed
# If a NDNDData file has already been saved for a precise configuration you don't need to rerun the entire dataprep

# Question :
test_ndnd_data<-askYesNo("Do you have already a NDNDData file for the configuration setups you wish?")

# If you have the NDNDData file you can directly load it
if(test_ndnd_data==F){
  source('./ndnd_code/Dataprep.r')
} else {load(file=paste(dir_tab$data_dir,"/NDNDData.RData",sep=""))}

# 5. run the simulation ------------------------------------------------------

# Here it run a the simulation for one simulation

T1<-Sys.time()                                       # Timer on
NDNDSimulation=SimNDND(NDNDData)
T2<-Sys.time()                                       # Timer off

Tdif<-difftime(T2,T1)                                # Simulation time
Tdif

# 6. save simulation outputs  -------------------------------------------

save(NDNDSimulation,file=paste(NDNDData$directories$outputs_dir,"/NDNDSim_",format(NDNDSimulation$Simulation.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))

# 7. Plots ---------------------------------------------------------------

Fig<-plot.NDND(NDNDSimulation = NDNDSimulation)
Fig[[2]]
