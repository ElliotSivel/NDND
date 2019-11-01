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
library(arrayhelpers)                       # Loading the arrayhelpers package -- Functions for array modification and shaping
library(cluster)
library(cpgsR)
library(dendextend)
library(devtools)                           # Loading the devtools package -- Tools for package development
library(factoextra)
library(ggdendro)                           # Loading the ggdendro package -- Package to plot dendrograms with ggplot
library(ggplot2)                            # Loading the ggplot package -- Graphics package
library(inline)                             # Loading the inline package -- Read C++ functions called in R
library(LIM)                                # Loading the LIM package -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)
library(linprog)                            # Loading the linprog package -- Linear programming and optimization
library(memisc)                             # Loading the memisc package -- 
library(NbClust)                            # Loading the NbClust package -- Package to define the optimal number of cluster in a dataset
library(plyr)                               # Loading the plyr package -- Dataframe manipulation
library(pracma)                             # Loading the pracma package -- Mathematical tools and functions
library(Rcpp)                               # Loading the Rcpp package -- Integration of R and C++
library(RcppEigen)                          # Loading the RcppEigen package -- Linear algebra on C++ integrated in R
library(rgl)                                # Loading the rgl package -- 3d graphs
library(rstudioapi)                         # Loading the rstudioapi package -- Directory and file interactive choice
library(scales)                             # Loading the scales package -- Graphical attributes package
library(svglite)
library(tidyverse)                          # Loading the tidyverse package -- Set of packages for data manipulations

# 2. Initialization -------------------------------------------------------

setwd(selectDirectory(caption = "Select your work directory : "))                  # Chose your initial work directory
wd<-getwd()                                                                        # Save initial work directory

s_length<-as.numeric(readline("How many trajectories?    "))                       # Set the number of simulations

# 3. Load functions -------------------------------------------------------

NDNDfunctions<-list.files(paste(wd,"/ndnd_functions",sep=""))                      # List the functions in the folder ndnd_functions
# Source the functions in the folder ndnd_functions
for (i in 1:length(NDNDfunctions)){
  function2source<-paste(wd,"/ndnd_functions","/",NDNDfunctions[i],sep="")
  source(function2source)
}

# 4. Load directories -------------------------------------------------

dir_tab<-directories(wd=wd)                    # Save all directories available from the initial directory in a vector

# 5. Load data ---------------------------------------------------

# The Data file is not varying from unless the initial data has been changed
# If a NDNDData file has already been saved for a precise configuration you don't need to rerun the entire dataprep

# Question :
test_ndnd_data<-askYesNo("Existing NDNDData file with wished setups?")

# If you have the NDNDData file you can directly load it
if(test_ndnd_data==F){
  source('./ndnd_code/Dataprep.r')
} else {load(file=paste(dir_tab$data_dir,"/NDNDData.RData",sep=""))}

# 6. run the simulation ------------------------------------------------------

# # Here it run a the simulation for one simulation
# 
# T1<-Sys.time()                                       # Timer on
# NDNDSimulation=SimNDND(NDNDData)
# T2<-Sys.time()                                       # Timer off
# 
# Tdif<-difftime(T2,T1)                                # Simulation time
# Tdif

# Make it run for several simulations

sim_bio<-NULL                 # Create the object that will contain the biomass simulation outputs
sim_flow<-NULL                # Create the object that will contain the flows simulation outputs
sim_phi<-NULL                 # Create the object that will contain the Fishing losses simulation outputs

for (i in 1:s_length){
  print(i)
  sim<-SimNDND(NDNDData)
  sim_bio<-rbind(sim_bio,sim$Output$BiomassSeries)                     # Binds the matrices of biomass time series in a thrid dimension 
  sim_flow<-rbind(sim_flow,sim$Output$FlowSeries)                      # Binds the matrices of flows time series in a third dimension
  sim_phi<-rbind(sim_phi,sim$Output$PhiSeries)                         # Binds the matrices of the fishing losses in a third dimension
}

sim_out<-NDND.nb.sim(Data = NDNDData,
                     s_length = s_length,
                     Biomass = sim_bio,
                     Flows = sim_flow,
                     Phi = sim_phi)
sim_out_NaN<-NDND.NaN.to.zero(simulation_output = sim_out)                          # We only keep the simulations for which we have the total number of flows 
Output<-list(Sim.tag=sim$Simulation.tag,Data=sim$Data,Output = sim_out, Output_NaN_0 = sim_out_NaN)

# 6. save simulation outputs  -------------------------------------------

save(Output,file=paste(NDNDData$directories$outputs_dir,"/NDNDSim_",format(sim$Simulation.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
save(Output,file=paste(NDNDData$directories$outputs_dir,"/NDNDsim_out.RData",sep = ""))

# Load the output data

load(file=paste(NDNDData$directories$outputs_dir,"/NDNDsim_out.RData",sep = ""))

## ==> Output exploration