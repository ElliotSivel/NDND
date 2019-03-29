###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model simulation
#### Version v1.1.4
#### 29.03.19
#### Author : Elliot Sivel
###############################################################################################

# 1. Initialization ----------------------------------------------------------

# Initialize the simulation : Clear your work environment

rm(list=ls())                                                                      #clear the work environment

# Load specific libraries

library(LIM)                                           # loading the library LIM -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)

# Set a time and date tag for the flags

date_time_name<-paste(format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),"_",sep = "")         # Creates a tag to identify the run in file names

# Set your work directory

setwd(choose.dir())                                         # Choose interactively your work directory
wd<-getwd()                                                 # Get the new work directory and save your work directory in your work environment for later

# 2. Source the functions -------------------------------------------------

setwd('./ndnd_functions')                             # Sets the directory to the folder where function files are located
source("readDATA_f.r")                                # Load the readDATA function -- Read the data and create a list of all data NDNDData
source("NDNDConfigreport_f.r")                        # Load the NDNDConfigreport function -- Creates a report of the input data and parameters
source("ComputeA_f.r")                                # Load the ComputeA function -- Computes the matrix of constraints on the flows
source("Computeb_f.r")                                # Load the Computeb function -- Computes the matrix of constraints on the biomasses
source("possibleAb_f.r")                              # Load the possibleAb function -- Deletes the irrelevant constraints in A and b matrices
source("ComputeBiomass_f.r")                          # Load the ComputeBiomass function -- Computes the biomass at time step t+1
setwd(wd)

# 3. load files --------------------------------------------------------------

# 6 .txt files + the function files needed to run the model
# Data files : species.txt, fluxes.txt, coefs.txt, import.txt, export.txt
# Configuration file : NDNDConfig -- Names of the data files, length of simulation and plotting parameter

# Load the Configuration file
# Choice of the file is set as interactive

setwd('./ndnd_config')                                                            # Sets the directory to the folder where configuration files are located
NDNDconfig<-choose.files(caption = "Choose your configuration file!")             # Opens a window to choose the configuration file your want to implement.
setwd(wd)

# Loading the data with readDATA

setwd('./ndnd_files')                                             # Sets the directory to the folder where data files are located
NDNDData<-readDATA(NDNDconfig = NDNDconfig)                       # Applies the readDATA function
setwd(wd)

# Saving the NDNDData.RData in a particular file

setwd("./ndnd_data")                                                       # Set the directory to the folder where NDNDData is saved
save(NDNDData,file=paste("NDNDData_",date_time_name,".RData"))             # Save NDNDData
setwd(wd)

# 4. Report Model Configuration --------------------------------------------

# Create a report of the input data and parametrization

setwd("./ndnd_configreport")                                  # Set the directory to the folder where configuration reports are saved 
NDNDConfigreport(NDNDData = NDNDData)                         # Applies the NDNDConfigreport function
                                                              # Creates the configuration report
setwd(wd)

# 5. Compute A : Matrix of constraint on flows ----------------------------

# NDND model is based on a random sampling of flows in a restricted range of possibilities
# Range of possibilites restricted by constraints expressed as inequalities
# Ax =< b
# A is the set of constraints on the fluxes
# b is the set of constraints set on biomasses
# x is a matric of fluxes and sampled randomly

# Compute matrix of constraints of flows -- Constant over the entire simulation

setwd('./ndnd_outputs/ndnd_computeA')
A<-ComputeA(Gama=NDNDData$Gama,Kapa = NDNDData$Kapa,ns=NDNDData$ns,nn=NDNDData$nn)
setwd(wd)

# 6. Initialization of the simulation -------------------------------------

# Creates the elements in which the flows, new biomasses are saved

BiomassSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = NDNDData$ns)               # Creates a matrix of dimension Tmax vs number of species full of 0. It is thought be filled with the biomass obtained during the calculations
FlowSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = sum(NDNDData$PFv))            # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
BiomassSeries[1,]<-NDNDData$Biomass                                                     # We give that the biomass for t = 1 is the initial biomass
sim_b<-NULL
sim_pb<-NULL

## Not implemented yet, needs to be done for further investigation
# T=1                                                                                     # Set initial time to 1 
# Tcrash=0                                                                                # Set Tcrash to 0

# 9. Main loop -----------------------------------------------------------------

# The model run over Tmax years

for (t in 1:(NDNDData$Tmax-1)) {
  print(t)                                        # We want to be able to see the state of the simulation
  
  # Select the import and export values for the simulation year in Importall and Exportall
  
  Import<-NDNDData$Importall[t,]
  Export<-NDNDData$Exportall[t,]
  
  # Compute b
  
  b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu,NDNDData$nn)
  sim_b<-cbind(sim_b,b)
  
  # Delete irrelevant constraints using possibleAb
  
  Abp<-possibleAb(A,b,NDNDData$PFv)
  pA<-Abp[[1]];pb<-Abp[[2]]               # Defines two matrices with the existing constraint on flows and biomasses
  sim_pb<-cbind(sim_pb,pb)
  
  ##### POSSIBLE FUNCTION
  # The LIM defines three possible constraints
  # For inequalities, the package formulates it Gx >= H
  # Transformation for sampling
  
  G<-as.matrix(-pA)                       # Defining A (constraints on fluxes) for the sampling
  h<-as.matrix(-pb)                       # Defining b (constraints on biomasses) for th sampling
  
  # Sampling Fluxes -- Sampling algorithm : mirror
  
  Fsample<-xsample(G=G,H=h,iter=100,burninlength=100,outputlength=100,type="mirror")

  #####
  
  # Fsample = list of four elements
  # First is sampled flows
  # We defined a number of iterations (100) -- We only need one vector
  
  F0<-Fsample[[1]][sample(1:nrow(Fsample[[1]]),1),]                   # Sample one random vector of flows among 100 sampled with the xsample function
  
  # Reattributing the flows values
  
  F<-rep(0,NDNDData$nn)                                               # Creating a vector of 0 and of length nn
  F[NDNDData$PFv==1]<-F0                                              # Attribute the flow values at the right place according to the vector PFv
  
  FlowSeries[t,]<-F0
  
  # Compute biomass at the next time step
  
  BiomassSeries[t+1,]=ComputeBiomass(Biomass = BiomassSeries[t,],F=F,Import = Import,Export = Export,Gama = NDNDData$Gama,Mu = NDNDData$Mu,Kapa = NDNDData$Kapa,ns = NDNDData$ns)
}


# 10. Plots ---------------------------------------------------------------



