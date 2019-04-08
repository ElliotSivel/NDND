###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model simulation
#### Version v1.1.5
#### 08.04.19
#### Author : Elliot Sivel
###############################################################################################

# 1. Initialization ----------------------------------------------------------

# Initialize the simulation : Clear your work environment

graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

# Load specific libraries

require(LIM)                                           # loading the library LIM -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)
require(ggplot2)                                       # ggplot

# Set a time and date tag for the flags
Sys.time.NDND=Sys.time()

# Set your work directories

### to set in function
# Choose interactively your work directory
c_info<-Sys.info()[1]
if (c_info == "Windows") {
  setwd(choose.dir(caption = "Choose main working directory (/NDND)"))
} else {setwd(tcltk::tk_choose.dir(caption = "Choose main working directory (/NDND)"))}
                                        
wd<-getwd()                                                 # Get the new work directory and save your work directory in your work environment for later
config_dir=paste(wd,'/ndnd_config',sep="")                                             # Sets the directory to the folder where data files are located
configreport_dir=paste(wd,'/ndnd_configreport',sep="")                                             # Sets the directory to the folder where data files are located
data_dir=paste(wd,'/ndnd_data',sep="")                                             # Sets the directory to the folder where data files are located
files_dir=paste(wd,'/ndnd_files',sep="")                                             # Sets the directory to the folder where data files are located
functions_dir=paste(wd,'/ndnd_functions',sep="")                                             # Sets the directory to the folder where data files are located
outputs_dir=paste(wd,'/ndnd_outputs',sep="")                                             # Sets the directory to the folder where data files are located


# 2. Source the functions -------------------------------------------------

## source functions
NDNDfunctions<-list.files(functions_dir)
for (i in 1:length(NDNDfunctions)){
  paste('./ndnd_functions/',NDNDfunctions[i],sep="")
  source(paste('./ndnd_functions/',NDNDfunctions[i],sep=""))
}

# 3. load files --------------------------------------------------------------

# There are 6 .txt files needed to run the model
# Data files (5) : species.txt, fluxes.txt, coefs.txt, import.txt, export.txt
# Configuration file : NDNDConfig -- Names of the data files, length of simulation and plotting parameter

# Load the Configuration file
# Choice of the file is set as interactive

setwd(config_dir)                                                            # Sets the directory to the folder where configuration files are located
config_file<-file.choose()             # Opens a window to choose the configuration file your want to implement.

# Loading the data with readDATA

NDNDData<-readDATA(config_file=config_file,files_dir=files_dir,Sys.time.NDND = Sys.time.NDND)                       # Applies the readDATA function

# Saving the NDNDData.RData in a particular file

save(NDNDData,file=paste(data_dir,"/NDNDData_",paste(format(Sys.time.NDND,"%Y_%m_%d_%H_%M_%S"),sep = ""),".RData",sep=""))             # Save NDNDData

# 4. Report Model Configuration --------------------------------------------

# Create a report of the input data and parametrization

setwd(configreport_dir)                                  # Set the directory to the folder where configuration reports are saved 
NDNDConfigreport(NDNDData = NDNDData)                         # Applies the NDNDConfigreport function
                                                              # Creates the configuration report

# 5. Compute A : Matrix of constraint on flows ----------------------------

# NDND model is based on a random sampling of flows in a restricted range of possibilities
# Range of possibilites restricted by constraints expressed as inequalities
# Ax =< b
# A is the set of constraints on the fluxes
# b is the set of constraints set on biomasses
# x is a matric of fluxes and sampled randomly

# Compute matrix of constraints of flows -- Constant over the entire simulation

setwd(outputs_dir)
A<-ComputeA(NDNDData)

# 6. Initialization of the simulation -------------------------------------

# Creates the elements in which the flows, new biomasses are saved

BiomassSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = NDNDData$ns)               # Creates a matrix of dimension Tmax vs number of species full of 0. It is thought be filled with the biomass obtained during the calculations
FlowSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = sum(NDNDData$PFv))            # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
BiomassSeries[1,]<-NDNDData$Biomass                                                     # We give that the biomass for t = 1 is the initial biomass
#sim_b<-NULL #to be removed?
#sim_pb<-NULL #to be removed?

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
  
  #b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu,NDNDData$nn)
  b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu,NDNDData$nn)
  #sim_b<-cbind(sim_b,b)#to be removed?
  
  # Delete irrelevant constraints using possibleAb
  
  Abp<-possibleAb(A,b,NDNDData$PFv)
  pA<-Abp[[1]];pb<-Abp[[2]]               # Defines two matrices with the existing constraint on flows and biomasses
  #sim_pb<-cbind(sim_pb,pb)#to be removed?
  
  ##### POSSIBLE FUNCTION
  # The LIM defines three possible constraints
  # For inequalities, the package formulates it Gx >= H
  # Transformation for sampling
  
  G<-as.matrix(-pA)                       # Defining A (constraints on fluxes) for the sampling
  h<-as.matrix(-pb)                       # Defining b (constraints on biomasses) for th sampling
  
  # Sampling Fluxes -- Sampling algorithm : mirror
  # TRY=try()
  # if(inherits(TRY, "try-error")==FALSE){
  # } else {
  # }  
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

Fig<-Plot.NDND(BiomassSeries,NDNDData$Tmax,NDNDData$Species,NDNDData$ns,NDNDData$Plotting)
Fig

