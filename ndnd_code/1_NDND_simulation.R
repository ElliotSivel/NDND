###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model simulation
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################


# 0. Libaries -------------------------------------------------------------

# Needed libraries
###
require(yesno)
require(abind)
###

# Source the initialization file
source(file = "C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_code/0_2_NDND_Data.R")

# 1. Simulation -----------------------------------------------------------

### We test if an output is already available
test_ndnd_Output<-yesno2("Do you have a NDNDSim file ?", yes = "y", no="n")

if (test_ndnd_Output==F){ # If simulation outputs not available, we run the simulation
  
  s_length<-as.numeric(readline("How many trajectories?    "))                       # Set the number of simulations
  
  # Estimate the time of simulation, T1 & T2, T1-T2 method
  T1<-Sys.time()                                # Timer on
  
  # Define empty arrays that are going to be filled during the simulation
  sim_bio<-NULL # Create the object that will contain the biomass simulation outputs
  sim_flow<-NULL # Create the object that will contain the flows simulation outputs
  ExpF<-NULL
  Tcrash<-NULL # Create an empty vector where the number of crashes per run is saved
  AllA<-NULL
  ALLsimb<-NULL
  
  # Run the simulation
  for (i in 1:s_length){
    print(i)
    
    #####
    
    sim<-SimNDND(NDNDData)                                               # SIMULATION
    
    #####
    
    # Filling the output objects
    sim_bio<-rbind(sim_bio,sim$Output$BiomassSeries)                     # Binds the matrices of biomass time series in a thrid dimension 
    sim_flow<-rbind(sim_flow,sim$Output$FlowSeries)                      # Binds the matrices of flows time series in a third dimension
    ExpF<-rbind(ExpF,sim$Output$Fish)
    Tcrash[i]<-sim$Tcrash                                                # Fills the Tcrash vector
    AllA<-abind(AllA,sim$A,along = 3)
    ALLsimb<-abind(ALLsimb,sim$b, along=3)
  }
  
  T2<-Sys.time()                                       # Timer off
  Tdif<-difftime(T2,T1)                                # Simulation time
  
  # For the rest of the computation, we need to transform NaN into 0
  sim_bio[is.na(sim_bio)]<-0
  sim_flow[is.na(sim_flow)]<-0
  
  # Gather outputs in one list
  Output<-list(Sim.tag = sim$Simulation.tag,               # Saves the output with a unique ID based on the date and time of simulation
               Data = sim$Data,                            # Saves the data used to perform the simulation -- Copy of the NDNDData file
               Tcrash = Tcrash,                            # Saves the Tcrash table -- gives an information of the quality of the simulation
               Biomass = sim_bio,                          # Saves the Biomass time series
               Flows = sim_flow,                           # Saves the Flows time series
               Fish = ExpF,                              # Saves the fishing mortality time series
               Code = sim$Code,
               A=AllA,
               b=ALLsimb)                            # Saves the code used to run the simulation -- Code + functions
  # )
  
  # 4. save simulation outputs  -------------------------------------------
  
  save(Output,file=paste(NDNDData$directories$outputs_dir,"/NDNDSim_",format(sim$Simulation.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
  save(Output,file=paste(NDNDData$directories$outputs_dir,"/NDNDsim_out.RData",sep = ""))
  
} else {load(file=paste(NDNDData$directories$outputs_dir,"/NDNDsim_out.RData",sep = ""))} # If there is a simulation output available, we load it

