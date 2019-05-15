###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Main simulation function
#### Version v1.1.5
#### 09.04.19
#### Author : Elliot Sivel & Co
#### License ???
###############################################################################################

SimNDND <- function(NDNDData){
  # Set a time and date tag for the flags
  
  Simulation.tag=Sys.time()
  
  # 1. Initialization of the simulation -------------------------------------

# NDND model is based on a random sampling of flows in a restricted range of possibilities
# Range of possibilites restricted by constraints expressed as inequalities
# Ax =< b
# A is the set of constraints on the fluxes
# b is the set of constraints set on biomasses
# x is a matric of fluxes and sampled randomly

# Compute matrix of constraints of flows -- Constant over the entire simulation

A<-ComputeA(NDNDData)

# Creates the elements in which the flows, new biomasses are saved

BiomassSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = NDNDData$ns)               # Creates a matrix of dimension Tmax vs number of species full of 0. It is thought be filled with the biomass obtained during the calculations
BiomassSeries[1,]<-NDNDData$Biomass                                                     # We give that the biomass for t = 1 is the initial biomass
FlowSeries<-matrix(data = 0, nrow = (NDNDData$Tmax-1), ncol = sum(NDNDData$PFv))            # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
PhiSeries<-matrix(data = 0, nrow = (NDNDData$Tmax-1), ncol = sum(NDNDData$ns))                                   # Series of fishing mortalities

## Not implemented yet, needs to be done for further investigation
# T=1                                                                                     # Set initial time to 1 
# Tcrash=0                                                                                # Set Tcrash to 0

# 2. Main loop -----------------------------------------------------------------

# The model run over Tmax years

for (t in 1:(NDNDData$Tmax-1)) {
  # print(t)                                        # We want to be able to see the state of the simulation
  
  NDNDData$Phi=ComputePhi(NDNDData,BiomassSeries[t,]) # Compute fishing mortality
  
  b<-Computeb(NDNDData,BiomassSeries[t,],t)   # Compute b
  
  # Delete irrelevant constraints using possibleAb
  
  Abp<-possibleAb(A,b,NDNDData$PFv)
  pA<-Abp[[1]];pb<-Abp[[2]]               # Defines two matrices with the existing constraint on flows and biomasses
  
  TRY=try(Fsample<-Sampling(pA,pb,"cpgs2"))
  if(inherits(TRY, "try-error")==FALSE){
    Fsample<-Sampling(pA,pb,"cpgs2")
  } else {
    print("No polytope solution")
    break
  }
  
  # Reattributing the flows values
  
  Fluxes<-rep(0,NDNDData$nn)                                               # Creating a vector of 0 and of length nn
  Fluxes[NDNDData$PFv==1]<-Fsample                                              # Attribute the flow values at the right place according to the vector PFv
  
  # Compute biomass at the next time step
  
  BiomassSeries[t+1,]=ComputeBiomass(NDNDData,BiomassSeries[t,],Fluxes,t)
  PhiSeries[t,]=NDNDData$Phi  # store fishing mortalities
  FlowSeries[t,]<-Fsample # store trophic fluxes
  
}

colnames(BiomassSeries)<-NDNDData$Species
rownames(BiomassSeries)<-1:NDNDData$Tmax
colnames(FlowSeries)<-NDNDData$flows
rownames(FlowSeries)<-1:(NDNDData$Tmax-1)
colnames(PhiSeries)<-NDNDData$Species
rownames(PhiSeries)<-1:(NDNDData$Tmax-1)

BiomassSeries<-as.data.frame(BiomassSeries)
FlowSeries<-as.data.frame(FlowSeries)
PhiSeries<-as.data.frame(PhiSeries)

# 3. return simulation outputs  -------------------------------------------

NDNDOutput=list(BiomassSeries=BiomassSeries,FlowSeries=FlowSeries,PhiSeries = PhiSeries)
NDNDCode=NULL
maincodefile=paste(NDNDData$directories$code_dir,'/NDND_main.r',sep='')
if (file.exists(maincodefile)==TRUE){
  NDNDCode$main <- scan(maincodefile,what="",sep="\n")
  for (i in 1:length(NDNDfunctions)){
    function2scan<-paste(NDNDData$directories$functions_dir,"/",NDNDfunctions[i],sep="")
    NDNDCode$functions[[i]]=scan(function2scan,what="",sep="\n")
  }
  NDNDSimulation=list(Simulation.tag=Simulation.tag,
                      Data=NDNDData,
                      Output=NDNDOutput,
                      Code=NDNDCode)
  
  return(NDNDSimulation)
}
}

# 
#   ##### POSSIBLE FUNCTION
#   # The LIM defines three possible constraints
#   # For inequalities, the package formulates it Gx >= H
#   # Transformation for sampling
#   
#   G<-as.matrix(-pA)                       # Defining A (constraints on fluxes) for the sampling
#   h<-as.matrix(-pb)                       # Defining b (constraints on biomasses) for th sampling
#   
#   # Sampling Fluxes -- Sampling algorithm : mirror
#   # TRY=try()
#   # if(inherits(TRY, "try-error")==FALSE){
#   # } else {
#   # }  
#   Fsample<-xsample(G=G,H=h,iter=100,burninlength=100,outputlength=100,type="mirror")
#   # Fsample<-cpgs(100,-G,-h)
# 
#   #####
#   
# Fsample = list of four elements
# First is sampled flows
# We defined a number of iterations (100) -- We only need one vector

# F0<-Fsample[[1]][sample(1:nrow(Fsample[[1]]),1),]                   # Sample one random vector of flows among 100 sampled with the xsample function
# F0<-Fsample[sample(1:nrow(Fsample),1),]                   # Sample one random vector of flows among 100 sampled with the xsample function
