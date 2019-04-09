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
FlowSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = sum(NDNDData$PFv))            # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
BiomassSeries[1,]<-NDNDData$Biomass                                                     # We give that the biomass for t = 1 is the initial biomass

## Not implemented yet, needs to be done for further investigation
# T=1                                                                                     # Set initial time to 1 
# Tcrash=0                                                                                # Set Tcrash to 0

# 2. Main loop -----------------------------------------------------------------

# The model run over Tmax years

for (t in 1:(NDNDData$Tmax-1)) {
  print(t)                                        # We want to be able to see the state of the simulation
  
  # Select the import and export values for the simulation year in Importall and Exportall
  
  Import<-NDNDData$Importall[t,]
  Export<-NDNDData$Exportall[t,]
  
  # Compute b
  
  #b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu,NDNDData$nn)
  b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sgma,NDNDData$Bta,NDNDData$Mu,NDNDData$nn)
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


# 3. return simulation outputs  -------------------------------------------

NDNDOutput=list(BiomassSeries=BiomassSeries,FlowSeries=FlowSeries)
NDNDCode=NULL
maincodefile=paste(NDNDData$directories$code_dir,'/NDND_main.r',sep='')
if (file.exists(maincodefile)==TRUE){
  NDNDCode$main <- scan(maincodefile,what="",sep="\n")
  }  # reads the current file and store it into the variable 'code'
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

