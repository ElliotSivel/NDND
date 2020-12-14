###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Main simulation function
#### Version v1.1.5
#### 24.06.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

##### SimNDND function
# Runs 1 simulation of the NDND model
# Computes the constraints, samples the flows and estimates the biomass at time step t+1
# Gives times series of biomass, flows and fishing mortality

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
  FlowSeries<-matrix(data = 0, nrow = (NDNDData$Tmax-1), ncol = sum(NDNDData$PFv))        # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
  tcrash_m<-matrix(0,nrow=1,ncol=s_length)
 
  # 2. Main loop -----------------------------------------------------------------
  
  # The model run over Tmax years
  tcrash=0                           # Set the number of crashes to 0
  t=1                                # Set time step to 1
  while (t<=(NDNDData$Tmax-1)) {                    # the following is applied as long as t is lower than Tmax-1
    print(t)                                      # We want to be able to see the state of the simulation
    
    # Compute fishing mortality that will be taken from the biomass
    NDNDData$Exportall[t,]=ExportFish(NDNDData,BiomassSeries[t,])
    
    # Compute b
    b<-Computeb(NDNDData,BiomassSeries[t,],t)
    
    # Delete irrelevant constraints using possibleAb
    Abp<-possibleAb(A,b,NDNDData$PFv)
    pA<-Abp[[1]];pb<-Abp[[2]]               # Defines two matrices with the existing constraint on flows and biomasses
    
    TRY=try(Fsample<-Sampling(pA,pb))
    TRYAB<-try(pA%*%Fsample<=pb)
    
    if(inherits(TRY, "try-error")==FALSE && length(TRYAB[TRYAB==FALSE])==0){
       # Reattributing the flows values
       Fluxes<-rep(0,NDNDData$nn)             # Creating a vector of 0 and of length nn
       Fluxes[NDNDData$PFv==1]<-Fsample       # Attribute the flow values at the right place according to the vector PFv
       
       BiomassSeries[t+1,]=ComputeBiomass(NDNDData,BiomassSeries[t,],Fluxes,t) # Compute biomass at the next time step
       FlowSeries[t,]<-Fsample # store trophic fluxes
       t=t+1 # Incrementation
    } else {
      if (t>=10) { # if there is no solution, one of the condition is not fullfilled, then it means that there is no polytope or that sampling is outside the polytope, then we stop the sampling and go back 5 time steps.  
        print("crashed - go back 5 time steps")
        t=t-5
        tcrash=tcrash+1 # When this happens, we add 1 to the number of crashs
      }
    }
  }
  
  # Save the number of crashs
  tcrash_m[i]<-tcrash
  
  # reattributes species names for columns and year numbers for rows
  colnames(BiomassSeries)<-NDNDData$Species
  rownames(BiomassSeries)<-1:NDNDData$Tmax
  colnames(FlowSeries)<-NDNDData$flows
  rownames(FlowSeries)<-1:(NDNDData$Tmax-1)

  # Transform the time series in data frames 
  BiomassSeries<-cbind(as.data.frame(BiomassSeries),"Simulation"=rep(i,times=nrow(BiomassSeries)))
  FlowSeries<-cbind(as.data.frame(FlowSeries),"Simulation"=rep(i,times=nrow(FlowSeries)))
  NDNDData$Exportall<-cbind(as.data.frame(NDNDData$Exportall),"Simulation"=rep(i,times=nrow(NDNDData$Exportall)))
  
  # 3. return simulation outputs  -------------------------------------------
  
  NDNDOutput=list(BiomassSeries=BiomassSeries,FlowSeries=FlowSeries)         # Creates a list with all outputs
  
  # Load the code to insert it into the NDNDSimulation file
  NDNDCode=NULL
  initcodefile=paste(NDNDData$directories$code_dir,'/0_1_NDND_initialization.r',sep='')
  datacodefile=paste(NDNDData$directories$code_dir,'/0_2_NDND_data.r',sep='')
  dataprepcodefile=paste(NDNDData$directories$code_dir,'/0_3_NDND_dataprep.r',sep='')
  simcodefile=paste(NDNDData$directories$code_dir,'/1_NDND_simulation.r',sep='')
  if (file.exists(initcodefile)==TRUE && file.exists(datacodefile)==TRUE && file.exists(dataprepcodefile)==TRUE){
    NDNDCode$init <- scan(initcodefile,what="",sep="\n")
    NDNDCode$data <- scan(datacodefile,what="",sep="\n")
    NDNDCode$dataprep <- scan(dataprepcodefile,what="",sep="\n")
    NDNDCode$sim <- scan(simcodefile,what="",sep="\n")
    for (i in 1:length(NDNDfunctions)){
      function2scan<-paste(NDNDData$directories$functions_dir,"/",NDNDfunctions[i],sep="")
      NDNDCode$functions[[i]]=scan(function2scan,what="",sep="\n")
    }
  }
    # Create a list containing the input data, output and the code used to run the simulation  
    NDNDSimulation=list(Simulation.tag=Simulation.tag,
                        Data=NDNDData,
                        Tcrash=tcrash_m,
                        Output=NDNDOutput,
                        Code=NDNDCode)
    
    # Object to return from the function
    return(NDNDSimulation)
  
}