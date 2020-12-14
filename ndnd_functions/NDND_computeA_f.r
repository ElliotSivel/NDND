###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### ComputeA function
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel
###############################################################################################

##### ComputeA function
# Principle of the model is to determine trophic fluxes randomly in a set of possibilities delimited by constraints
# Constraints expressed in terms of inequalities
# Ax =< b
# ComputeA compute the matrix A : matrix of constraints on flows
# Need 2 parameter vectors : Gama (assimilation efficiency parameter) and Kapa (digestability parameter)
# Need 2 elements : number of species (ns) and number of potential flows (nn)

ComputeA<-function(NDNDData){
  # extract necessary coefficients from the NDNDdata object
  Gama <- NDNDData$Gama
  Kapa <- NDNDData$Kapa
  ns <- NDNDData$ns
  nn <- NDNDData$nn
  
  # We build three temporary matrix that are going to be used to build the final matrices
  
  c<-matrix(1:ns)                               # Create a vector from 1 to ns
  d<-matrix(1,ns,1)                             # Create a vector of ones, dimensions of 1:ns     
  e<-matrix(1,1,nn)                             # Create a vector of ones, dimensions of 1:nn
  
  # Computing the part of the master equation that correspond to the fluxes
  # Gama(i)*sum(Fji*Kapa(i))-sum(Fij)
  
  # Vector VI --> the possible flows towards predators for each species
  VI<-as.matrix(rep(1:ns,byrow=F,each=ns))      # Vector from 1 to ns
                                                # Replication each value ns times
                                                # Vector dimension 1:nn
  
  # Vector VJ --> The possible flows incoming for each species
  VJ<-as.matrix(rep(1:ns,times=ns))             # Vector from 1 to ns
                                                # Replication of the vector ns times
                                                # Vector of dimension 1:nn
  
  # Compute Fij matrix --> SumFij
  # Matrix of flows outgoing from the species i to the species j
  Fij<-c%*%e==d%*%t(VI)                         # Matrix product of c-e and d-VI
                                                # Fij is the matrix were the two computed matricial product are equal
  SumFij<-apply(Fij,c(1,2),as.numeric)          # Transformation from bolean to numerical values
  
  # Compute Fji matrix --> SumFji
  # Matrix of flows incoming species i from species j 
  Fji<-c%*%e==d%*%t(VJ)                         # Matrix product of c-e and d-VJ
                                                # Fji is the matrix were the two computed matricial product are equal
  SumFji<-apply(Fji,c(1,2),as.numeric)          # Transformation from bolean to numerical values
 
  # Compute Kapa matrix
  # For all parameters, values are in vectors
  # for matricial computation, we need to replicate the Kapa values --> Kapa2 
  Kapa2<-matrix(rep(as.numeric(Kapa), each = ns), nrow = ns)             # Replication of the kapa vector in rows ns times
  Kapa2<-t(apply(Kapa2,1,rep,each=ns,byrow=TRUE))                        # Replication of the columns of Kapa2 ns times

  # Compute SumKiFji matrix
  # (Vectorized) Matrix of flows from the species j to the species i -- predation flows 
  SumKjFji<-SumFji*Kapa2                          # Multiplying each flow of the SumFji matrix to the corresponding value of Kapa
 
  # Implementing constraints on flows
  # There are 5 constraints -- A1,A2,A3,A4,A5
  
  # First constraint : Biomass are bounded below - Refuge Biomass (Beta)
  A1<-SumFij-(Gama%*%e)*SumKjFji
  
  # Second constraint : Biomass increases are bounded above - Inertia (Rho)
  A2=-A1
  
  # Third constraint : Flows are bounded above - satiation (Sigma)
  A3<-SumFji
 
  # Fourth constraint : Flows are positive
  A4<--diag(nn)
 
  # Fifth constraint : Biomass decreases are bounded below (-(Rho))
  A5=A1
  
  # Create an object where all constraints matrix are just bound --> A
  A<-rbind(A1,A2,A3,A4,A5)
  
  return(A)
}