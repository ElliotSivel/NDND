###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model simulation
#### Version v1.1.4
#### 27.03.19
#### Author : Elliot Sivel
###############################################################################################

# 1. Initialization ----------------------------------------------------------

## The first step of the code is to initialize the simulation by clean the  work environment and defining the work directory

rm(list=ls())                                        #clear the work environment

# This code is not generalized. You need to enter your work directory just below.
wd<-getwd()                                              #Set the new work directory

date_time_name<-paste(format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),"_",sep = "")
# 2. load files --------------------------------------------------------------

# To run the entire code, you need a set of 6 files.
## The species.txt file contains the names of the trophospecies, in the present case 8.
## The fluxes.txt file contains the matrix of possible fluxes.
### It is a matrix of dimension number of species by number of species.
### When there is trophic relationship between two species, the value is 1.
### When there no trophic relationship, the value is 0.
## The coefs.txt file contain the values defined for each trophospecies of each coefficient.
### The first coefficient is the initial Biomass
### The second is Gama, the potential assimilation efficiency rate
### The third is Kapa, the food digestability
### The fourth is Mu, the metabolic losses
### The fifth is Rho, the inertia (the maximum change rate of biomasses)
### The sixth is Sigma, the satiation
### The seventh is Bta (Beta), the refuge biomass (threshold under which the prey is not available for the predator)
## The import.txt file is the vector of values of biomass import for each trophospecies at each time step
### Here the matrix is limited to a vector as we assume the import constant over time.
## The export.txt file is the vector of values of biomass export for each trophospecies at each time step
### Here the matrix is null. There is no export.
## The sixth file is the configuration file.
### It contains the names of the files listed above, the number of years for simulation and a plotting parameter.
### If the parameter equals 1, we plot teh time series.
### If the parameter equals 0, we don't plot the time series.

# You have to create a .txt file following the description above to be able to run the rest of the code

setwd("./ndnd_config")
NDNDconfig<-file.choose()                                   # We are opening a window to where you can get your configuration file.
                                                            # You can name as you want, You just have to select it and its directory will automaticaly be loaded in R
setwd(wd)

# We load the data by applying the function
## All data for the rest of the simulation are contained in the NDNDData list
NDNDData<-readDATA(NDNDconfig = NDNDconfig)



# To be able to give the initial condition of each simulation we are running, we are saving the NDNDData file with date and time
setwd("./ndnd_data")
save(NDNDData,file=paste("NDNDData_",date_time_name,".RData"))
setwd(wd)

# 3. Report Model Configuration --------------------------------------------

# This function is not a function that is major for the model
## It is a function that creates a report of the data input for a simulation


# We save the output in a file NDNDConfigreport in a txt format
setwd("./ndnd_configreport")
capture.output(NDNDConfigreport(NDNDData = NDNDData),file = paste("NDND_Config_",date_time_name,".txt",sep = ""),append = T)
setwd(wd)

#####

# The entire model relies on the definition of constraint that are expressed in the form of inequalities
## Ax<b
## A is the set of constraints on the fluxes
## x is a matric of fluxes 
## b is the set of constraints set on biomasses

#####
# 4. Compute A : Matrix of constraint on flows ----------------------------

# This functiona is calculating the constranints on fluxes
## It computes the matrix A that is going to multiplying the flux vector
## It depends on two parameters :
## Gama : the Assimilation efficiency of preys by the predators, and Kapa : the food degistability
## A is computed only once in the simulation and is constant overs the years
ComputeA<-function(Gama,Kapa){

  # As done in the previous function, we define a flag that allow us to identify the simulation
  cat("NDND model ComputeA : ",as.character(Sys.time()),"\n")             # Print Date and Time
  a<-R.Version()
  cat("R Version : ", a$version.string,"\n")                              # Print the R version with which the function was run
  b<-Sys.info()
  cat("Computer Information : ", b[1],b[2],b[3],b[5],"\n")                # Print the Computer software version (Windows, Mac or Linux)
  cat("","\n")
  
  # We set a list of initialization parameter to make the function general and not only applicable to the single set of parameters and input we have now
  ne<-length(Gama)                              # Get the number of species, Gama is vector of length number of species
  nn<-ne^2                                      # Get the number of possible flows
  c<-matrix(1:ne)                               # Create a vector from 1 to 8
  d<-matrix(1,ne,1)                             # Create a vector of ones, dimensions of 1:ne     
  e<-matrix(1,1,nn)                             # Create a vector of ones, dimensions of 1:nn
  
  # We are noe recreating the part of the master equation that correspond to the fluxes
  ## Gama(i)*sum(Fji*Kapa(i))-sum(Fij)
  
  # We create the matrix VI --> the possible going out from one species to the other
  ## It is vectorized (dimension 1:nn) the same way the vector of existing trophic fluxes is. (Mullon's trick for calculation)
  ## It should look like this :
  ## 11111111222222223333333344444444...
  VI<-as.matrix(rep(1:8,byrow=F,each=ne))       # We build a vector from 1 to ne then we replicate the each value ne times
  cat("","\n")
  cat("VI : ","\n")
  cat("","\n")
  print(t(VI))
  cat("","\n")

  cat("------------------------------","\n")

  # We create the matrix VJ --> The possible incoming fluxes 
  ## Here again it is vectorized (Mullon's trick for calculation)
  ## It should look like this :
  ## 12345678123456781234567812345678...
  VJ<-as.matrix(rep(1:8,times=ne))              # We build a vector from 1 to ne and we replicate the vector ne times
  cat("","\n")
  cat("VJ : ","\n")
  cat("","\n")
  print(t(VJ))
  cat("","\n")

  cat("------------------------------","\n")
  
  # We now calculate the matrix SumFij
  ## We are building the matrix of fluxes going out of one prey
  ## The output should look like this :
  ## 11111111000000000000000000000000...
  ## 00000000111111110000000000000000...
  ## 00000000000000001111111100000000...
  ## 00000000000000000000000011111111...
  ## ................................
  ## ................................
  ## ................................
  fij<-c%*%e==d%*%t(VI)                         
  SumFij<-apply(fij,c(1,2),as.numeric)
  cat("","\n")
  cat("SumFij : ","\n")
  cat("","\n")
  print(SumFij)
  cat("","\n")

  cat("------------------------------","\n")

  # We now calculate the matrix SumFji
  ## It is the matrix of possible fluxes coming in from the preys i to a predator j
  ## The output should look like this :
  ## 10000000100000001000000010000000...
  ## 01000000010000000100000001000000...
  ## 00100000001000000010000000100000...
  ## 00010000000100000001000000010000...
  ## ................................
  ## ................................
  ## ................................
  fji<-c%*%e==d%*%t(VJ)
  SumFji<-apply(fji,c(1,2),as.numeric)
  cat("","\n")
  cat("SumFji : ","\n")
  cat("","\n")
  print(SumFji)
  cat("","\n")

  cat("------------------------------","\n")
  
  # We create the Kapa2 matrix
  ## It is a temporary matrix that is created to fit to the dimension of the fluxes matrix
  ## It is used to apply the kapa values to the possible fluxes
  cat("","\n")
  Kapa2<-matrix(rep(as.numeric(Kapa), each = ne), nrow = ne)
  Kapa2<-t(apply(Kapa2,1,rep,each=ne,byrow=TRUE))
  cat("Kapa2 :","\n")
  cat("","\n")
  print(Kapa2)
  cat("","\n")
  
  # We create the matrix SumKiFji
  ## It correspond to the sum of the fluxes weighted by the digestability coefficient of each prey
  cat("SumKjFji :","\n")
  cat("","\n")
  SumKjFji<-SumFji*Kapa2
  print(SumKjFji)
  cat("","\n")
  
  cat("------------------------------","\n")
  cat("------------------------------","\n")
  # Implementing constraints on flows
  ## We determine the constraints by resolving Ax<b inequalities
  ## The term we give here corresponde to the A element of the inequalities
  ## There is five constraints, so 5 A terms.
  
  # First constraint : Biomass are bounded below - Refuge Biomass (Beta)
  cat("A1","\n")
  cat("","\n")
  A1<-SumFij-(Gama%*%e)*SumKjFji
  print(A1)
  cat("","\n")
   
  # Second constraint : Biomass increases are bounded above - Inertia (Rho)
  cat("A2","\n")
  cat("","\n")
  A2=-A1
  print(A2)
  cat("","\n")
   
  # Third constraint : Flows are bounded above - satiation (Sigma)
  cat("A3","\n")
  cat("","\n")
  A3<-SumFji
  print(A3)
  cat("","\n")
  
  # Fourth constraint : Flows are positive
  cat("A4","\n")
  cat("","\n")
  A4<--diag(nn)
  print(A4)
  cat("","\n")

  # Fifth constraint : Biomass decreases are bounded below (-(Rho))
  cat("A5","\n")
  cat("","\n")
  A5=A1
  print(A5)
  cat("","\n")
  
  cat("------------------------------","\n")

  cat("","\n")
  A<-rbind(A1,A2,A3,A4,A5)
  cat("Matrix A : ","\n")
  cat("","\n")
  print(A)
  # print(Gama)
  # print(Gama%*%e)
  # print(SumFij-(Gama%*%e))
  # print((SumFij-(Gama%*%e))*SumKjFji)
  # # A1<-
  return(A)
}

setwd("./ndnd_A")
capture.output(ComputeA(Gama=NDNDData$Gama,Kapa = NDNDData$Kapa),file = paste("NDND_ComputeA_",date_time_name,".txt",sep = ""),append = T)
setwd(wd)

# 5. Compute b : Matrix of constraints on biomasses -----------------------

# Now we are computing the other part of the inequality : b
## It represents the constraints on the biomasses
## Since all biomasses are changing at each time step, the matrix b is always recomputed at each loop
## To run it, 7 parameters are need
## Initial biomass
## The Import data
## The Export data
## Rho, the inertia parameter
## Sigma, the satiation parameter
## Bta (Beta), the refuge biomass
## Mu, the metabolic losses parameter
Computeb<-function(Biomass,Import,Export,Rho,Sigma,Bta,Mu){
  
  # Set the flag
  cat("NDND model Computeb : ",as.character(Sys.time()), t, "\n")
  a<-R.Version()
  cat("R Version : ", a$version.string,"\n")
  b<-Sys.info()
  cat("Computer Information : ", b[1],b[2],b[3],b[5],"\n")
  cat("","\n")

  ne<-length(Rho)                               # Get the number of species
  nn<-ne^2                                      # Get the number of possible flows
  
  C=as.matrix((1-exp(-Mu))/Mu)
  cat("","\n")
  cat("C : ","\n")
  cat("","\n")
  print(C)
  cat("","\n")
  
  cat("------------------------------","\n")
  
  D=as.matrix(exp(-Mu))
  cat("","\n")
  cat("D : ","\n")
  cat("","\n")
  print(D)
  cat("","\n")
  
  cat("------------------------------","\n")
  
  # First constraint : Biomasses are bounded below - Refuge Biomass
  b1<-as.matrix((1/C)*(D*Biomass-Bta)+Import-Export)
  cat("","\n")
  cat("b1 : ","\n")
  cat("","\n")
  print(b1)
  cat("","\n")
  
  cat("------------------------------","\n")
  
  # Second constraint : Biomass increases are bounded above - inertia (Rho)
  b2<-as.matrix(((1/C)*(exp(Rho)-D)*Biomass)+Export)
  cat("","\n")
  cat("b2 : ","\n")
  cat("","\n")
  print(b2)
  cat("","\n")
  
  cat("------------------------------","\n")
  
  # Third constraint : Flows are bounded above - satiation (Sigma)
  b3<-as.matrix(Biomass*Sigma)
  cat("","\n")
  cat("b3 : ","\n")
  cat("","\n")
  print(b3)
  cat("","\n")
  
  cat("------------------------------","\n")
  
  # Fourth constraint, flows are positive
  b4<-matrix(0,nn,1)
  cat("","\n")
  cat("b4 : ","\n")
  cat("","\n")
  print(b4)
  cat("","\n")
  
  cat("------------------------------","\n")
  
  # Fifth constraint, Biomass decreases are bounded below - Inertia (-Rho)
  b5<-as.matrix(((1/C)*(D-exp(-Rho))*Biomass)+Import)
  cat("","\n")
  cat("b5 : ","\n")
  cat("","\n")
  print(b5)
  cat("","\n")
  
  cat("------------------------------","\n")
  cat("------------------------------","\n")
  
  b<-rbind(b1,b2,b3,b4,b5)
  cat("","\n")
  cat("b : ","\n")
  cat("","\n")
  print(b)
  
  return(b)
}

# 6. Delete unpossible flows ----------------------------------------------

# The matrix of possible flows we defined contains number of species^2
## This matrix doesn't represent the topology of the food web
## We need to suppress from the flows matrix the flows that are not possible
### To do so we only keep the flows for the value is PFv is 1. 
possibleAb<-function(A,b,PFv){
  
  # Set the flag
  cat("NDND model possibleAb : ",as.character(Sys.time()),"\n")
  a<-R.Version()
  cat("R Version : ", a$version.string,"\n")
  c<-Sys.info()
  cat("Computer Information : ", c[1],c[2],c[3],c[5],"\n")
  cat("Number of possible flows : ", length(PFv),"\n")
  cat("Number of flows in the given topology : ", sum(PFv),"\n")
  cat("","\n")

  Ap<-A[,which(PFv==1)]                               # We delete the columns in the matrix A for which PFv = 0
  if (ncol(Ap)>1) {                                   # We test if the number of remaining columns is bigger than 1
                                                      # If it is the case then...
    Alines<-apply(abs(Ap), 1,sum)>0                   # We sum up the columns so we only have one vector
                                                      # And we identify the values that are bigger than 0
    Ap<-Ap[which(Alines==TRUE),]                      # We delete the lines that are = 0 
  } else if (ncol(Ap)==1) {                           # If the remaining number of columns equals 1 
    Alines<-Ap>0                                      # Then we test if the values in the remaining vector are = 0 
    Ap<-Ap[which(Alines==TRUE)]                       # We delete the values of the vector that are = 0 
  } else {stop("no topology given")}                  # If there is no remaining columns, it means that no trophic link has been implemented, you have no topology
                                                      # The simulation stops and an error message appears
  
  cat("Ap : ","\n")
  cat("","\n")
  print(Ap)
  cat("","\n")

  cat("------------------------------","\n")
  
  cat("bp : ","\n")
  cat("","\n")
  bp<-as.matrix(b[which(Alines==T)])                  # We delete all elements in the matrix that have been deleted in the matrix A
  print(bp)
  
  return(list(Ap,bp))
}

# 7. Compute Biomass : Calculate the biomass for t+1 ----------------------

# Now that we have defined how to compute A and b of the inequality and that x is going to be draw randomly, we are defining the function that is computing the biomass at time step t+1
## We apply the master equation : B(t+1)=D*B(t)+C*[Gama(i)*sum(Fji*Kapa(i))+Import(i)-sum(Fij)-Export(i)]
### Here we only calculated the biomass for the next time step, not for more
ComputeBiomass<-function(Biomass,F,Import,Export,Gama,Mu,Kapa){
  C=as.matrix((1-exp(-Mu))/Mu)                          # Defining C
  D=as.matrix(exp(-Mu))                                 # Defining D
  ne<-length(Gama)                                      # Defining the number of species
  dim(F)<-c(ne,ne)                                      # F is going to drawn as a vector. Here we resahpe it in a matrix of dimension ne*ne
  FM<-t(F)                                              # and transpose it
  NewBiomass<-matrix(0,ne)                              # We defin a new empty object of ne rows
  for (i in 1:ne){                                      # We fill up the new object with the new estimated biomasses
    NewBiomass[i]=D[i]*Biomass[i]+C[i]*(Gama[i]*sum(FM[,i]*Kapa)-sum(FM[i,])+Import[i]-Export[i])
  }
  return(NewBiomass)
}

# 8. Initialization of the simulation -------------------------------------

# Now that we have defined every functioned needed to run the model, we are defining the object and elements that are needed to run the model
BiomassSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = NDNDData$ns)               # Creates a matrix of dimension Tmax vs number of species full of 0. It is thought be filled with the biomass obtained during the calculations
FlowSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = sum(NDNDData$PFv))            # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
CurrentBiomass=NDNDData$Biomass                                                         # Set the initial biomass as the vector of biomass given in the NDNDData file
BiomassSeries[1,]<-CurrentBiomass                                                       # We give that the biomass for t = 1 is the initial biomass
A<-ComputeA(NDNDData$Gama,NDNDData$Kapa)                                                # We compute A (which constant for the entire simulation)

setwd("./ndnd_A/Rfiles")
save(A,file=paste("NDND_ A_",date_time_name,".RData"))
setwd(wd)

## Not implemented yet, needs to be done for further investigation
# T=1                                                                                     # Set initial time to 1 
# Tcrash=0                                                                                # Set Tcrash to 0

# 9. Main loop -----------------------------------------------------------------

# Here are runing the model over Tmax steps
# We sample the fluxes under a uniform probabilistic law
# Several sample algorithms are existing
## After performing tests with few of them, we have decided to perform the sampling with the Sampling algorithm described in the LIM package
## It is a package developed by Dick van Oevelen and Karline Soetaert (2015)
library(LIM)                                           # loading the library LIM

setwd("./ndnd_b")
dir.create(paste(date_time_name,sep = ""))
setwd(paste(getwd(),"/",date_time_name,sep=""))
for (t in 1:NDNDData$Tmax) {
  print(t)
  # We have defined Importall and Exportall as the matrix for all the years that we are simulating
  # Here we take out one year after another the value of Import and Export to estimate b each year
  Import<-NDNDData$Importall[t,]
  Export<-NDNDData$Exportall[t,]
  # We compute b
  b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu)
  capture.output(Computeb(NDNDData$Biomass,NDNDData$Import,NDNDData$Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu),file = paste("NDND_Computeb",date_time_name,"_",t,".txt",sep = ""),append = T)
  # We delimit the flow matrix to the possible flows by applying the possibleAb function
  Abp<-possibleAb(A,b,NDNDData$PFv)
  pA<-Abp[[1]];pb<-Abp[[2]]
  G<-as.matrix(-pA)                       # Defining A (constraints on fluxes) for the sampling
  h<-as.matrix(-pb)                       # Defining b (constraints on biomasses) for th sampling
  
  # Sampling Fluxes
  Fsample<-xsample(G=G,H=h,iter=100,burninlength=100,outputlength=100,type="mirror")
  # The xsample function returns a list of four elements of which the first on is the matrix of fluxes that has been sampled
  
  # We take the matrix of sampled fluxes in the sampling output
  # The output has the dimension iter*number of fluxes
  ## We sample randomly one vector of fluxes in the matrix
  F0<-Fsample[[1]][sample(1:nrow(Fsample[[1]]),1),]
  F<-rep(0,64)
  F[NDNDData$PFv==1]<-F0
  FlowSeries[t,]<-F0
  BiomassSeries[t+1,]=ComputeBiomass(Biomass = BiomassSeries[t,],F=F,Import = Import,Export = Export,Gama = NDNDData$Gama,Mu = NDNDData$Mu,Kapa = NDNDData$Kapa)
  for (c in 1:ncol(BiomassSeries)) {
    if (BiomassSeries[t+1,c]<0){
      BiomassSeries[t+1,c]<-0
    }
  }
}


# 10. Plots ---------------------------------------------------------------

library(ggplot2)
library(gridExtra)
ggplot(as.data.frame(BiomassSeries[,1]))+aes(1:NDNDData$Tmax,BiomassSeries[,1])+geom_line()

grid.arrange(ggplot(as.data.frame(BiomassSeries[,1]))+aes(1:NDNDData$Tmax,BiomassSeries[,1])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[1]),
             ggplot(as.data.frame(BiomassSeries[,2]))+aes(1:NDNDData$Tmax,BiomassSeries[,2])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[2]),
             ggplot(as.data.frame(BiomassSeries[,3]))+aes(1:NDNDData$Tmax,BiomassSeries[,3])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[3]),
             ggplot(as.data.frame(BiomassSeries[,4]))+aes(1:NDNDData$Tmax,BiomassSeries[,4])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[4]),
             ggplot(as.data.frame(BiomassSeries[,5]))+aes(1:NDNDData$Tmax,BiomassSeries[,5])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[5]),
             ggplot(as.data.frame(BiomassSeries[,6]))+aes(1:NDNDData$Tmax,BiomassSeries[,6])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[6]),
             ggplot(as.data.frame(BiomassSeries[,7]))+aes(1:NDNDData$Tmax,BiomassSeries[,7])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[7]),
             ggplot(as.data.frame(BiomassSeries[,8]))+aes(1:NDNDData$Tmax,BiomassSeries[,8])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[8]),
             ncol=3,nrow=3)

# F<-rep(0,64)
# F[NDNDData$PFv==1]<-f1
# 
# test<-ComputeBiomass(BiomassSeries[1,],F,as.numeric(Import),as.numeric(Export),NDNDData$Gama,NDNDData$Mu,NDNDData$Kapa)

capture.output(Computeb(NDNDData$Biomass,NDNDData$Import,NDNDData$Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu),file = paste("NDND_Computeb",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".txt",sep = ""),append = T)
save(b,file=paste("NDND_ b_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
capture.output(possibleAb(A,b,NDNDData$PFv),file = paste("NDND_possibleAb",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".txt",sep = ""),append = T)
save(Ap,file=paste("NDND_pA_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
save(bp,file=paste("NDND_pb_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
