###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Version v1.0
#### 10.12.18
#### Author : Elliot Sivel
###############################################################################################


# 1. Initialization ----------------------------------------------------------

rm(list=ls())                                        #clear the work environment                                                      
getwd()                                              #Search for the work environment
setwd("C:/Users/a22073/Documents/Model_NDND")        #Set the new work directory
wd<-getwd()                                          #Save the work directory

# 2. load files --------------------------------------------------------------

# To load the data, we define a function for reading the data. Aim of the readDATA function is to provide an object containing all data and parameters needed to run the model.

NDNDconfig<-file.choose()             # We get the file directory of the file containing the names of the files, the number of runs and the parameters defining if we plot or not at the end

readDATA<-function(NDNDconfig){             # The variable NDNDconfig correspond to the file directory of the file NDNDconfig containing the file names
  
  file.names<-scan(NDNDconfig,what=character())             # Get the file names as characters : species, coefficients, flux, import, export. The file also contain the number of runs and the plot option parameter
  
  # Verifying if the files containg the data are existing
  exist.file<-vector(length=5)            # Create an empty vector. The length of the vector corresponds to the number of file needed
  # For each file testing if it exist in the work directory
  for (i in (1:(length(file.names)-2))) {
    exist.file[i]<-file.exists(file.names[i])
  }
  
  #####
  
  # If all file are available in the work directory, then load all files
  if (sum(exist.file)==5){
    
      Species<-scan(file=file.names[1],what=character(),quote="")             # Read the species names
      PF<-read.table(file.names[2], header=F, sep="\t",stringsAsFactors = F)            # Load the matrix of possible flows
      Coefficients<-as.matrix(read.table(file.names[3], header=F, sep="\t",stringsAsFactors = F))             # Load the values of coefficients 
      Import<-scan(file.names[4], what=numeric())             # Read the import values
      Export<-scan(file.names[5], what=numeric())             # Read the export values
  }  else {stop("Missing file")}         # If one file is missing, then stop the calculations and return "missing file"
  
  #####
  
  # Once all files have been loaded, we extract and format the data in the shape we want it for running the model
  
  # We extract the number of runs that we are running from the NDNDConfig file
  Tmax<-as.numeric(file.names[6])
  
  # We extract the names of species
  Speciesnames<-Species
  
  # We estimate the number of possible flowes. It corresponds to nb of species*nb of species
  ns<-as.numeric(length(Species))             # Number of species
  nn<-ns*ns             # Number of possible flows in the food web
  
  # We set the column names and the row names of the coefficient matrix
  colnames(Coefficients)<-Species
  rownames(Coefficients)<-c("Biomass","Gama","Kapa","Mu","Rho","Sigma","Bta")
  # We export then the data in single vectors or matrix
  Biomass<-Coefficients[1,]               # Row 1 is the initial biomass for each species, it is given in tons per km2
  Gama<-Coefficients[2,]                  # Row 2 is Gama, the potential assimiliation efficiency for each species
  Kapa<-Coefficients[3,]                  # Row 3 is Kapa, the digestability of each species
  Mu<-Coefficients[4,]                    # Row 4 is Mu, the metabolic losses for each species
  Rho<-Coefficients[5,]                   # Row 5 is Rho, the parameter accounting for inertia
  Sgma<-Coefficients[6,]                  # Row 6 is Sigma, the parameter accounting for satiation. It has been renamed Sgma because Sigma is already a function in R
  Bta<-Coefficients[7,]                   # Row 7 is Beta, the value for refuge biomass
  
  # Import and Export matrix must have Tmax*Nb of species for dimensions
  # If the import/export file has only onr row, the values in the row are copied for all time step
  if (nrow(as.matrix(t(Import)))==1){
    Importall<-matrix(rep(as.matrix(t(Import)),each=Tmax),nrow=Tmax)                        # If the nb of rows in Import equals 1, then we duplicate the vector times Tmax
  } else if (nrow(as.matrix(t(Import)))==Tmax){
    Importall<-Import                                                                       # If the nb of rows in Import equals Tmax, then the matrix Import equals Importall
  } else if (nrow(as.matrix(t(Import))) > Tmax){                      
    Importall<-Import[1:Tmax,]                                                              # If the nb of rows in Import is greater than Tmax, then the matrix Import from row 1 to Tmax equals Importall
  } else if (nrow(as.matrix(t(Import))) < Tmax && nrow(as.matrix(t(Import))) != 1){
    stop("1 < number of rows < Tmax")                                                       # If the nb of rows in Import is between one and Tmax, then the calculation stops
  } else {stop("Import data missing")}                                                      # If none of the condition above is fulfilled, then the nb of rows is equals to 0, then the calculation stops because import data is missing
  
  
  if (nrow(as.matrix(t(Export)))==1){
    Exportall<-matrix(rep(as.matrix(t(Export)),each=Tmax),nrow=Tmax)                        # If the nb of rows in Export equals 1, then we duplicate the vector times Tmax
  } else if (nrow(as.matrix(t(Export)))==1){
    Exportall<-Export                                                                       # If the nb of rows in Export equals Tmax, then the matrix Export equals Exportall
  } else if (nrow(as.matrix(t(Export))) > Tmax){
    Exportall<-Export[1:Tmax,]                                                              # If the nb of rows in Export is greater than Tmax, then the matrix Export from row 1 to Tmax equals Exportall
  } else if (nrow(as.matrix(t(Export))) < Tmax && nrow(as.matrix(t(Export))) != 1){
    stop("1 < number of rows < Tmax")                                                       # If the nb of rows in Export is between one and Tmax, then the calculation stops
  } else {stop("Export data missing")}                                                      # If none of the condition above is fulfilled, then the nb of rows is equals to 0, then the calculation stops because export data is missing
  
  # Then we vectorize the matrix of possible flows : trick used by Christian mullon to facilitate the calculation of the model
  PFv=as.vector(t(PF))
  
  # Finally we set the plotting parameter which a bolean value for plotting or not the graphs
  Plotting<- file.names[7]                          # If Plotting equals 1, we plot the graphs. If it equals 0 we don't               
  
  return(list(Tmax,Speciesnames,ns,nn,Biomass,Gama,Kapa,Mu,Rho,Sgma,Bta,Importall,Import,Exportall,Export,PF,PFv,Plotting))
}

NDNDData<-readDATA(NDNDconfig = NDNDconfig)

save(NDNDData,file=paste("NDNDData_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))

# We add the species names to the tables where they are needed
names(NDNDData)<-c("Tmax","Speciesnames","ns","nn","Biomass","Gama","Kapa","Mu","Rho","Sigma","Beta","Importall","Import","Exportall","Export","PF","PFv","Plotting")
colnames(NDNDData$PF)<-NDNDData$Speciesnames
rownames(NDNDData$PF)<-NDNDData$Speciesnames
colnames(NDNDData$Importall)<-NDNDData$Speciesnames
colnames(NDNDData$Exportall)<-NDNDData$Speciesnames

# 3. Report Model Configuration --------------------------------------------

NDNDConfigreport<-function(NDNDData){
  
  #####
  # Display data implemented in the model
  
  cat("NDND model configuration : ",as.character(Sys.time()),"\n")
  a<-R.Version()
  cat("R Version : ", a$version.string,"\n")
  b<-Sys.info()
  cat("Computer Information : ", b[1],b[2],b[3],b[5],"\n")
  cat("Number of Species : ",NDNDData$ns,"\n")
  
  cat("------------------------------","\n")
  
  cat("Species names : ","\n")
  cat("","\n")
  print(NDNDData$Speciesnames)
  
  cat("------------------------------","\n")
  
  cat("Starting Biomasses : ","\n")
  cat("","\n")
  print(NDNDData$Biomass)
  
  cat("------------------------------","\n")
  
  cat("Potential assimilation efficienties (gamma) :","\n")
  cat("","\n")
  print(NDNDData$Gama)
  
  cat("------------------------------","\n")
  
  cat("Food quality :","\n")
  cat("","\n")
  print(NDNDData$Kapa)
  
  cat("------------------------------","\n")
  
  cat("Other intial metabolic losses :","\n")
  cat("","\n")
  print(NDNDData$Mu)
  
  cat("------------------------------","\n")
  
  cat("Intertia :","\n")
  cat("","\n")
  print(NDNDData$Rho)
  
  cat("------------------------------","\n")
  
  cat("Maximum biomass growth rate :","\n")
  cat("","\n")
  print(exp(NDNDData$Rho))
  
  cat("------------------------------","\n")
  
  cat("Maximum biomass decline rate :","\n")
  cat("","\n")
  print(exp(-(NDNDData$Rho)))
  
  cat("------------------------------","\n")
  
  cat("Satiation :","\n")
  cat("","\n")
  print(NDNDData$Sigma)
  
  cat("------------------------------","\n")
  
  cat("Food requirements with stables biomass (no predation) : ","\n")
  cat("","\n")
  print(NDNDData$Mu/NDNDData$Gama)
  
  cat("------------------------------","\n")
  
  cat("Food requirements in year 2 with stable biomass (no predation) : ","\n")
  cat("","\n")
  print((NDNDData$Mu/NDNDData$Gama)*NDNDData$Biomass)
  
  cat("------------------------------","\n")
  
  cat("Food requirements with max decrease in biomass (no predation) : ","\n")
  cat("","\n")
  print((NDNDData$Mu/NDNDData$Gama)*(exp(-NDNDData$Rho)-exp(-NDNDData$Mu))/(1-exp(-NDNDData$Mu)))
  
  cat("------------------------------","\n")
  
  cat("","\n")
  # Test to check if the satiation has been set high enough for the calculations.
  # In case it is too low for one or more species, then a warning message should be printed.
  # This test does not work for the phytoplankton which does not have a value for satiation since it does not predate.
  satiation.setup<-NDNDData$Mu/NDNDData$Gama>NDNDData$Sigma
  F<-satiation.setup[which(satiation.setup==T)]
  if (length(F) > 0)  {
    cat("Warning : Satiation is set too low for the following species","\n")
    print(names(F))
  } else { cat("All satiation values high enough for calculation","\n")}
  cat("","\n")
    
  cat("------------------------------","\n")
  
  cat("","\n")
  cat("Refuge Biomass (Beta) : ","\n")
  cat("","\n")
  print(NDNDData$Beta)
  
  cat("------------------------------","\n")
  
  cat("","\n")
  cat("Biomass import in the system : ","\n")
  cat("","\n")
  print(head(NDNDData$Importall))
  cat("","\n")
  cat("Mean Import : ","\n")
  print(apply(NDNDData$Importall,2,mean))
  cat("","\n")
  cat("Variance of import : ","\n")
  print(apply(NDNDData$Importall,2,var))
  
  cat("------------------------------","\n")
  
  cat("","\n")
  cat("Biomass export out of the system : ","\n")
  cat("","\n")
  print(head(NDNDData$Exportall))
  cat("","\n")
  cat("Mean Export : ","\n")
  print(apply(NDNDData$Exportall,2,mean))
  cat("","\n")
  cat("Variance of export : ","\n")
  print(apply(NDNDData$Exportall,2,var))
  
  cat("------------------------------","\n")
  
  cat("","\n")
  cat("Who eats whom? : Trophic flows matrix (Preys in rows / Predators in columns)","\n")
  cat("","\n")
  print(NDNDData$Speciesnames)
  cat("","\n")
  print(NDNDData$PF)
  
}

# capture.output(NDNDConfigreport(NDNDData = NDNDData),file = paste("NDND_Config_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".txt",sep = ""),append = T)     # Save a file containing all parameters in the work directory

# 4. Compute A : Matrix of constraint on flows ----------------------------

ComputeA<-function(Gama,Kapa){

  #####
  # Set the flag
  # cat("NDND model ComputeA : ",as.character(Sys.time()),"\n")
  # a<-R.Version()
  # cat("R Version : ", a$version.string,"\n")
  # b<-Sys.info()
  # cat("Computer Information : ", b[1],b[2],b[3],b[5],"\n")
  # cat("","\n")
  
  ne<-length(Gama)                              # Get the number of species
  nn<-ne^2                                      # Get the number of possible flows
  c<-matrix(1:ne)                                # Create a vector from 1 to 8
  d<-matrix(1,ne,1)                              # Create a vector of ones, dimensions of 1:ne     
  e<-matrix(1,1,nn)                             # Create a vector of ones, dimensions of 1:nn
  
  # Creating matrix VI
  VI<-as.matrix(rep(1:8,byrow=F,each=ne))
  # cat("","\n")
  # cat("VI : ","\n")
  # cat("","\n")
  # print(t(VI))
  # cat("","\n")
  # 
  # cat("------------------------------","\n")
  # 
  # # Creating matrix VJ
  VJ<-as.matrix(rep(1:8,times=ne))
  # cat("","\n")
  # cat("VJ : ","\n")
  # cat("","\n")
  # print(t(VJ))
  # cat("","\n")
  
  # cat("------------------------------","\n")
  # 
  fij<-c%*%e==d%*%t(VI)
  SumFij<-apply(fij,c(1,2),as.numeric)
  # cat("","\n")
  # cat("SumFij : ","\n")
  # cat("","\n")
  # print(SumFij)
  # cat("","\n")
  # 
  # cat("------------------------------","\n")
  # 
  fji<-c%*%e==d%*%t(VJ)
  SumFji<-apply(fji,c(1,2),as.numeric)
  # cat("","\n")
  # cat("SumFji : ","\n")
  # cat("","\n")
  # print(SumFji)
  # cat("","\n")
  # 
  # cat("------------------------------","\n")
  # 
  # cat("","\n")
  Kapa2<-matrix(rep(as.numeric(Kapa), each = ne), nrow = ne)
  Kapa2<-t(apply(Kapa2,1,rep,each=ne,byrow=TRUE))
  # cat("Kapa2 :","\n")
  # cat("","\n")
  # print(Kapa2)
  # cat("","\n")
  # 
  # cat("SumKjFji :","\n")
  # cat("","\n")
  SumKjFji<-SumFji*Kapa2
  # print(SumKjFji)
  # cat("","\n")
  # 
  #####
  # Implementing constraints on flows
  
  # cat("------------------------------","\n")
  
  # First constraint : Biomass are bounded below - Refuge Biomass (Beta)
  # cat("A1","\n")
  # cat("","\n")
  A1<-SumFij-(Gama%*%e)*SumKjFji
  # print(A1)
  # cat("","\n")
  # 
  # Second constraint : Biomass increases are bounded above - Inertia (Rho)
  # cat("A2","\n")
  # cat("","\n")
  A2<--A1
  # print(A2)
  # cat("","\n")
  # 
  # Third constraint : Flows are bounded above - satiation (Sigma)
  # cat("A3","\n")
  # cat("","\n")
  A3<-SumFji
  # print(A3)
  # cat("","\n")
  # 
  # Fourth constraint : Flows are positive
  # cat("A4","\n")
  # cat("","\n")
  A4<--diag(nn)
  # print(A4)
  # cat("","\n")
  # 
  # Fifth constraint : Biomass decreases are bounded below (-(Rho))
  # cat("A5","\n")
  # cat("","\n")
  A5=A1
  # print(A5)
  # cat("","\n")
  
  # cat("------------------------------","\n")
  # 
  # cat("","\n")
  A<-rbind(A1,A2,A3,A4,A5)
  # cat("Matrix A : ","\n")
  # cat("","\n")
  # print(A)
  # print(Gama)
  # print(Gama%*%e)
  # print(SumFij-(Gama%*%e))
  # print((SumFij-(Gama%*%e))*SumKjFji)
  # # A1<-
  return(A)
}

# 5. Compute b : Matrix of constraints on biomasses -----------------------

Computeb<-function(Biomass,Import,Export,Rho,Sigma,Bta,Mu){
  
  #####
  # Set the flag
  # cat("NDND model Computeb : ",as.character(Sys.time()),"\n")
  # a<-R.Version()
  # cat("R Version : ", a$version.string,"\n")
  # b<-Sys.info()
  # cat("Computer Information : ", b[1],b[2],b[3],b[5],"\n")
  # cat("","\n")
  # 
  ne<-length(Rho)                              # Get the number of species
  nn<-ne^2                                      # Get the number of possible flows
  
  #####
  C=as.matrix((1-exp(-Mu))/Mu)
  # print(C)
  D=as.matrix(exp(-Mu))
  # print(D)
  
  # First constraint : Biomasses are bounded below - Refuge Biomass
  b1<-as.matrix((1/C)*(D*Biomass-Bta)+Import-Export)
  # print(b1)
  
  # Second constraint : Biomass increases are bounded above - inertia (Rho)
  b2<-as.matrix((1/C)*(exp(Rho)-D)*Biomass+Import)

  # Third constraint : Flows are bounded above - satiation (Sigma)
  b3<-as.matrix(Biomass*Sigma)
  # print(b3)
  
  # Fourth constraint, flows are positive
  b4<-matrix(0,nn,1)
  # print(b4)
  
  # Fifth constraint, Biomass decreases are bounded below - Inertia (-Rho)
  b5<-as.matrix((1/C)*(D-exp(-Rho))*Biomass+Export)
  # print(b5)
  
  b<-rbind(b1,b2,b3,b4,b5)
  return(b)
}

# 6. Delete unpossible flows ----------------------------------------------

possibleAb<-function(A,b,PFv){
  
  # Set the flag
  # cat("NDND model possibleAb : ",as.character(Sys.time()),"\n")
  # a<-R.Version()
  # cat("R Version : ", a$version.string,"\n")
  # c<-Sys.info()
  # cat("Computer Information : ", c[1],c[2],c[3],c[5],"\n")
  # cat("","\n")
  # 
  Ap<-A[,which(PFv==1)]
  if (ncol(Ap>1)) {
    Alines<-apply(abs(Ap), 1,sum)>0
    Ap<-Ap[which(Alines==TRUE),]
  } else {
    Alines<-Ap>0
    Ap<-Ap[which(Alines==TRUE)]
  }
  
  # cat("Ap : ","\n")
  # cat("","\n")
  # print(Ap)
  # cat("","\n")
  # 
  # cat("bp : ","\n")
  # cat("","\n")
  bp<-as.matrix(b[which(Alines==T)])
  # print(bp)
  return(list(Ap,bp))
}

# 7. Initialization of the simulation -------------------------------------

BiomassSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = NDNDData$ns)               # Creates a matrix of dimension Tmax vs number of species full of 0. It is thought be filled with the biomass obtained during the calculations
FlowSeries<-matrix(data = 0, nrow = NDNDData$Tmax, ncol = sum(NDNDData$PFv))            # Creates a matrix of dimension Tmax vs the number of possible flows. We kept only the links for which we had a 1 in PF. There are 18 links in our topology.
T=1                                                                                     # Set initial time to 1 
Tcrash=0                                                                                # Set Tcrash to 0
CurrentBiomass=NDNDData$Biomass                                                         # Set the initial biomass as the vector of biomass given in the NDNDData file

# capture.output(ComputeA(NDNDData$Gama,NDNDData$Kapa),file = paste("NDND_ComputeA",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".txt",sep = ""),append = T)
A<-ComputeA(NDNDData$Gama,NDNDData$Kapa)
# save(A,file=paste("NDND_ A_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))

BiomassSeries[1,]<-CurrentBiomass

# 8. Compute Biomass : Calculate the biomass for t+1 ----------------------

ComputeBiomass<-function(Biomass,F,Import,Export,Gama,Mu,Kapa){
  C=as.matrix((1-exp(-Mu))/Mu)
  D=as.matrix(exp(-Mu))
  ne<-length(Gama)
  dim(F)<-c(ne,ne)
  FM<-t(F)
  NewBiomass<-matrix(0,ne)
  for (i in 1:ne){
    NewBiomass[i]=D[i]*Biomass[i]+C[i]*(Gama[i]*sum(FM[,i]*Kapa)-sum(FM[i,])+Import[i]-Export[i])
  }
  return(NewBiomass)
}

# 9. Main loop -----------------------------------------------------------------

library(LIM)

for (t in 1:NDNDData$Tmax) {
  while (t<NDNDData$Tmax) {
    print(t)
    Import<-NDNDData$Importall[t,]
    Export<-NDNDData$Exportall[t,]
    b<-Computeb(BiomassSeries[t,],Import,Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu)
    Abp<-possibleAb(A,b,NDNDData$PFv)
    pA<-Abp[[1]];pb<-Abp[[2]]
    G<-as.matrix(-pA)
    h<-as.matrix(-pb)
    Fsample<-xsample(G=G,H=h,iter=100,burninlength=100,outputlength=100,type="mirror")
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
}

library(ggplot2)
ggplot(as.data.frame(BiomassSeries[,8]))+aes(1:NDNDData$Tmax,BiomassSeries[,8])+geom_line()

# F<-rep(0,64)
# F[NDNDData$PFv==1]<-f1
# 
# test<-ComputeBiomass(BiomassSeries[1,],F,as.numeric(Import),as.numeric(Export),NDNDData$Gama,NDNDData$Mu,NDNDData$Kapa)

capture.output(Computeb(NDNDData$Biomass,NDNDData$Import,NDNDData$Export,NDNDData$Rho,NDNDData$Sigma,NDNDData$Beta,NDNDData$Mu),file = paste("NDND_Computeb",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".txt",sep = ""),append = T)
save(b,file=paste("NDND_ b_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
capture.output(possibleAb(A,b,NDNDData$PFv),file = paste("NDND_possibleAb",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".txt",sep = ""),append = T)
save(Ap,file=paste("NDND_pA_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
save(bp,file=paste("NDND_pb_",format(Sys.time(),"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
