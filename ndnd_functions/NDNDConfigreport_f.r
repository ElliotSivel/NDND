###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### NDNDConfigreport function
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel
###############################################################################################

##### NDNDConfigreport function
# This function is not having any role in the model
# It is providing a .txt file reporting alls the setups of the model
# The list of species, the time of simulation, the parameter values and some other information about the system
# To run it needs the NDNDData file created by the function readDATA

NDNDConfigreport<-function(directory=".",NDNDData){
  
  sink(paste(directory,"/NDND_Config_",paste(format(NDNDData$Data.tag,"%Y_%m_%d_%H_%M_%S"),sep = ""),".txt",sep=""))        # Opens a file with the defined name, in the directory indicated in the simulation code
                                                          # Print all the wanted elements of function without printing them in the console
                                                          # File accessible and created after sink()
  
  # Loading the information about the running session -- systems information, packages, R version, OS,...
  
  session_info<-sessionInfo()                     
  
  # Extract the version of loaded packages
  
  Package_version<-packages.version.NDND(Session_info = session_info)
  
  # Set identification flag
  
  cat("NDND model configuration : ",as.character(NDNDData$Data.tag),"\n")                # Prints out the date and time of running
  a<-R.Version()                                                                  
  cat("R Version : ", a$version.string,"\n")                                      # Prints out the version of R used to run the model
  b<-Sys.info()                                                                     
  cat("Computer Information : ", b[1],b[2],b[3],b[5],"\n")                        # Prints out the software system (Windows, Mac or Linux) version
  cat("Package information : ","\n")
  for (i in 1:dim(Package_version)[1]){
    cat(Package_version[i,1]," : ",Package_version[i,2],"\n")
    }
  cat("Number of Species : ",NDNDData$ns,"\n")                                    # Prints out the number of species
  cat("Length of the simulation : ", NDNDData$Tmax," years","\n")                 # Prints out the length of the simulation
  cat("Sample algorithm : ", NDNDData$Sampling_algorithm, "\n")                   # Prints out the chosen sampling algorithm
  cat("Plotting parameter :", NDNDData$Plotting,"\n")                             # Prints out if the time series are going to be plotted
  cat("------------------------------","\n")
  
  # Prints out the species names implemented in the model
  
  cat("Species names : ","\n")
  cat("","\n")
  print(NDNDData$Species)                                                         # Prints out the species vector in NDNDData
  cat("------------------------------","\n")
  
  # Prints out the initial biomasses implemented in the model
  
  cat("Starting Biomasses : ","\n")
  cat("","\n")
  print(NDNDData$Biomass)                                                         # Prints out the vector of initial biomass
  cat("------------------------------","\n")
  
  # Prints out the values for the Gama parameter implemented in the model
  
  cat("Potential assimilation efficienties (gamma) :","\n")
  cat("","\n")
  print(NDNDData$Gama)                                                            # Prints out the vector of Gama parameter values 
  cat("------------------------------","\n")
  
  # Prints out the values for the Kapa parameter implemented in the model
  
  cat("Food quality :","\n")
  cat("","\n")
  print(NDNDData$Kapa)                                                            # Prints out the vector of Kapa parameter values 
  cat("------------------------------","\n")
  
  # Prints out the values for the Mu parameter implemented in the model
  
  cat("Other intial metabolic losses :","\n")
  cat("","\n")
  print(NDNDData$Mu)                                                              # Prints out the vector of Mu parameter values 
  cat("------------------------------","\n")
  
  # Prints out the values for the Rho parameter implemented in the model
  
  cat("Inertia :","\n")
  cat("","\n")
  print(NDNDData$Rho)                                                             # Prints out the vector of Rho parameter values 
  cat("------------------------------","\n")
  
  # We calculate and print out the maximum growth rate as the exponential value of Rho
  
  cat("Maximum biomass growth rate :","\n")
  cat("","\n")
  print(exp(NDNDData$Rho))                                                        # Computes and prints out the vector of exp(Rho) parameter values 
  cat("------------------------------","\n")
  
  # We calculate and print out the maximum decline rate as the exponential value of -Rho
  
  cat("Maximum biomass decline rate :","\n")
  cat("","\n")
  print(exp(-(NDNDData$Rho)))                                                     # Computes and prints out the vector of exp(-Rho) parameter values
  cat("------------------------------","\n")
  
  # Prints out the values of the Sigma parameter implemented in the model
  
  cat("Satiation :","\n")
  cat("","\n")
  print(NDNDData$Sgma)
  cat("------------------------------","\n")                                      # Prints out the vector of Sigma parameter values 
  
  # We estimate and print out the food requirements with stable biomass for 2 years
  
  cat("Food requirements with stables biomass (no predation) : ","\n")
  cat("","\n")
  print(NDNDData$Mu/NDNDData$Gama)                                                # Prints the estimated values of food needed for each species at t+1 when biomass is constant
  cat("------------------------------","\n")
  cat("Food requirements in year 2 with stable biomass (no predation) : ","\n")
  cat("","\n")
  print((NDNDData$Mu/NDNDData$Gama)*NDNDData$Biomass)                             # Prints the estimated values of food needed for each species at t+2 when biomass is constant
  cat("------------------------------","\n")
  
  # We estimate and print out the food requirements when the biomass decrease at maximum rate
  
  cat("Food requirements with max decrease in biomass (no predation) : ","\n")
  cat("","\n")
  print((NDNDData$Mu/NDNDData$Gama)*(exp(-NDNDData$Rho)-exp(-NDNDData$Mu))/(1-exp(-NDNDData$Mu)))          # Prints the estimated values of food needed for each species when biomass is decreasing at its maximum rate
  cat("------------------------------","\n")
  
  
  # Is satiation values high enough for calculations?
  # In case it is too low for one or more species, then a warning message should be printed.
  # This test does not work for the phytoplankton which does not have a value for satiation since it does not predate.
  
  cat("","\n")
  satiation.setup<-NDNDData$Mu/NDNDData$Gama>NDNDData$Sigma                       # Logical test : is the satiation lower than the values of ingestion related parameters?
  F<-satiation.setup[which(satiation.setup==T)]                                   # Bolean output of the test
  if (length(F) > 0)  {                                                           # If one T --> satiation values too low
    cat("Warning : Satiation is set too low for the following species","\n")      # Prints out a warning message
    print(names(F))                                                               
  } else { cat("All satiation values high enough for calculation","\n")}          # Else it prints an information message
  cat("","\n")
  cat("------------------------------","\n")
  
  # Prints out the values of refuge biomass for each species
  
  cat("Refuge Biomass (Beta) : ","\n")
  cat("","\n")
  print(NDNDData$Bta)                                                            # Prints out the vector of refuge biomass for all species
  cat("------------------------------","\n")
  
  # Prints out the Import values
  
  cat("Biomass import in the system : ","\n")
  cat("","\n")
  print(head(NDNDData$Importall))                                                 # Prints out the Import values for each year of the simulation
  cat("","\n")
  cat("Mean Import : ","\n")
  print(apply(NDNDData$Importall,2,mean))                                         # Prints out the mean Import values for each species
  cat("","\n")
  cat("Variance of import : ","\n")
  print(apply(NDNDData$Importall,2,var))                                          # Prints out the variance of Import values for each species
  cat("------------------------------","\n")
  
  # Export values are presented the same way
  
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
  
  # Prints out the topology of the food web : the existing trophic links
  
  cat("","\n")
  cat("Who eats whom? : Trophic flows matrix (Preys in rows / Predators in columns)","\n")
  cat("","\n")
  print(NDNDData$Species)                                                    # Prints out the species names
  cat("","\n")
  print(NDNDData$PF)                                                              # Prints out the flows matrix
  
  flows_2<-NULL
  for (i in 1:NDNDData$ns){
    for (j in 1:NDNDData$ns) {
      flux<-paste(rownames(as.matrix(NDNDData$PF[i,])),"-->",colnames(as.matrix(NDNDData$PF[i,]))[j])
      flows_2<-cbind(flows_2,flux)
    }
  }
  flows<-flows_2[which(NDNDData$PFv==1)]
  
  cat("","\n")
  for (i in 1:length(flows)) {
    cat(flows[i],"\n")
  }
  
  sink()                                                                          # Closing the file containing the report
}