###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### readDATA function
#### Version v1.0
#### 27.03.19
#### Author : Elliot Sivel
###############################################################################################

##### readData function
# This function takes the NDNDConfig file to find the data files and load them into R
# After loading the data, it assemble them in a NDNDData object

readDATA<-function(config_file,files_dir){
  file.names<-scan(config_file,what=character())                     # Reading the elements given by the NDNDConfig file
                                                                    # file.names is object containing the elements of the NDNDConfig file as character             
  current_wd=getwd()
  setwd(files_dir)
  # Test the presence of the data files needed in the work directory
  
  exist.file<-vector(length=length(file.names)-2)                   # Create an empty vector of length of number of files
                                                                    # "-2" because the two last elements of the NDNDConfig file are the length of simulation and the plotting parameter
                                                                    # We don't need them now
  for (i in (1:(length(file.names)-2))) {                           # Loop on file.names to test if each file is existing in the work directory                
    exist.file[i]<-file.exists(file.names[i])                       # Output a vector of bolean values (FALSE : the file is missing ; TRUE : the file is there)
  }
  
  # Load the files
  # We want the files to be loaded if all files are available
  # We want the function to return an error message if at least one of the files is missing
  
  if (sum(exist.file)==length(file.names)-2){                                                               # Condition for the model to run : the sum of exist.file = length of file.names
                                                                                                            # If the condition is verified :
    Species<-scan(file=file.names[1],what=character(),quote="")                                             # Read species.txt
    PF<-read.table(file.names[2], header=F, sep="\t",stringsAsFactors = F)                                  # Read fluxes.txt
    Coefficients<-as.matrix(read.table(file.names[3], header=F, sep="\t",stringsAsFactors = F))             # Read coefs.txt
    Import<-scan(file.names[4], what=numeric())                                                             # Read import.txt
    Export<-scan(file.names[5], what=numeric())                                                             # Read export.txt
  }  else {stop("Missing file")}                                                                            # If the condition is not verified : Print error message       
  
  # Files are now loaded
  # Defining parameters and vectors from NDNDConfig and the 5 data files
  
  Tmax<-as.numeric(file.names[6])             # Defining the length of the simulation : Tmax
  ns<-as.numeric(length(Species))             # Defining the number of species
  nn<-ns*ns                                   # Defining the number of possible flows in the food web
  
  # Formating Data : give names to the data frames 
  
  colnames(Coefficients)<-Species                                                        # Set the species names as column names in the coefficient dataframe
  rownames(Coefficients)<-c("Biomass","Gama","Kapa","Mu","Rho","Sigma","Bta")            # Set the names of the parameters as row names in the coefficient dataframe
  colnames(PF)<-Species                                                                  # Set the species names as column names in the flows matrix
  rownames(PF)<-Species                                                                  # Set the species names as row names in the flows matrix

  # Spliting the Coefficient matrix into single parameters values
  
  Biomass<-Coefficients[1,]               # Row 1 is the initial biomass for each species, it is given in tons per km2
  Gama<-Coefficients[2,]                  # Row 2 is Gama, the potential assimiliation efficiency for each species
  Kapa<-Coefficients[3,]                  # Row 3 is Kapa, the digestability of each species
  Mu<-Coefficients[4,]                    # Row 4 is Mu, the metabolic losses for each species
  Rho<-Coefficients[5,]                   # Row 5 is Rho, the parameter accounting for inertia
  Sgma<-Coefficients[6,]                  # Row 6 is Sigma, the parameter accounting for satiation. It has been renamed Sgma because Sigma is already a function in R
  Bta<-Coefficients[7,]                   # Row 7 is Beta, the value for refuge biomass
  
  # Import and Export values can be in vectors or matrices
  # If in a vector : replication of the vector for the length of the simulation (Tmax)
  # If in a matrix with nrow > Tmax : Importall = Import[1:Tmax,]
  # If in a matrix with nrow < Tmax : Return an error message
  # If in a matrix with nrow = Tmax : Importall = Import
  # Identical for Export
 
  if (nrow(as.matrix(t(Import)))==1){
    Importall<-matrix(rep(as.matrix(t(Import)),each=Tmax),nrow=Tmax)                       
  } else if (nrow(as.matrix(t(Import)))==Tmax){
    Importall<-Import                                                                       
  } else if (nrow(as.matrix(t(Import))) > Tmax){                      
    Importall<-Import[1:Tmax,]                                                              
    stop("1 < number of rows < Tmax")                                                       
  } else {stop("Import data missing")}                                               
  
  if (nrow(as.matrix(t(Export)))==1){
    Exportall<-matrix(rep(as.matrix(t(Export)),each=Tmax),nrow=Tmax)                        
  } else if (nrow(as.matrix(t(Export)))==1){
    Exportall<-Export                                                                       
  } else if (nrow(as.matrix(t(Export))) > Tmax){
    Exportall<-Export[1:Tmax,]                                                             
  } else if (nrow(as.matrix(t(Export))) < Tmax && nrow(as.matrix(t(Export))) != 1){
    stop("1 < number of rows < Tmax")                                                       
  } else {stop("Export data missing")}                                                      
  
  # Formating Importall and Exportall : setting colnames and rownames
  
  colnames(Importall)<-Species               # Set the species names as column names in Importall matrix
  rownames(Importall)<-1:Tmax                # Set the years as row names in the Importall matrix
  colnames(Exportall)<-Species               # Set the species names as column names in Exportall matrix
  rownames(Exportall)<-1:Tmax                # Set the years as row names in the Exportall matrix
  
  # Rest of the model uses matrix calculation
  # Easier to work with vectors than with matrices
  # Vectorization of the possible flows matrix
  
  PFv=as.vector(t(PF))                  # Transform the possible flows matrix in a single vector of length nn
  
  # Defining if the time series are going to be plotted
  
  Plotting<- file.names[7]                   # Binary parameter (0 = no plot ; 1 = plot)               
  
  # Merging all created elements na list
  
  NDNDData<-list(Tmax,Species,ns,nn,Biomass,Gama,Kapa,Mu,Rho,Sgma,Bta,Importall,Import,Exportall,Export,PF,PFv,Plotting)
  
  # Formating the list NDNDData : Adding the names of the elements in the list
  
  names(NDNDData)<-c("Tmax","Species","ns","nn","Biomass","Gama","Kapa","Mu","Rho","Sigma","Beta","Importall","Import","Exportall","Export","PF","PFv","Plotting")
  
  # Give the elements the function should return
  setwd(current_wd)
  return(NDNDData)
}
