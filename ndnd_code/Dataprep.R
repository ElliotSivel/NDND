# Temporary script to prepare NDND data


# 1. Time stamp --------------------------------------------------------------

Data.tag=Sys.time()

# 2. working directories -----------------------------------------------------

directories=NULL
wd<-selectDirectory(caption = "Select your work directory : ")
directories$code_dir=paste(wd,'/ndnd_code',sep="")                                             # Sets the directory to the folder where data files are located
directories$config_dir=paste(wd,'/ndnd_config',sep="")                                             # Sets the directory to the folder where data files are located
directories$configreport_dir=paste(wd,'/ndnd_configreport',sep="")                                             # Sets the directory to the folder where data files are located
directories$data_dir=paste(wd,'/ndnd_data',sep="")                                             # Sets the directory to the folder where data files are located
directories$files_dir=paste(wd,'/ndnd_files',sep="")                                             # Sets the directory to the folder where data files are located
directories$functions_dir=paste(wd,'/ndnd_functions',sep="")                                             # Sets the directory to the folder where data files are located
directories$outputs_dir=paste(wd,'/ndnd_outputs',sep="")                                             # Sets the directory to the folder where data files are located

# 3. source functions -----------------------------------------------------

NDNDfunctions<-list.files(directories$functions_dir)
for (i in 1:length(NDNDfunctions)){
  function2source<-paste(directories$functions_dir,"/",NDNDfunctions[i],sep="")
  source(function2source)
}

# 4. load data ---------------------------------------------------------------

# There are 6 .txt files needed to run the model
# Data files (5) : species.txt, fluxes.txt, coefs.txt, import.txt, export.txt
# Configuration file : NDNDConfig -- Names of the data files, length of simulation and plotting parameter

# Load the Configuration file
# Choice of the file is set as interactive

setwd(directories$config_dir)                                                            # Sets the directory to the folder where configuration files are located
config_file<-selectFile(caption = "Select your configuration file : ")             # Opens a window to choose the configuration file your want to implement.

# Loading the data with readDATA

NDNDData<-readDATA(directories = directories,config_file = config_file,Data.tag = Data.tag)                       # Applies the readDATA function

# 5. save data ------------------------------------------------------------

save(NDNDData,file=paste(directories$data_dir,"/NDNDData_",
                         paste(format(Data.tag,"%Y_%m_%d_%H_%M_%S"),sep = ""),
                         ".RData",sep=""))             # Save NDNDData with time stamp
save(NDNDData,file="NDNDData.RData")                   # Save NDNDData without time stamp

# 6. Report Model Configuration --------------------------------------------

# Create a report of the input data and parametrization

NDNDConfigreport(directory = directories$configreport_dir, NDNDData = NDNDData)                         # Applies the NDNDConfigreport function
# Creates the configuration report

