###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Data preparation
#### Version v1.1.5
#### 15.05.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

# 1. Time stamp --------------------------------------------------------------

Data.tag=Sys.time()               # Define a tag / ID for the data -- Date and time of creation

# 2. load data ---------------------------------------------------------------

# There are 7 .txt files needed to run the model
# Data files () : species.txt, fluxes.txt, coefs.txt, coefs_MP.txt, import.txt, export.txt
# Configuration file : NDNDConfig -- Names of the data files, length of simulation and plotting parameter

# Load the Configuration file
# Choice of the file is set as interactive

setwd(dir_tab$config_dir)                                                     # Sets the directory to the folder where configuration files are located
config_file<-selectFile(caption = "Select your configuration file : ")        # Opens a window to choose the configuration file your want to implement.

# Loading the data with readDATA

NDNDData<-readDATA(directories = dir_tab,config_file = config_file,Data.tag = Data.tag)           # Applies the readDATA function

# 3. save data ------------------------------------------------------------

# Two NDNDData files are saved, one with stamp for archive and one without for the remaining simulation (Simplier to load)

save(NDNDData,file=paste(dir_tab$data_dir,"/NDNDData_",
                         paste(format(Data.tag,"%Y_%m_%d_%H_%M_%S"),sep = ""),
                         ".RData",sep=""))                                             # Save NDNDData with time stamp
save(NDNDData,file=paste(dir_tab$data_dir,"/NDNDData.RData",sep=""))                   # Save NDNDData without time stamp

# 4. Report Model Configuration --------------------------------------------

# Create a report of the input data and parametrization

NDNDConfigreport(directory = dir_tab$configreport_dir, NDNDData = NDNDData)            # Applies the NDNDConfigreport function
                                                                                       # Creates the configuration report

