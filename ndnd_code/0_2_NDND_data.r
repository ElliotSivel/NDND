###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Model Data
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

# 0. Libaries -------------------------------------------------------------

# Needed libraries
###
require(yesno)
###

# 1. Initialization -------------------------------------------------------

# Source the initialization file
source(file = "C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_code/0_1_NDND_initialization.R") # Runs the initialization file

# 2. Load or create NDNDData file -----------------------------------------

# To run the simulation, you need a NDNDData file computed from the Dataprep.r file
# If you change parameter valeurs, you need to rerun Dataprep.r with the new parameter values
test_ndnd_data<-yesno2("Do you have a NDNDDate file containing the setup for the simulation?", yes = "y", no="n")

# Load the NDNDData or run the Dataprep.r code
if(test_ndnd_data==F){                                                        # Do you have a NDNDData file
  source('./ndnd_code/0_3_NDND_dataprep.r')                                   # If not, run the Dataprep.r file
} else {load(file=paste(dir_tab$data_dir,"/NDNDData.RData",sep=""))}          # If yes, just load the existing NDNDData file
