###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Simulation initialization
#### Version v1.0
#### 02.12.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

# 0. Libraries ------------------------------------------------------------

# Needed libraries
###
require(inline)                             # Loading the inline package -- Read C++ functions called in R
require(Rcpp)                               # Loading the Rcpp package -- Integration of R and C++
require(RcppEigen)                          # Loading the RcppEigen package -- Linear algebra on C++ integrated in R
require(rstudioapi)                         # Loading the rstudioapi package -- Directory and file interactive choice
###

# 1. Locate Rtools --------------------------------------------------------

# The code needs to run parts in C++ to increase the speed of simulation
# Therefore, the some functions are coded in c++ and need to be compiled
# If you are working with Linux or Mac, the compiler is included, you don't need to include any other package
# If you are working on Windows, Rtools needs to be downloaded and installed on the computer
# https://cran.r-project.org/bin/windows/Rtools/
# Once installed you need to add the path to Rtools to the environment path with the function sys.setenv()

# Needs to be generalized
Sys.setenv(PATH = "C:/Rtools/bin;
           C:/Program Files/R-3.6.1/bin/x64;
           C:/Program Files(x86)/Common Files/Oracle/Java/javapath;
           C:/WINDOWS/system32;
           C:/WINDOWS; 
           C:/WINDOWS/System32/Wbem;
           C:/WINDOWS/System32/WindowsPowerShell/v1.0/;
           C:/Program Files/Intel/WiFi/bin/;
           C:/Program Files/Common Files/Intel/WirelessCommon/;
           C:/Program Files/Git/cmd;
           C:/Users/Administrator/AppData/Local/Microsoft/WindowsApps;
           C:/Users/Administrator/AppData/Local/Box/Box Edit/")

# 2. Clear graphic window, console and work environment -------------------

# Initialize the simulation : Clear your work environment
graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

# 3. Define the initial work directory ------------------------------------

# Choose your work directory
# In a project, the work directory is automatically set to the location of the project file
# Choose the directory where all NDND folders are located
setwd(selectDirectory(caption = "Select your work directory : "))                  # Chose your initial work directory
wd<-getwd()                                                                        # Save initial work directory

# 4. Load NDND functions --------------------------------------------------

# Load all functions
NDNDfunctions<-list.files(paste(wd,"/ndnd_functions",sep=""))                      # List the functions in the folder ndnd_functions

# Source the functions in the folder ndnd_functions
for (i in 1:length(NDNDfunctions)){
  function2source<-paste(wd,"/ndnd_functions","/",NDNDfunctions[i],sep="")
  source(function2source)
}

# 5. Load all directories -------------------------------------------------

# Create the table with all directories containing files and where outputs are saved
# It used later to navigate through the NDND directories
dir_tab<-directories(wd=wd)                    # Save all directories available from the initial directory in a vector
