###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### packages.version.NDND function
#### Version v1.0
#### 08.04.19
#### Author : Elliot Sivel
###############################################################################################

#### Package.version.NDND function
# Loading the version of all packages that have been loaded
# Function needs the session information report

packages.version.NDND<-function(Session_info=session_info){
  added_packages<-Session_info$otherPkgs                                 # Isolates the non-basic packages loaded in the R session
  packages_version<-matrix(0,ncol=2,nrow=length(added_packages))         # creates a matrix of nrow = the number of packages loaded and 2 columns
  
  # Fills the packages_version matrix with the names and the corresponding packages version number
  
  for (p in 1:length(added_packages)) {
    packages_version[p,1]<-added_packages[[p]][[1]]
    packages_version[p,2]<-added_packages[[p]][["Version"]]
  }
  
  return(packages_version)
}