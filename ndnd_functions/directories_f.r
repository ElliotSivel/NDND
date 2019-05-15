###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Directories function
#### Version v1.0
#### 07.05.19
#### Author : Elliot Sivel
###############################################################################################

##### Directories function
# Loads the directories available in the initial wd and gives a vector containing all directories and their names

directories<-function(wd){
  dir_tab=NULL
  dir_tab$code_dir=paste(wd,'/ndnd_code',sep="")                                      # Sets the directory to the folder where data files are located
  dir_tab$config_dir=paste(wd,'/ndnd_config',sep="")                                  # Sets the directory to the folder where data files are located
  dir_tab$configreport_dir=paste(wd,'/ndnd_configreport',sep="")                      # Sets the directory to the folder where data files are located
  dir_tab$data_dir=paste(wd,'/ndnd_data',sep="")                                      # Sets the directory to the folder where data files are located
  dir_tab$files_dir=paste(wd,'/ndnd_files',sep="")                                    # Sets the directory to the folder where data files are located
  dir_tab$functions_dir=paste(wd,'/ndnd_functions',sep="")                            # Sets the directory to the folder where data files are located
  dir_tab$outputs_dir=paste(wd,'/ndnd_outputs',sep="")                                # Sets the directory to the folder where data files are located
  return(dir_tab)
}
