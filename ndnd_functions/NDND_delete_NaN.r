###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Shift NaN to 0 function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# This function replaces NaN with 0 in the output data

NDND.NaN.to.zero<-function(simulation_output){
  simulation_output$Biomass[is.na(simulation_output$Biomass)]<-0
  simulation_output$Flows[is.na(simulation_output$Flows)]<-0
  simulation_output$Phi[is.na(simulation_output$Phi)]<-0
  return(simulation_output)
}