###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Shifts Number function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# This function computes the Number of shifts between system configurations

NDND.Nb.Shifts<-function(Output,Out.Cluster){
  
  Nb_shifts<-0                                                   # We set a counter, It is set to 0 at the beginning
  
  # We only take rows for wihich we have biomass values
  a<-drop_na(as.data.frame(Out.Cluster$Biomass.cluster))                                # We delete the rows where we have NAs
  
  for (i in 1:(nrow(a)-1)){
    if (a[i,9]==1 && a[i+1,9]==2){                                  # Test for shift from configuration 1 to configuration 2
      Nb_shifts<-Nb_shifts+1
    } else if (a[i,9]==2 && a[i+1,9]==1){                          # Test for shift from configuration 2 to configuration 1
      Nb_shifts<-Nb_shifts+1
    } else {Nb_shifts=Nb_shifts}                             # If no shift, the counter don't change
  }

  a<-array(Out.Cluster$Biomass.cluster[,9],dim=c(dim(Output$with.NaN.BI$Biomass)[1],dim(Output$with.NaN.BI$Biomass)[3]))
  b<-abind(Output$with.NaN.BI$Biomass,a,along=2)
  
  Nb_shifts_sim<-matrix(NA,nrow = 1, ncol = dim(Output$with.NaN.BI$Biomass)[3])
  
  for (j in 1:dim(Output$with.NaN.BI$Biomass)[3]){
    for (k in 1:(dim(Output$with.NaN.BI$Biomass)[1]-1)){
      if (b[k,9,j]==1 & b[k+1,9,j]==2){                                  # Test for shift from configuration 1 to configuration 2
        Nb_shifts_sim<-Nb_shifts_sim+1
      } else if (b[k,9,j]==2 && b[k+1,9,j]==1){                          # Test for shift from configuration 2 to configuration 1
        Nb_shifts_sim<-Nb_shifts_sim+1
      } else {Nb_shifts_sim=Nb_shifts_sim}                             # If no shift, the counter don't change
    }
  }
  
  return(list(Nb.Shifts=Nb_shifts,Cluster=b,Nb.Shifts.Sim=Nb_shifts_sim))
}