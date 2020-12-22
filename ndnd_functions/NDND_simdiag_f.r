###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Simulation diagnose function
#### Version v1.1
#### 22.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

##### Simulation diagnose function

require(purrr)

# We estimate the fishing mortality based on estimated values of fishing mortality rates and the Baranov equations

Sim.Diagnose<-function(Out.file){
  
# 1. Loading the data and parameters --------------------------------------

  # Extract the data needed
  Biomass<-as.matrix(Out.file$Biomass)                                   # Biomass output
  Flows<-Out.file$Flows                                                  # Flows output
  sim_l<-as.numeric(length(table(Biomass[,ncol(Biomass)])))              # Number of simulations
  Simulation<-Biomass[,ncol(Biomass)]                                    # Vector of simulation ID
  Tmax<-Out.file$Data$Tmax                                               # Length of one simulation
  Years<-rep(1:(Tmax-1),times=as.numeric(length(table(Biomass[,ncol(Biomass)])))) # Creating the vector of years
  Mu<-Out.file$Data$Mu                                                   # Metabolic losses
  Sp<-Out.file$Data$Species                                              # List of species names
  Kapa<-Out.file$Data$Kapa                                               # Prey digestability
  Gama<-Out.file$Data$Gama                                               # Predator assimilation capacity
  Imp<-Out.file$Data$Import                                              # Import data
  Export<-as.matrix(Out.file$Fish)                                           # Export data
  
  # Create subsets of the initial biomass outputs
  s_b<-Biomass %>% split(as.factor(Biomass[,ncol(Biomass)]))
  sub_bio<-map(s_b,function(x){
    sub<-matrix(x,ncol = ncol(Biomass))
    sub<-sub[1:(Tmax-1),]
  })
  sub_bio<-do.call(rbind.data.frame,sub_bio)
  sub_bio<-cbind(sub_bio,"Years"=Years)
  colnames(sub_bio)<-c(Sp,"Simulation","Years")
  
  # 2. Computing losses due to fishing --------------------------------------
  
  # Mu correspond to one value for each species
  # We replicate these values in order to have a matrix of the same size as the sub matrix
  Mu_m<-Mu %>% rep(times=nrow(sub_bio)) %>% matrix(nrow = nrow(sub_bio),ncol = length(Sp),byrow = T)
  
  # Computing metabolic losses
  E<-cbind((1-exp(-(Mu_m)))*sub_bio[,1:8],sub_bio[,9:10])
  
  # 4. Format parameter matrices  ------------------------------------------
  
  Kapa_m<-Kapa %>% rep(times=nrow(sub_bio)) %>% matrix(nrow = nrow(sub_bio),ncol = length(Sp),byrow = T)
  colnames(Kapa_m)<-Sp
  
  Gama_m<-Gama %>% rep(times=nrow(sub_bio)) %>% matrix(nrow = nrow(sub_bio),ncol = length(Sp),byrow = T)
  colnames(Gama_m)<-Sp
  
  Imp_m<-Imp %>% rep(times=nrow(sub_bio)) %>% matrix(nrow = nrow(sub_bio),ncol = length(Sp),byrow = T)
  
  s_e<-Export %>% split(as.factor(Export[,ncol(Export)]))
  sub_exp<-map(s_e,function(x){
    sub<-matrix(x,ncol = ncol(Export))
    sub<-sub[1:(Tmax-1),]
  })
  sub_exp<-do.call(rbind.data.frame,sub_exp)
  sub_exp<-cbind(sub_exp,"Years"=Years)
  colnames(sub_exp)<-c(Sp,"Simulation","Years")
  
  inv_k<-1-Kapa_m
  inv_g<-1-Gama_m
  
  # 5. Computing income of Biomass per Species per year ---------------------
  
  gains<-NULL
  for (s in 1:length(Sp)) {
    sub_s<-as.matrix(Flows[,grep(pattern=paste("* --> ",Sp[s],sep=""),x=colnames(Flows))])
    sum_sp<-as.numeric(apply(sub_s, 1, sum))
    gains<-cbind(gains,sum_sp)
  }
  
  gains<-gains+Imp_m
  colnames(gains)<-Sp
  
  # 6. Computing transfert of Biomass from prey to predator per Spec --------
  
  transfert<-NULL
  for (s in 1:length(Sp)) {
    sub_t<-as.matrix(Flows[,grep(pattern=paste(Sp[s]," --> *",sep=""),x=colnames(Flows))])
    sum_t<-as.numeric(apply(sub_t, 1, sum))
    transfert<-cbind(transfert,sum_t)
  }
  
  transfert<-transfert+sub_exp[,-c(9,10)]
  colnames(transfert)<-Sp
  
  # 7. Computing losses due to digestability or assimilation of prey --------
  
  # Compute losses due to non-assimilation of the prey
  Fk<-NULL
  for (s in 1:length(Sp)) {
    sub_l<-Flows[,grep(pattern=paste("*--> ",Sp[s],sep=""),x=colnames(Flows))]
    na_co<-colnames(Flows)[grep(pattern=paste("*--> ",Sp[s],sep=""),x=colnames(Flows))]
    na_co<-sub("-->.*","",na_co)
    kapa_sub<-Kapa_m[,colnames(Kapa_m) %in% trimws(c(na_co))]
    sub_l_k<-sub_l*kapa_sub
    sum_sp<-as.numeric(apply(as.matrix(sub_l_k), 1, sum))
    Fk<-cbind(Fk,sum_sp)
  }
  
  losses_assim<-Fk*inv_g
  
  # Compute losses due to non digestion of the prey
  Fg<-NULL
  for (s in 1:length(Sp)) {
    sub_l<-Flows[,grep(pattern=paste("*--> ",Sp[s],sep=""),x=colnames(Flows))]
    na_co<-sub("-->.*","",colnames(sub_l))
    kapa_sub<-Kapa_m[,colnames(Kapa_m) %in% trimws(c(na_co))]
    sub_l_k<-sub_l*(1-kapa_sub)
    sum_sp<-as.numeric(apply(sub_l_k, 1, sum))
    Fg<-cbind(Fg,sum_sp)
  }
  
  losses_dig<-Fg
  
# 8a. Computing Variation of Biomass based on the master equation ----------
  
  C<-(1-exp(-Mu_m))/Mu_m
  delta_B<-((gains-transfert-losses_assim-losses_dig)*C)-E[,-c(9,10)]

# 8b. Computing variation of Biomass using the simulation outputs ----------

  dif_B<-NULL
  for (j in 1:sim_l){
    sub_j<-Biomass[Biomass[,ncol(Biomass)]==as.matrix(1:sim_l)[j],]
    for (k in (1:nrow(sub_j))-1){
      dlta<-sub_j[k+1,]-sub_j[k,]
      dif_B<-rbind(dif_B,dlta)
    }
  }

# 9. Difference Simulation vs. estimated ----------------------------------

  DIF<-delta_B-dif_B

# 11. Return output of function -------------------------------------------
  
  return(list(Sim.tag=Out.file$Sim.tag,
              bio=sub_bio,
              fish=sub_exp,
              Other_losses=E,
              Gains=gains,
              Transfers=transfert,
              Digestion_losses=losses_dig,
              Assimilation_losses=losses_assim,
              ME_delta=delta_B,
              Sim_delta=dif_B,
              dif=DIF))
}