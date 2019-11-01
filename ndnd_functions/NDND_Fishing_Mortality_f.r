###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Fishing Mortality function
#### Version v1.1
#### 22.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

##### Fishing Mortality function
# We estimate the fishing mortality based on estimated values of fishing mortality rates and the Baranov equations

Phi.Mortality<-function(Out.file){
  
  # Extract the data needed to estimate the landings
  Biomass<-as.matrix(Out.file$Output$Biomass)                                # Biomass output
  sim_l<-as.numeric(length(table(Biomass[,ncol(Biomass)])))                     # Number of simulations
  Mu<-Out.file$Data$Mu                                                          # Metabolic losses
  Sp<-Out.file$Data$Species                                                     # List of species names
  Phi<-Out.file$Output$Phi                                                   # Fishing mortality rates output
  Tmax<-Out.file$Data$Tmax                                                      # Length of one simulation
  BI<-Out.file$Burn.In                                                          # Burn-in value 
  
  # Create subsets of the initial biomass outputs
  sub<-NULL                                                                     # Creates an empty object : sub
  
  # We don't take the last biomass values for each simulation ; We don't have fishing mortality rates for these biomass values'
  # sub is a subset of the initial biomass data that doesn't take into account the last year of biomass
  for (i in 1:sim_l) {
    a<-Biomass[as.matrix(Biomass[,ncol(Biomass)])==as.matrix(1:sim_l)[i],]
    b<-a[1:(nrow(a)-1),]
    sub<-rbind(sub,b)
  }
  
  # Mu is given as a vector
  # We repeat the vector so we get a matrix of the same dimension as sub
  Mu<-rep(Mu,times=nrow(sub))
  Mu<-matrix(Mu,nrow = nrow(sub),ncol = length(Sp),byrow = T)
  
  # Computing the fishing mortality based on the Baranov equations
  Mu_p_Phi<-abind(Mu,Phi[,1:length(Sp)],along=3)
  Mu_p_Phi<-apply(Mu_p_Phi,c(1,2),sum)
  
  c<-Phi[,1:length(Sp)]/Mu_p_Phi
  
  d<-(1-exp(-(Mu_p_Phi)))
  
  e<-c*d
  
  F.Mortality<-e*sub[,1:length(Sp)]
  
  # Adding a column containing the number fo simulation per statistical individual
  Simulation<-sub[,ncol(sub)]
  Years<-rep(1:((Tmax-1)-BI),times=as.numeric(length(table(Biomass[,ncol(Biomass)]))))
  
  F.Mortality<-cbind(F.Mortality,Simulation,Years)
  
  # Keep only the species for which there is fishing
  F.Mortality<-F.Mortality[,which(apply(F.Mortality, 2, sum)>0)]
  FM_no_ice<-cbind((F.Mortality[,1:2]*1400000),F.Mortality[,ncol(F.Mortality)-1],F.Mortality[,ncol(F.Mortality)])
  FM_ice<-cbind((F.Mortality[,1:2]*1600000),F.Mortality[,ncol(F.Mortality)-1],F.Mortality[,ncol(F.Mortality)])
  # colnames(FM_no_ice[,(ncol(FM_no_ice)-1):ncol(FM_no_ice)])<-c("Simulation","Years")
  
  # Data manipulation for plotting
  FM_gg<-gather(as.data.frame(FM_no_ice),key="Species",value="F.Mortality",1:(ncol(F.Mortality)-2))
  # FM_gg$V3<-as.factor(FM_gg$V3)
  
  # FM_gg$Simulation<-factor(FM_gg$Simulation,levels=1:sim_l)
  
  # Plotting fishing mortality 
  gg_FM<-ggplot(FM_gg,aes(x=as.factor(V4),y=F.Mortality,group=V3),label=FALSE)+
    geom_line()+
    facet_wrap(.~factor(Species,levels = Out.file$Data$Species[5:6]),ncol=1,scale="free")+
    scale_x_discrete(breaks = seq(min(FM_gg$V4),max(FM_gg$V4),50))+
    # scale_color_brewer(palette = "Paired")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "gray"))
  
  FM_gg_ice<-gather(as.data.frame(FM_ice),key="Species",value="F.Mortality",1:(ncol(F.Mortality)-2))
  
  gg_FM_ice<-ggplot(FM_gg_ice,aes(x=as.factor(V4),y=F.Mortality,group=V3))+
    geom_line(aes(linetype=as.factor(V3),color=Species),size=1)+
    facet_wrap(.~factor(Species,levels = Out.file$Data$Species[5:6]),ncol=1,scale="free")+
    scale_x_discrete(breaks = seq(min(FM_gg_ice$V4),max(FM_gg_ice$V4),50))+
    # scale_color_brewer(palette = "Paired")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "gray"))
  
  return(list(Sim.tag=Out.file$Sim.tag,
              Outputs=Out.file$Output,
              Fishing.Mortality=FM_no_ice,
              Fishing.Mortality.Ice=FM_ice,
              TimeSeries.FM=gg_FM))
}