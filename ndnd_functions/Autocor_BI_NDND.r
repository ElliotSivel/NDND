###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Autocorrelation function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

##### Autocorrelation exploration function
# We estimate the autocorrelation for the time series to define the right burn-in

NDND.auto.cor<-function(Out.file){
  
  # Extract information we need from the file Output
  Biomass<-Out.file$Output_NaN_0$Biomass
  Tmax<-Out.file$Data$Tmax
  Sp<-Out.file$Data$Species
  sim_l<-as.numeric(length(table(Biomass[,ncol(Biomass)])))
  
  # Create an three dimensional array that should contain the autocorrelation values for each species and each simulation 
  AC<-NULL
  for (i in 1:(ncol(Biomass)-1)){
    correl.AC<-by(data=Biomass[,i],
                  INDICES = Biomass[,ncol(Biomass)],
                  FUN = acf,
                  lag.max=nrow(Biomass),
                  type = "correlation",
                  plot=F)
    ac<-NULL
    for (j in 1:sim_l){
      ac<-rbind(ac,as.matrix(correl.AC[[j]]$acf))
    }
    AC<-cbind(AC,ac)
  }
  
  AC<-cbind(AC,
            as.matrix(rep(1:sim_l,byrow=F,each=Tmax)),
            as.matrix(rep(1:Tmax,times=sim_l)))
  colnames(AC)<-c(Sp,"Simulation","Years")
  
  # Array manipulation for plotting correlogram
  AC_gg<-gather(as.data.frame(AC),key = "Species",value = "Correlation",1:8)
  
  # Estimation of significance level for the dataset
  sig_AC<-qnorm((1+0.975)/2)/sqrt(Tmax)
  
  # plotting correlograms by species
  AC_plot<-ggplot(AC_gg,aes(x=Years,y=Correlation,group=Species))+
    geom_line(aes(color=Species),size=0.1)+
    geom_line(y=0)+
    geom_line(y=sig_AC,col="blue",linetype=2)+
    geom_line(y=-sig_AC,col="blue",linetype=2)+
    facet_wrap(.~factor(Species,levels = Sp),ncol=3)+
    scale_x_discrete(breaks = seq(0,Tmax,Tmax/10))+
    scale_color_brewer(palette = "Paired")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "gray"))+
    ggtitle("Temporal autocorrelation correlograms")
  
  # Defining for which lag the autocorrelation gets non-significant
  # Creating an array in which we will fill the first lag value for which the autocorrelation is not significant
  SIG_val<-array(NA,
                 dim=c(sim_l,
                       length(Sp)),
                 dimnames=list(1:sim_l,
                               Sp))
  
  # We fill the SIG_val matrix
  for (k in 1:sim_l) {
    sub<-AC_gg[as.matrix(AC_gg$Simulation)==as.matrix(1:sim_l)[k],]
    for (l in 1:length(Sp)) {
      sub2<-sub[sub$Species==Sp[l],]
      SIG_val[k,l]<-min(sub2$Years[sub2$Correlation<sig_AC])
    }
  }
  
  # We define the burn-in as the highest lag value for which the autocorrelation is not significant
  Burn_In <- max(SIG_val)
  
  # Create the list of outputs
  NDND_auto_cor<-list(Sim.tag=Output$Sim.tag,                                     # Simulation tag for identifying the simulation
                      Autocorrelation = AC,                                       # The autocorrelation values array
                      Autocorrelation_significance_level = sig_AC,                # Value of significance of autocorrelation
                                                                                  # It's a single value because the significance levels depends only on the number of statistical samples
                      Burn_In = Burn_In,                                          # Burn-In value 
                      Correlogram = AC_plot)                                      # ggplots correlograms
  
  return(NDND_auto_cor)
}