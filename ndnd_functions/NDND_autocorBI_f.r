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

# 0. Needed libraries -----------------------------------------------------
###
require("reshape2")
require("ggplot2")
require("dplyr")
require("tidyr")
require("ggpubr")

###

NDND.auto.cor<-function(Output){

# 1. Extraction of data ---------------------------------------------------

  # Extract information we need from the file Output
  Biomass<-log10(Output$Biomass[apply(Output$Biomass,1,sum)-Output$Biomass[,ncol(Output$Biomass)]!=0,])                                         # Extract the Biomass time series
  Tmax<-Output$Data$Tmax                                          # Extract the length of the time series
  Sp<-Output$Data$Species                                         # Extract the vector of Species names
  # Biomass<-cbind(apply(Biomass[,1:length(Sp)], c(1,2), log10),
  #                "Simulation"= Biomass[,length(Sp)+1])
  sim_l<-as.numeric(length(table(Biomass[,ncol(Biomass)])))         # Extract the number of simulations

# 2. Compute the autocorrelation of the time series -----------------------

  # Create an three dimensional array that should contain the autocorrelation values for each species and each simulation 
  AC<-NULL                                                    # Create an empty array for saving autocorrelation values
  for (i in 1:(ncol(Biomass)-1)){                             # Loop for autocorrelation
    correl.AC<-by(data=Biomass[,i],                           # Performing autocrrelation computation per species   
                  INDICES = Biomass[,ncol(Biomass)],
                  FUN = acf,
                  lag.max=(Tmax/2)-1,                         # The lag can only be on half of the simulation
                  type = "correlation",
                  plot=F)
    ac<-NULL
    for (j in 1:sim_l){
      ac<-rbind(ac,as.matrix(correl.AC[[j]]$acf))
    }
    AC<-cbind(AC,ac)
  }
  
  AC<-as.data.frame(cbind(AC,
            as.matrix(rep(1:sim_l,byrow=F,each=(Tmax/2))),
            as.matrix(rep(1:(Tmax/2),times=sim_l))))
  colnames(AC)<-c(Sp,"Simulation","Years")
  

# 3. Plotting correlograms ------------------------------------------------

  # Array manipulation for plotting correlogram
  # AC_gg<-gather(AC,key = "Species",value = "Correlation",1:8)
  AC_gg<-pivot_longer(AC[,-1],cols=-c(Simulation, Years),names_to = "Species",values_to = "Correlation")
  
  # Estimation of significance level for the dataset
  sig_AC_up<-NULL                            
  for (j in 1:(Tmax/2)){
    a<-qnorm((1+0.975)/2)/sqrt((Tmax/2)-j)
    sig_AC_up<-rbind(sig_AC_up,a)
  }

  sig_AC_up<-rep(sig_AC_up,times=length(Sp)-1)
  c<-rep(Sp[-1],byrow=F,each=(Tmax/2))
  d<-rep(1:(Tmax/2),times=length(Sp)-1)
  sig_AC_up<-cbind("Significance"=sig_AC_up,"Species"=c, "Years"=d)
  sig_AC_up<-as.data.frame(sig_AC_up)
  
  AC_gg_bis<-merge(AC_gg,sig_AC_up,by=c("Species","Years"),all.x=T)
  AC_gg_bis[is.infinite(AC_gg_bis$Significance)]<-NA

# 4. Identifying the significant values of correlation  -------------------

  # Defining for which lag the autocorrelation gets non-significant
  # Creating an array in which we will fill the first lag value for which the autocorrelation is not significant
  SIG_val<-array(NA,
                 dim=c(sim_l,
                       length(Sp)-1),
                 dimnames=list(1:sim_l,
                               Sp[-1]))
  
  # We fill the SIG_val matrix
  for (k in 1:sim_l) {
    sub<-AC_gg_bis[as.matrix(AC_gg_bis$Simulation)==as.matrix(1:sim_l)[k],]
    for (l in 1:(length(Sp)-1)) {
      sub2<-as.data.frame(sub[sub$Species==Sp[l+1],])
      sub3<-sub2[as.numeric(as.character(sub2$Correlation))<as.numeric(as.character(sub2$Significance)),]
      SIG_val[k,l]<-min(sub3$Years)
    }
  }

# 5. Estimate the burn-in -------------------------------------------------

  # We define the burn-in as the highest lag value for which the autocorrelation is not significant
  Burn_In <- round(max(apply(SIG_val,2,quantile,0.5,na.rm=T)),digits=0)
  
  # plotting correlograms by species
  listAC<-list()
  for (i in 1:(length(Sp)-1)){
    subggAC<-AC_gg_bis[AC_gg_bis$Species %in% Sp[i+1],]
    AC_plot<-ggplot(data=AC_gg_bis,aes(x=Years,group=Simulation))+
      geom_line(aes(y=Correlation),color="turquoise",size=1)+
      geom_hline(aes(yintercept=0))+
      geom_vline(aes(xintercept=Burn_In),linetype="dashed",color="black")+
      geom_line(aes(y=as.numeric(as.character(Significance))),col="blue",linetype=2,size=0.5)+
      geom_line(aes(y=-as.numeric(as.character(Significance))),col="blue",linetype=2,size=0.5)+
      ylim(-1.1,1.1)+
      theme_bw()+
      theme(legend.position = "none")
      labs(y="Correlation")
    listAC[[i]]<-AC_plot
    }
  
  ACPLOTGG<-ggarrange(listAC[[1]],listAC[[2]],listAC[[3]],listAC[[4]],listAC[[5]],listAC[[6]],listAC[[7]],labels = c("A","B","C","D","E","F","G"),nrow=4,ncol=2)
  
  # Create the list of outputs
  NDND_auto_cor<-list(Sim.tag=Output$Sim.tag,                                     # Simulation tag for identifying the simulation
                      Autocorrelation = AC,                                       # The autocorrelation values array
                      Autocorrelation_significance_level = sig_AC_up,             # Value of significance of autocorrelation
                      Burn_mat = SIG_val,                                         # It's a single value because the significance levels depends only on the number of statistical samples
                      Burn_In = Burn_In,                                          # Burn-In value 
                      Correlogram = ACPLOTGG)                                      # ggplots correlograms
  
  return(NDND_auto_cor)
}