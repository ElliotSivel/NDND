###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Plot function
#### Version v1.0
#### 05.04.19
#### Author : Elliot Sivel
###############################################################################################

#### Plot.NDND function
# Function for plotting the time series of biomass
# For plotting, it uses the ggplot2 package

plot.NDND<-function(NDNDSimulation){

  df_b<-gather(NDNDSimulation$Output$BiomassSeries,key = "Species" , value = "Biomass")
  gg2_b<-ggplot(df_b,aes(x=rep(1:NDNDSimulation$Data$Tmax,times=NDNDSimulation$Data$ns),y=Biomass))+
    geom_line()+
    facet_wrap(~factor(Species,levels = NDNDSimulation$Data$Species),ncol=3,scale="free")+
    xlab("Years")
  
  df_flow<-gather(NDNDSimulation$Output$FlowSeries,key = "Flow" , value = "Biomass")
  gg2_flow<-ggplot(df_flow,aes(x=rep(1:(NDNDSimulation$Data$Tmax-1),times=length(NDNDSimulation$Data$flows)),y=Biomass))+
    geom_line()+
    facet_wrap(~factor(Flow,levels = NDNDSimulation$Data$flows),ncol=5,scale="free")+
    xlab("Years")
  
  df_fu<-gather(NDNDSimulation$Output$FuSeries,key = "Species" , value = "Fu")
  gg2_fu<-ggplot(df_fu,aes(x=rep(1:(NDNDSimulation$Data$Tmax-1),times=NDNDSimulation$Data$ns),y=Fu))+
    geom_line()+
    facet_wrap(~factor(Species,levels = NDNDSimulation$Data$Species),ncol=3)+
    xlab("Years")

 gg<-list(gg2_b,gg2_flow,gg2_fu)
 return(gg)
}







# library(gridExtra)
# grid.arrange(ggplot(as.data.frame(BiomassSeries[,1]))+aes(1:NDNDData$Tmax,BiomassSeries[,1])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[1]),
#              ggplot(as.data.frame(BiomassSeries[,2]))+aes(1:NDNDData$Tmax,BiomassSeries[,2])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[2]),
#              ggplot(as.data.frame(BiomassSeries[,3]))+aes(1:NDNDData$Tmax,BiomassSeries[,3])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[3]),
#              ggplot(as.data.frame(BiomassSeries[,4]))+aes(1:NDNDData$Tmax,BiomassSeries[,4])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[4]),
#              ggplot(as.data.frame(BiomassSeries[,5]))+aes(1:NDNDData$Tmax,BiomassSeries[,5])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[5]),
#              ggplot(as.data.frame(BiomassSeries[,6]))+aes(1:NDNDData$Tmax,BiomassSeries[,6])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[6]),
#              ggplot(as.data.frame(BiomassSeries[,7]))+aes(1:NDNDData$Tmax,BiomassSeries[,7])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[7]),
#              ggplot(as.data.frame(BiomassSeries[,8]))+aes(1:NDNDData$Tmax,BiomassSeries[,8])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[8]),
#              ncol=3,nrow=3)
