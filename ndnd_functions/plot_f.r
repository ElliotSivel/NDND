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

library(ggplot2)

plot.NDND<-function(BiomassSeries,Tmax,Species,ns,Plotting){
  if (Plotting==1){
    plot_biom_temp<-as.numeric(as.vector(BiomassSeries))
    plot_species<-rep(Species,byrow=F,each=Tmax)
    plot_years<-rep(1:Tmax,times=ns)
    plot_biom<-as.data.frame(cbind(plot_years,plot_biom_temp,plot_species))
    plot<-ggplot(plot_biom,aes(x=as.numeric(plot_years),y=as.numeric(plot_biom_temp)))+geom_line()+facet_wrap(~plot_species,ncol=3)
  } else {print("Plotting parameter =0 ; no plotting configurated")}
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
