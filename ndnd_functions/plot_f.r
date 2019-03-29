BiomassSeries_v<-as.vector(BiomassSeries)
Species_plot<-rep(NDNDData$Species,byrow=F,each=NDNDData$Tmax)
Tmax_v<-as.numeric(rep(1:NDNDData$Tmax,by=NDNDData$ns))
Biomass_plot<-as.data.frame(cbind(Tmax_v,BiomassSeries_v,Species_plot))



library(ggplot2)
library(gridExtra)
ggplot(Biomass_plot,aes(x=Tmax_v,y=BiomassSeries_v))+geom_line()+facet_grid(BiomassSeries_v~Species_plot)

grid.arrange(ggplot(as.data.frame(BiomassSeries[,1]))+aes(1:NDNDData$Tmax,BiomassSeries[,1])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[1]),
             ggplot(as.data.frame(BiomassSeries[,2]))+aes(1:NDNDData$Tmax,BiomassSeries[,2])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[2]),
             ggplot(as.data.frame(BiomassSeries[,3]))+aes(1:NDNDData$Tmax,BiomassSeries[,3])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[3]),
             ggplot(as.data.frame(BiomassSeries[,4]))+aes(1:NDNDData$Tmax,BiomassSeries[,4])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[4]),
             ggplot(as.data.frame(BiomassSeries[,5]))+aes(1:NDNDData$Tmax,BiomassSeries[,5])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[5]),
             ggplot(as.data.frame(BiomassSeries[,6]))+aes(1:NDNDData$Tmax,BiomassSeries[,6])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[6]),
             ggplot(as.data.frame(BiomassSeries[,7]))+aes(1:NDNDData$Tmax,BiomassSeries[,7])+geom_line()+xlab("Years")+ylab(expression(Biomass (tons/ ~ km^2 )))+ggtitle(NDNDData$Speciesnames[7]),
             ggplot(as.data.frame(BiomassSeries[,8]))+aes(1:NDNDData$Tmax,BiomassSeries[,8])+geom_line()+xlab("Years")+ylab("")+ggtitle(NDNDData$Speciesnames[8]),
             ncol=3,nrow=3)