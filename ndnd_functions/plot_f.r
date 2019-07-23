###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Plot function
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel
###############################################################################################

#### Plot.NDND function
# Function for plotting the time series of biomass
# For plotting, it uses the ggplot2 package

plot.NDND<-function(Out.file){
  
  Biomass<-as.data.frame(Out.file$Output.BI$Biomass)
  Flows<-as.data.frame(Out.file$Output.BI$Flows)
  Phi<-as.data.frame(Out.file$Output.BI$Phi)
  Sp<-Out.file$Data$Species
  Tmax<-Out.file$Data$Tmax
  BI<-Out.file$Burn.In
  
  Years<-rep(1:(Tmax-BI),times=as.numeric(length(table(Biomass[,ncol(Biomass)]))))
  Biomass<-cbind(Years,Biomass)
  bio_gg<-gather(Biomass,key="Species",value="Biomass",2:(ncol(Biomass)-1))
 
  gg_bio<-ggplot(bio_gg,aes(x=as.factor(Years),y=Biomass,group=Simulation))+
    geom_line(aes(color=Species))+
    facet_wrap(.~factor(Species,levels = Out.file$Data$Species),ncol=3,scale="free")+
    scale_x_discrete(breaks = seq(min(bio_gg$Years),max(bio_gg$Years),50))+
    scale_color_brewer(palette = "Paired")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "gray"))
 
  Years<-rep(1:((Tmax-1)-BI),times=as.numeric(length(table(Biomass[,ncol(Biomass)]))))
  Flows<-cbind(Years,Flows)
  f_gg<-gather(Flows,key="Flows",value="Biomass",2:(ncol(Flows)-1))
 
  gg_flow<-ggplot(f_gg,aes(x=as.factor(Years),y=Biomass,group=Simulation))+
    geom_line(aes(color=Flows))+
    facet_wrap(.~factor(Flows,levels = Out.file$Data$flows),ncol=3,scale="free")+
    scale_x_discrete(breaks = seq(min(f_gg$Years),max(f_gg$Years),50))+
    # scale_color_brewer(palette = "Paired")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "gray"))
 
  Years<-rep(1:((Tmax-1)-BI),times=as.numeric(length(table(Biomass[,ncol(Biomass)]))))
  Phi<-cbind(Years,Phi)
  Phi<-Phi[,c("Years","Pelagics","Demersals","Simulation")]
  Phi_gg<-gather(Phi,key="F.Mortality",value="Biomass",2:(ncol(Phi)-1))
 
  gg_phi<-ggplot(Phi_gg,aes(x=as.factor(Years),y=Biomass,group=Simulation))+
    geom_line(aes(color=F.Mortality))+
    facet_wrap(.~factor(F.Mortality,levels = Out.file$Data$Species[5:6]),ncol=1,scale="free")+
    scale_x_discrete(breaks = seq(min(Phi_gg$Years),max(Phi_gg$Years),50))+
    # scale_color_brewer(palette = "Paired")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "gray"))
 
  gg<-list(Biomass=gg_bio,Flows=gg_flow,Phi=gg_phi)
 
  return(gg)
}

# cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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

# Output_BI$Output.BI$Biomass<-cbind(Output_BI$Output.BI$Biomass,rep(1:(NDNDData$Tmax-AutoCor$Burn_In),times=as.numeric(length(table(Output_BI$Output.BI$Biomass[,ncol(Output_BI$Output.BI$Biomass)])))))
