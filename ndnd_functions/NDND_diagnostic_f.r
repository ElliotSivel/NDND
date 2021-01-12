###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Create a diagnostic pdf with figures
#### Version v1.0
#### 05.01.2021
#### Author : Elliot Sivel
###############################################################################################

##### Simulation diagnose function
require(scales)

ndnddiagnose<-function(Output){
  
  sim<-sample(1:s_length,30)
  
  Biomass<-Output$Biomass[Output$Biomass$Simulation %in% sim,]
  Flows<-Output$Flows[Output$Flows$Simulation %in% sim,]
  Fish<-Output$Fish[Output$Fish$Simulation %in% sim,]
  
  # Check for fisheries and HCRs
  subfish<-cbind(Fish[,apply(Fish,2,sum)!=0],"Years"=rep(1:Output$Data$Tmax,times=length(sim))) # We consider only the species for which we consider fisheries
  subbiomass<-cbind(Biomass[,apply(Fish,2,sum)!=0],"Years"=rep(1:Output$Data$Tmax,times=length(sim))) # We get the biomass time-series for those species
  
  Fish.gg<-pivot_longer(subfish,c(1,2),values_to = "Landings", names_to = "Species") # Format the data for ggplot
  Bio.gg<-pivot_longer(subbiomass,c(1,2),values_to = "Biomass", names_to = "Species")
  
  Fish.gg<-cbind(Fish.gg,"Biomass"=Bio.gg$Biomass) # Merge landings and biomass in the same table 
  Fish.gg$Ratio<-Fish.gg$Landings/Fish.gg$Biomass # HCR are computed as a ratio of landings per unit of biomass, we reconstruct this ratio
  Fish.gg<-Fish.gg[Fish.gg$Years %in% 1:(Output$Data$Tmax-1),]
  
  print(ggplot(data = Fish.gg,aes(x=Biomass,y=Ratio))+geom_point(aes(color=Species))+facet_wrap(.~factor(Species,levels = unique(Fish.gg$Species)),ncol=2,scales="free")+ggtitle("Ratio of catch per biomass unit")) # plot the HCRs
  print(ggplot(data = Fish.gg,aes(x=Biomass,y=Landings))+geom_point(aes(color=Species))+facet_wrap(.~factor(Species,levels = unique(Fish.gg$Species)),ncol=2,scales="free")+ggtitle("Landings in biomass")) # plot the landings
  
  # Check for differences between master equation estimation and simulations
  SD<-Sim.Diagnose(Output,sim) # This step computes all separate block of the master equation of the NDND model, computes the differences between biomass at each time steps in the simulations and computes the difference between the added blocks and the variations in the simulations.
  dif<-pivot_longer(SD$dif,cols= c(1:8),names_to = "Species", values_to = "Differences")
  dif$Species<-factor(dif$Species,levels=Output$Data$Species)
  print(ggplot(data = dif,aes(x=Species,y=Differences))+geom_boxplot(fill="turquoise")+theme_bw()+scale_y_continuous(limits = c(-1,1))+ggtitle("Differences between master equation variations and simulations"))
  
  # Check growth
  s_b<-Biomass %>% split(as.factor(Biomass[,ncol(Biomass)]))
  sub_bio<-map(s_b,function(x){
    sub<-as.matrix(x)
    sub<-sub[1:(Output$Data$Tmax-1),]
  })
  sub_bio<-do.call(rbind.data.frame,sub_bio)
  Biomass$Years<-rep(1:(Output$Data$Tmax),times=length(sim))
  
  growth<-NULL
  for (s in 1:length(sim)){
    subbio<-Biomass[Biomass$Simulation %in% sim[s],-c(9,10)]
    for (i in 1:(Output$Data$Tmax-1)) {
      a<-subbio[i+1,]/subbio[i,]
      growth<-rbind(growth,a)
    }
  }
  growth$Simulation<-sub_bio$Simulation;growth$Years<-rep(1:(Output$Data$Tmax-1),times=length(sim))

  gg_growth<-pivot_longer(growth,-c(9,10),names_to = "Species", values_to = "Growth")
  gg_growth$upinertiaconstraint<-Output$Data$Rho %>% exp() %>% rep(times=(Output$Data$Tmax-1)*length(sim))
  gg_growth$lowinertiaconstraint<-Output$Data$Rho %>% rep(times=(Output$Data$Tmax-1)*length(sim)) ; gg_growth$lowinertiaconstraint<-exp(-gg_growth$lowinertiaconstraint)
  gg_growth$Biomass<-sub_bio[,-c(9,10)] %>% t() %>% as.vector()
  
  print(ggplot(data = gg_growth,aes(x=Years,y=Growth))+geom_point(color="turquoise")+geom_line(aes(y=upinertiaconstraint),linetype="dashed")+geom_line(aes(y=lowinertiaconstraint))+facet_wrap(.~factor(Species,levels=unique(gg_growth$Species)),ncol=2,scales="free")+theme_bw()+ggtitle("Growth per year (Inertia)"))
  # ggplot(data = gg_growth,aes(x=Biomass,y=Growth))+geom_point(color="turquoise")+geom_hline(yintercept=0,linetype="dashed")+facet_wrap(.~factor(Species,levels=unique(gg_growth$Species)),ncol=2,scales="free")+theme_bw()
  print(ggplot(data = gg_growth,aes(x=Biomass,y=Growth))+geom_point(color="turquoise")+geom_line(aes(y=upinertiaconstraint))+geom_line(aes(y=lowinertiaconstraint))+geom_hline(yintercept = 1,linetype="dashed",alpha=0.7)+facet_wrap(.~factor(Species,levels=unique(gg_growth$Species)),ncol=2,scales="free")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+theme_bw()+ggtitle("Growth (in log10 scale) per unit of biomass (Inertia)"))
  
  # Check for satiation
  gains<-NULL
  for (s in 1:length(Output$Data$Species)) {
    sub_s<-as.matrix(Flows[,grep(pattern=paste("* --> ",Output$Data$Species[s],sep=""),x=colnames(Output$Flows))])
    sum_sp<-as.numeric(apply(sub_s, 1, sum))
    gains<-cbind(gains,sum_sp)
  }
  gains<-cbind(gains,"Simulation"=rep(sim,byrow=F,each=(Output$Data$Tmax-1)),"Years"=rep(1:(Output$Data$Tmax-1),times=length(sim)))
  colnames(gains)<-c(Output$Data$Species,"Simulation","Years")
  
  gg.gains<-pivot_longer(as.data.frame(gains[,-1]),-c(8,9),names_to="Species",values_to="Gains")
  gg.gains$Satiation<-Output$Data$Sgma[-1] %>% rep(times=length(sim)*(Output$Data$Tmax-1))
  gg.gains$Biomass<-sub_bio[,-c(1,9,10)] %>% t() %>% as.vector()
  gg.gains$ratio<-gg.gains$Gains/gg.gains$Biomass
  
  print(ggplot(data = gg.gains, aes(x=Biomass,y=ratio))+geom_point(color="turquoise")+geom_line(aes(y=Satiation),linetype="dashed")+facet_wrap(.~factor(Species,levels=unique(gg.gains$Species)),ncol=2,scales="free")+theme_bw()+ggtitle("Inflow of biomass per unit of biomass (Satiation)"))
  
  # Check for distribution of flows
  Flows$Years<-rep(1:(Output$Data$Tmax-1),times=length(sim))
  gg.flows<-pivot_longer(Flows,-c(19,20),names_to = "Flows",values_to = "Values")
  
  print(ggplot(data = gg.flows)+geom_histogram(aes(x=Values),fill="darksalmon")+facet_wrap(.~factor(Flows,levels=unique(gg.flows$Flows)),ncol=3,scales="free")+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+theme_bw()+ggtitle("Distribution of flows (on log10 scale)"))
  
  # Check for distribution of biomasses
  gg.bio<-pivot_longer(Biomass,-c(9,10),names_to = "Species", values_to = "Biomass")
  
  print(ggplot(data = gg.bio)+geom_histogram(aes(x=Biomass),fill="darksalmon")+facet_wrap(.~factor(Species,levels=unique(gg.bio$Species)),ncol=2,scales="free")+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+theme_bw()+ggtitle("Distribution of biomass (on log10 scale)"))
  
  # Check for refuge biomass constraint
  gg.bio$refbio<-rep(Output$Data$Bta,times=Output$Data$Tmax*length(sim))
  
  print(ggplot(data = gg.bio,aes(x=Years,y=Biomass))+geom_point(color="turquoise")+geom_line(aes(y=refbio),linetype="dashed")+facet_wrap(.~factor(Species,levels=unique(gg.bio$Species)),ncol=2,scales="free")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+theme_bw()+ggtitle("Refuge biomass"))
  
  #Check for positive flows
  print(ggplot(data = gg.flows)+geom_point(aes(x=Years,y=Values),color="turquoise")+geom_hline(yintercept = 0,linetype="dashed")+facet_wrap(.~factor(Flows,levels=unique(gg.flows$Flows)),ncol=3,scales="free")+theme_bw()+ggtitle("Positive Flows"))

}

