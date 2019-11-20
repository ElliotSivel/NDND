###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Ecosystem configuration function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# This function computes the average and median profiles of the different clusters
# It performs pie charts of the distribution of the trophospecies for each cluster
NDND.Sys.Config<-function(Out.Cluster){
  
  # Create item containing the biomass for each cluster
  a1<-Out.Cluster$Biomass.cluster[Out.Cluster$Biomass.cluster[,9]==1,1:8]
  a2<-Out.Cluster$Biomass.cluster[Out.Cluster$Biomass.cluster[,9]==2,1:8]
  
  # Compute the mean biomass value for each trophospecies for the different clusters
  b1<-apply(a1,2, mean, na.rm=T)
  b2<-apply(a2,2, mean, na.rm=T)
  
  # Compute the proportions of the trophospecies in each clusters; for the mean values
  pb1<-b1*100/sum(b1)
  pb2<-b2*100/sum(b2)
  
  # Compute the median biomass value for each trophospecies for the different clusters
  c1<-apply(a1,2, median, na.rm=T)
  c2<-apply(a2,2, median, na.rm=T)
  
  # Compute the proportions of the trophospecies in each clusters; for the median values
  pc1<-c1*100/sum(c1)
  pc2<-c2*100/sum(c2)
  
  # Matrix manipulation for plotting -- Mean profiles
  d1<-as.data.frame(t(pb1))                                        # Transform matrix into dataframe
  d2<-gather(d1,key = "Species",value = "Prop")                    # Transform vector into a 2-columns dataframe : Proportions and Species
  
  d3<-as.data.frame(t(pb2))                                        # Transform matrix into dataframe
  d4<-gather(d3,key = "Species",value = "Prop")                    # Transform vector into a 2-columns dataframe : Proportions and Species

  # Matrix manipulation for plotting -- Median profiles
  e1<-as.data.frame(t(pc1))                                        # Transform matrix into dataframe
  e2<-gather(e1,key = "Species",value = "Prop")                    # Transform vector into a 2-columns dataframe : Proportions and Species
  
  e3<-as.data.frame(t(pc2))                                        # Transform matrix into dataframe
  e4<-gather(e3,key = "Species",value = "Prop")                    # Transform vector into a 2-columns dataframe : Proportions and Species
  
  # Plotting pie charts of mean system profiles
  gg<-ggplot(d2,aes(x="",y=Prop,fill=Species))+
    geom_bar(width=1,stat = "identity")+
    coord_polar("y",start = 0)+
    theme(axis.text.x = element_blank())+
    geom_text(aes(y=Prop/8+c(0,cumsum(Prop)[-length(Prop)]),label=percent(Prop/100)),size=5)
  
  gg1<-ggplot(d4,aes(x="",y=Prop,fill=Species))+
    geom_bar(width=1,stat = "identity")+
    coord_polar("y",start = 0)+
    theme(axis.text.x = element_blank())+
    geom_text(aes(y=Prop/8+c(0,cumsum(Prop)[-length(Prop)]),label=percent(Prop/100)),size=5)
  
  # Plotting pie charts of median system profiles
  gg2<-ggplot(e2,aes(x="",y=Prop,fill=Species))+
    geom_bar(width=1,stat = "identity")+
    coord_polar("y",start = 0)+
    theme(axis.text.x = element_blank())+
    geom_text(aes(y=Prop/8+c(0,cumsum(Prop)[-length(Prop)]),label=percent(Prop/100)),size=5)
  
  gg3<-ggplot(e4,aes(x="",y=Prop,fill=Species))+
    geom_bar(width=1,stat = "identity")+
    coord_polar("y",start = 0)+
    theme(axis.text.x = element_blank())+
    geom_text(aes(y=Prop/8+c(0,cumsum(Prop)[-length(Prop)]),label=percent(Prop/100)),size=5)
  
  output<-list(Sim.tag=Out.Cluster$Sim.tag,
               Cluster.1=list(Mean.Profile=b1,
                              Median.Profile=c1,
                              Prop.Mean.Profile=pb1,
                              Prop.Median.Profile=pc1,
                              Mean.Profile.Plot=gg,
                              Median.Profile.Plot=gg2),
               Cluster.1=list(Mean.Profile=b2,
                              Median.Profile=c2,
                              Prop.Mean.Profile=pb2,
                              Prop.Median.Profile=pc2,
                              Mean.Profile.Plot=gg1,
                              Median.Profile.Plot=gg3))
  
  return(output)
}
