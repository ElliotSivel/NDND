###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Boxplot function
#### Version v1.0
#### 07.10.19
#### Author : Elliot Sivel
###############################################################################################

##### Boxplot function
# We plot the boxplot for each trophospecies according to the clusters
# Phyotplankton   |H.zooplankton   |O.zooplankton   |Benthos        |...  
# c1  c2  c3  c4  |c1  c2  c3  c4  |c1  c2  c3  c4  |c1  c2  c3  c4 |...

NDND.Boxplot<-function(MDS.out){
  Part<-MDS.out$Partition
  g<-gather(Part,key = "Species",value="Biomass",-Cluster)
  # g$Biomass<-10^(g$Biomass)
  Sp<-unique(g$Species)
  
  bxp<-ggplot(g,aes(x=Cluster,y=Biomass,group=Cluster))+
    geom_boxplot(aes(fill=as.factor(Cluster)))+
    
    facet_wrap(.~factor(Species,levels = Sp),ncol = 3,scales='free')+
    scale_fill_manual(values = c('springgreen4','firebrick','darkblue','gold3'))+
    theme_bw()+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    labs(x="Clusters",y="log10(Biomass)",fill="Clusters")
  
  
  # for (i in 1:(ncol(Part)-1)) {
  #   e<-as.numeric(Part[,i])
  #   f<-as.numeric(Part[,ncol(Part)])
  #   bxp<-boxplot(e~f,xlab=NULL,ylab=NULL,main=names(Part)[i],col=c("red","orange","lightblue","purple"),cex.lab=1.5)
  # }
  # mtext("Cluster",1,outer = T,line = -2,cex = 1.25)
  # mtext("log10(Biomass)",2,outer = T,line=-1.75,cex = 1.25)
  out_bxp<-list(Sim.tag=MDS.out$Sim.tag,
                bxp=bxp)
  return(out_bxp)
}
