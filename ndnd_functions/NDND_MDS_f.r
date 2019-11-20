###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Multi-Dimensional Scaling function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# This function performs a MDS on the NDND outputs

NDND.MDS<-function(Out.file,kmeans,method){
  
  Biomass<-Out.file$Output$Biomass
  Biomass_0<-drop_na(as.data.frame(Biomass))
  Biomass<-log10(Biomass_0[which(apply(Biomass_0[,1:8],1,sum)>0),1:8])
  
  a<-as.numeric(Sys.time())
  set.seed(a)

  k<-kmeans(Biomass[,1:8],kmeans)
  
  Bio_k<-cbind(Biomass,Cluster=k$cluster)
  
  d<-dist(Biomass[,1:8],method = method)
  
  mds2d<-isoMDS(d,k=2)                          
                                                  
  mds2d$points<-as.data.frame(mds2d$points)
  
  mds_plot_2<-ggplot(mds2d$points,aes(mds2d$points[,1],mds2d$points[,2]))+
    geom_jitter(aes(color=as.factor(k$cluster)))+
    stat_ellipse(aes(fill=as.factor(k$cluster),color=as.factor(k$cluster)),type="norm",size=1)
  
  mds3d<-isoMDS(d,k=3)
  
  mds3d$points<-as.data.frame(mds3d$points)
  
  mds_plot_3<-plot3d(x=mds3d$points[,1],y=mds3d$points[,2],z=mds3d$points[,3], col=as.factor(k$cluster))
  
  
  output<-list(Sim.tag=Out.file$Sim.tag,
               Partition=Bio_k,
               MDS2=mds2d,
               MDS3=mds3d,
               MDS.plot2d=mds_plot_2,
               MDS.plot3d=mds_plot_3)
  
  return(output)
}