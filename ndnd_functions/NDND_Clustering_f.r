###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Clustering function
#### Version v1.1
#### 16.07.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# This function performs the clustering of the NDND outputs

NDND.cluster<-function(Out.file,D,Linkage){
  
  Sp<-names(Out.file$Partition[,1:length(NDNDData$Species)])
  Biomass<-Out.file$Partition                           # Extract the array of Biomass containing all simulations (Complete and incomplete)
  
  dend1<- Biomass[,1:length(Sp)] %>% 
    dist(method=D) %>% 
    hclust(method = Linkage) %>% 
    as.dendrogram
  
  dend1 %>%
    hang.dendrogram(hang = -1) %>%
    set("labels",NULL) %>%
    plot(xlab = "Configurations",ylab="Height",ylim = c(0,7))
  
  colored_bars(colors=Biomass$Cluster,dend = dend1,sort_by_labels_order = TRUE,rowLabels = "Cluster")
  #   
  # g_dend<-ggdendrogram(dend,labs=NULL)
  
  # df2<-data.frame(cluster=cutree(dend,4),kmeans=factor(dend$labels,levels=dend$labels))

  d_dend<-dendro_data(dend,type="rectangle")
  
  
  gg<-ggplot(segment(d_dend))+ 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    theme_bw() +
    ylab("Height")+
    xlab("Configurations")+
    geom_hline(yintercept=4.75,col="red",lty=2,lwd=2)+
    geom_text(x=0,y=5,label="k=4",cex=5,col="red")+
    geom_tile(Out.file)+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_line(),
          legend.position = "None") 
  
  # Opt_clust<-NbClust(c, diss=d, distance = NULL,method = "kmeans")       # Define the optimal number of cluster using the kmeans method
  # Tab_clust<-as.data.frame(table(Opt_clust$Best.nc[1,]))                 # Extracting the result of cluster given by all methods 
  
  # Gather the outputs in a list
  out_clust<-list(Sim.tag=Output$Sim.tag,
                  # log.trans=c,
                  # Dist.Matrix=d3,
                  # Cluster=e3,
                  # Dendrogram.stats.Gower=g3,
                  # Optimal.Clustering=Opt_clust,
                  # Cluster.Estimation=Tab_clust,
                  # Dendrogram=dend,
                  # Dendrogram_gg=dend_gg,
                  Dendrogram=gg)
  
  return(out_clust)
}