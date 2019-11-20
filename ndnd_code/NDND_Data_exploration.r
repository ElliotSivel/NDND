# 1. Autocorrelation analysis -----------------------------------------------------

AutoCor<-NDND.auto.cor(Out.file = Output)
x11()
AutoCor$Correlogram
ggsave(paste(wd,"/images/Autocor_plot.png",sep = ""), width =16, height = 12, units = "cm", pointsize=12)

# 2. Delete Burn-In -------------------------------------------------------

Output_BI<-NDND.BI(AutoCor = AutoCor, Out.file = Output)

# 3. Plot time series -----------------------------------------------------

Output_TS<-plot.NDND(Out.file = Output)
x11()
Output_TS$Biomass
ggsave(paste(wd,"/images/Timeseries_plot.png",sep = ""), width = 16, height = 12, units = "cm", pointsize=12)

# 4. Estimation of Fishing mortality --------------------------------------

Output_FM<-Phi.Mortality(Out.file = Output_BI)

# 5. Multi-Dimensional Scaling --------------------------------------------

Output_MDS<-NDND.MDS(Out.file = Output,kmeans = 4,method = "euclidean")
x11()
Output_MDS$MDS.plot2d

# 6. Perform Clustering ---------------------------------------------------

# NDNDCluster<-NDND.cluster(Out.file = Output_BI,D="euclidean",Linkage = "complete")

# 7. Boxplots -------------------------------------------------------------

Output_BXP<-NDND.Boxplot(MDS.out = Output_MDS4)
ggsave(paste(wd,"/images/Boxplot_plot.png",sep = ""), width = 16, height = 16, units = "cm", pointsize=12)
# 6. Adding Partition to Biomass matrix -----------------------------------

Output_Cluster<-NDND.out.partition(Data = NDNDData, Output = Output_BI, Cluster = NDNDCluster)

# 7. Defining System Configuration ----------------------------------------

Ecosys_Conf<-NDND.Sys.Config(Out.Cluster = Output_Cluster)

# 9. Estimating the number of shifts in the simulation --------------------

Nb_Shifts<-NDND.Nb.Shifts(Output = Output_BI, Out.Cluster = Output_Cluster)
