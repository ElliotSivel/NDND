###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Output preparation
#### Version v1.0
#### 03.12.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################

source(file = "C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_code/0_1_NDND_initialization.R")

# 1. Load data ------------------------------------------------------------

load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/NDNDsim_out.Rdata")

# 2. Estimate the Burn-in -------------------------------------------------

# To estimate the burn-in, we need to estimate the Autocorrelation of the single time series

AC<-NDND.auto.cor(Output) # Computes autocorrelation functions for each single simulation and estimates the burn-in period

save(AC,file=paste(dir_tab$outputs_dir,"/AC/AC_",format(Output$Sim.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
save(AC,file=paste(dir_tab$outputs_dir,"/AC/ACfile.RData",sep = ""))
load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/AC/ACfile.Rdata")

x11()
AC$Correlogram # Plot autocorrelation

ggsave("C:/Users/a22073/Desktop/Comparison NDND version/Updated/AC.tif", scale=1, width = 16, height = 12, units = "cm",dpi = 300,device = "tiff")

# Delete the Burn-in for the simulation outputs
O_BI<-NDND.BI(AutoCor = AC, Out.file = Output)

save(O_BI,file=paste(dir_tab$outputs_dir,"/Output_Burnin/NDNDSimBI_",format(O_BI$Sim.tag,"%Y_%m_%d_%H_%M_%S"),".RData",sep = ""))
save(O_BI,file=paste(dir_tab$outputs_dir,"/Output_Burnin/NDNDSimBI.RData",sep = ""))
