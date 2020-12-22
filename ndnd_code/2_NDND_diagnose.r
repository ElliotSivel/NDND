###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Output Diagnose
#### Version v1.0
#### 03.12.19
#### Author : Elliot Sivel, Benjamin Planque and Ulf Lindstrøm
#### License ???
###############################################################################################


# 0. Libraries ------------------------------------------------------------

# Needed Libraries
###
require(tidyverse)
require(ggplot2)
###

# Source the initialization file, functions, directories, etc,...
source(file = "C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_code/0_1_NDND_initialization.R") # Used to load all functions in the work environment

x11()

# 1. Load outputs ---------------------------------------------------------

load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/NDNDsim_out.RData") # Load the latest simulation
# load("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/NDNDSim_2020_11_18_19_11_20.RData") # Load the chose simulation after a specific time tag

# 2. Check for HCR -----------------------------------

# We check if the caught biomass follows the harvest control rule we have defined
Fish<-cbind(Output$Fish[,apply(Output$Fish,2,sum)!=0],"Years"=rep(1:Output$Data$Tmax,times=s_length)) # We consider only the species for which we consider fisheries
Bio<-cbind(Output$Biomass[,apply(Output$Fish,2,sum)!=0],"Years"=rep(1:Output$Data$Tmax,times=s_length)) # We get the biomass time-series for those species

Fish.gg<-pivot_longer(Fish,c(1,2),values_to = "Landings", names_to = "Species") # Format the data for ggplot
Bio.gg<-pivot_longer(Bio,c(1,2),values_to = "Biomass", names_to = "Species")

Fish.gg<-cbind(Fish.gg,"Biomass"=Bio.gg$Biomass) # Merge landings and biomass in the same table 
Fish.gg$Ratio<-Fish.gg$Landings/Fish.gg$Biomass # HCR are computed as a ratio of landings per unit of biomass, we reconstruct this ratio
Fish.gg<-Fish.gg[Fish.gg$Years %in% 1:(Output$Data$Tmax-1),]

ggplot(data = Fish.gg,aes(x=Biomass,y=Ratio))+geom_point(aes(color=Species))+facet_wrap(.~factor(Species,levels = unique(Fish.gg$Species)),ncol=2,scales="free") # plot the HCRs
ggplot(data = Fish.gg,aes(x=Biomass,y=Landings))+geom_point(aes(color=Species))+facet_wrap(.~factor(Species,levels = unique(Fish.gg$Species)),ncol=2,scales="free") # plot the actual estimated landings

# 3. Master equation ------------------------------------------------------

SD<-Sim.Diagnose(Out.file = Output) # This step computes all separate block of the master equation of the NDND model, computes the differences between biomass at each time steps in the simulations and computes the difference between the added blocks and the variations in the simulations.
dif<-pivot_longer(SD$dif,cols= c(1:8),names_to = "Species", values_to = "Differences")
dif$Species<-factor(dif$Species,levels=Output$Data$Species)
ggplot(data = dif,aes(x=Species,y=Differences))+geom_boxplot(fill="turquoise")+theme_bw()
# A recurrent issue on herbivorous zooplankton is found


# 4. Inertia constraint respected -----------------------------------------

up_rho<-cbind((exp(Output$Data$Rho)*Output$Biomass[,-9]),"Simulation"=Output$Biomass[,9],"Years"=rep(1:Output$Data$Tmax,times=length(unique(Output$Biomass[,9]))))
low_rho<-cbind((exp(-Output$Data$Rho)*Output$Biomass[,-9]-Output$Fish[,-9]),"Simulation"=Output$Biomass[,9],"Years"=rep(1:Output$Data$Tmax,times=length(unique(Output$Biomass[,9]))))

test_inertia_up<-array(0,dim = c(Output$Data$Tmax-1,length(Output$Data$Species),length(unique(Output$Biomass[,ncol(Output$Biomass)]))))
test_inertia_low<-array(0,dim = c(Output$Data$Tmax-1,length(Output$Data$Species),length(unique(Output$Biomass[,ncol(Output$Biomass)]))))

for (s in 1:length(unique(Output$Biomass[,ncol(Output$Biomass)]))) {
  sub<-Output$Biomass[Output$Biomass[,ncol(Output$Biomass)] %in% unique(Output$Biomass[,ncol(Output$Biomass)])[s],]
  sub_rho_up<-up_rho[Output$Biomass[,ncol(Output$Biomass)] %in% unique(Output$Biomass[,ncol(Output$Biomass)])[s],]
  sub_rho_low<-low_rho[Output$Biomass[,ncol(Output$Biomass)] %in% unique(Output$Biomass[,ncol(Output$Biomass)])[s],]
  for (i in 1:length(Output$Data$Species)){
    for (j in 1:(Output$Data$Tmax-1)) {
      test_inertia_up[j,i,s]<- sub[j,i] < sub_rho_up[j+1,i]
      test_inertia_low[j,i,s]<- sub[j,i] > sub_rho_low[j+1,i]
    }
  }
}

test_inertia_up<-as.logical(test_inertia_up);test_inertia_low<-as.logical(test_inertia_low)
test_inertia_up[test_inertia_up==FALSE]
test_inertia_low[test_inertia_low==FALSE]

up_rho.gg<-pivot_longer(up_rho,c(1:8),names_to = "Species",values_to = "Rho")
low_rho.gg<-pivot_longer(low_rho,c(1:8),names_to = "Species",values_to = "-Rho")
Rho.gg<-cbind(up_rho.gg,"-Rho"=low_rho.gg$`-Rho`)
Rho.gg$Rho<-sqrt(sqrt(Rho.gg$Rho));Rho.gg$`-Rho`<-sqrt(sqrt(Rho.gg$`-Rho`))

ggplot(data = Rho.gg, aes(x=Years))+geom_ribbon(aes(ymin=`-Rho`,ymax=Rho),fill="black",alpha=0.25)+facet_wrap(.~factor(Species,levels = Output$Data$Species)+Simulation,ncol=3,scales = "free")


# 5. Satiation constraint respected ---------------------------------------

IMP<-NULL
for(i in 1:length(unique(Output$Biomass$Simulation))){
  IMP<-rbind(IMP,Output$Data$Importall[-1,])
}

gainsonly<-as.data.frame(cbind(SD$Gains[,-1]-IMP[,-1],"Years"=rep(1:(Output$Data$Tmax-1),times=length(unique(Output$Biomass$Simulation))), "Simulation"=rep(1:length(unique(Output$Biomass$Simulation)),byrow=F,each=Output$Data$Tmax-1)))
gg.gainsonly<-pivot_longer(gainsonly,c(1:7),names_to="Species",values_to="Gains")
gg.gainsonly$Gains<-log10(gg.gainsonly$Gains)

Sat<-cbind(Output$Biomass[,-c(1,9)]*Output$Data$Sgma[-1],"Years"=rep(1:(Output$Data$Tmax),times=length(unique(Output$Biomass$Simulation))), "Simulation"=rep(1:length(unique(Output$Biomass$Simulation)),byrow=F,each=Output$Data$Tmax))
Sat<-Sat[Sat$Years %in% 2:Output$Data$Tmax,]
Sat$Years<-Sat$Years-1
gg.sat<-pivot_longer(Sat,c(1:7),names_to = "Species",values_to = "Sat_limit")
gg.sat$Sat_limit<-log10(gg.sat$Sat_limit)

ggplot()+geom_line(data = gg.sat,aes(x=Years,y=Sat_limit,group=Simulation),color="blue")+geom_line(data = gg.gainsonly,aes(x=Years,y=Gains,group=Simulation),color="darkgrey")+facet_wrap(.~factor(Species,levels = Output$Data$Species)+Simulation,ncol=3,scales="free")


# 6. Refuge biomas respected ----------------------------------------------

bta_m<-Output$Data$Bta %>% rep(times=length(unique(Output$Biomass$Simulation))*Output$Data$Tmax) %>% matrix(ncol=8,byrow = T) %>% cbind("Years"=rep(1:(Output$Data$Tmax),times=length(unique(Output$Biomass$Simulation))), "Simulation"=rep(1:length(unique(Output$Biomass$Simulation)),byrow=F,each=Output$Data$Tmax)) %>% as.data.frame()
colnames(bta_m)<-c(Output$Data$Species,"Years","Simulation")
gg.bta<-pivot_longer(bta_m,c(1:8),names_to="Species",values_to = "Beta")

bio<-cbind(Output$Biomass,"Years"=rep(1:(Output$Data$Tmax),times=length(unique(Output$Biomass$Simulation))))
gg.bio<-pivot_longer(bio,c(1:8),names_to = "Species",values_to = "Biomass")

ggplot()+geom_line(data = gg.bta, aes(x=Years,y=Beta,group=Simulation),color="blue",linetype="dashed")+geom_line(data=gg.bio,aes(x=Years,y=Biomass,group=Simulation),color="darkgrey")+facet_wrap(.~factor(Species,levels = Output$Data$Species),ncol=3,scales="free")


# 7. Test inequality ------------------------------------------------------

for (s in 1:length(unique(Output$Biomass$Simulation))){
  print(s)
  for (i in 1:(Output$Data$Tmax-1)){
    Atest<-Output$A[,,s]
    btest<-Output$b[,i,s]
    Ftest<-as.matrix(Output$Flows[i,-19])
    Test<-Atest%*%t(Ftest)<=btest
    if (sum(Test)!=nrow(Atest)){
      print("Inequality false")
    } else {
      print("TRUE")
    }
  }
}

Output$A[,,2]
Output$b[,3,2]
Output$Flows[i,-19]
