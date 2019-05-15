## code for first simulation 13.05.19

# Run for 200 years
# 8 Species
# 2000 simulations

T1<-Sys.time()

# 1. Reference run --------------------------------------------------------

# create the files where the times series are stocked 
bio_ref<-NULL
flow_ref<-NULL
phi_ref<-NULL

# simulations
for (s in 1:2000){
  print(s)
  
  ref_sim<-SimNDND(NDNDData)
  
  bio_ser<-cbind(ref_sim$Output$BiomassSeries)
  flow_ser<-cbind(ref_sim$Output$FlowSeries)
  phi_ser<-cbind(ref_sim$Output$PhiSeries)
  
  bio_ref<-abind(bio_ref,bio_ser,along=3)
  flow_ref<-abind(flow_ref,flow_ser,along=3)
  phi_ref<-abind(phi_ref,phi_ser,along = 3)
}

simref<-list(bio_ref,flow_ref,phi_ref)
save(simref,file=paste(dir_tab$outputs_dir,"/simref",Sys.Date(),".Rdata",sep=""))


# 2. Run for temperature increase by 1 degree ------------------------------

# create the files where the times series are stocked
bio_t1<-NULL
flow_t1<-NULL
phi_t1<-NULL

# Transform NDNDData
NDNDData_t1<-TEOMR(NDNDData,c(1,1,1,1,1,1,0,0),1)
save(NDNDData_t1,file =paste(dir_tab$data_dir,"/NDNDData_t1.Rdata",sep=""))

# simulations
for (s in 1:2000){
  print(s)
  
  t1_sim<-SimNDND(NDNDData_t1)
  
  bio_ser<-cbind(t1_sim$Output$BiomassSeries)
  flow_ser<-cbind(t1_sim$Output$FlowSeries)
  phi_ser<-cbind(t1_sim$Output$PhiSeries)
  
  bio_t1<-abind(bio_t1,bio_ser,along=3)
  flow_t1<-abind(flow_t1,flow_ser,along=3)
  phi_t1<-abind(phi_t1,phi_ser,along = 3)
}

simt1<-list(bio_t1,flow_t1,phi_t1)
save(simt1,file=paste(dir_tab$outputs_dir,"/simt1",Sys.Date(),".Rdata",sep=""))


# 3. Run for temperature increase by 2 degrees ----------------------------

# create the files where the times series are stocked
bio_t2<-NULL
flow_t2<-NULL
phi_t2<-NULL

# Transform NDNDData
NDNDData_t2<-TEOMR(NDNDData,c(1,1,1,1,1,1,0,0),2)
save(NDNDData_t2,file =paste(dir_tab$data_dir,"/NDNDData_t2.Rdata",sep=""))

# simulations
for (s in 1:2000){
  print(s)
  
  t2_sim<-SimNDND(NDNDData_t2)
  
  bio_ser<-cbind(t2_sim$Output$BiomassSeries)
  flow_ser<-cbind(t2_sim$Output$FlowSeries)
  phi_ser<-cbind(t2_sim$Output$PhiSeries)
  
  bio_t2<-abind(bio_t2,bio_ser,along=3)
  flow_t2<-abind(flow_t2,flow_ser,along=3)
  phi_t2<-abind(phi_t2,phi_ser,along = 3)
}

simt2<-list(bio_t2,flow_t2,phi_t2)
save(simt2,file=paste(dir_tab$outputs_dir,"/simt2",Sys.Date(),".Rdata",sep=""))

# 4. Run for reference temperature and less fishing mortality -------------

# BSECO project set a scenario to 25% less fishing mortality
# Implement variation in fishing mortality
NDNDData_Fm<-NDNDData
NDNDData_Fm$Fmp<-NDNDData_Fm$Fmp*0.75
save(NDNDData_Fm,file =paste(dir_tab$data_dir,"/NDNDData_Fm.Rdata",sep=""))

# create the files where the times series are stocked 
bio_Fm<-NULL
flow_Fm<-NULL
phi_Fm<-NULL

# simulations
for (s in 1:2000){
  print(s)
  
  Fm_sim<-SimNDND(NDNDData_Fm)
  
  bio_ser<-cbind(Fm_sim$Output$BiomassSeries)
  flow_ser<-cbind(Fm_sim$Output$FlowSeries)
  phi_ser<-cbind(Fm_sim$Output$PhiSeries)
  
  bio_Fm<-abind(bio_Fm,bio_ser,along=3)
  flow_Fm<-abind(flow_Fm,flow_ser,along=3)
  phi_Fm<-abind(phi_Fm,phi_ser,along = 3)
}

simFm<-list(bio_Fm,flow_Fm,phi_Fm)
save(simFm,file=paste(dir_tab$outputs_dir,"/simFM",Sys.Date(),".Rdata",sep=""))

# 5. Run for less fishing mortality and temperature increase by 1 degree  --------

# create the files where the times series are stocked
bio_Fmt1<-NULL
flow_Fmt1<-NULL
phi_Fmt1<-NULL

# Transform NDNDData
NDNDData_Fmt1<-TEOMR(NDNDData_Fm,c(1,1,1,1,1,1,0,0),1)
save(NDNDData_Fmt1,file =paste(dir_tab$data_dir,"/NDNDData_Fmt1.Rdata",sep=""))

# simulations
for (s in 1:2000){
  print(s)
  
  Fmt1_sim<-SimNDND(NDNDData_Fmt1)
  
  bio_ser<-cbind(Fmt1_sim$Output$BiomassSeries)
  flow_ser<-cbind(Fmt1_sim$Output$FlowSeries)
  phi_ser<-cbind(Fmt1_sim$Output$PhiSeries)
  
  bio_Fmt1<-abind(bio_Fmt1,bio_ser,along=3)
  flow_Fmt1<-abind(flow_Fmt1,flow_ser,along=3)
  phi_Fmt1<-abind(phi_Fmt1,phi_ser,along = 3)
}

simFmt1<-list(bio_Fmt1,flow_Fmt1,phi_Fmt1)
save(simFmt1,file=paste(dir_tab$outputs_dir,"/simFmt1",Sys.Date(),".Rdata",sep=""))

# 6. Run for less fishing mortality and temperature increase by 2 degrees --------

# create the files where the times series are stocked
bio_Fmt2<-NULL
flow_Fmt2<-NULL
phi_Fmt2<-NULL

# Transform NDNDData
NDNDData_Fmt2<-TEOMR(NDNDData,c(1,1,1,1,1,1,0,0),2)
save(NDNDData_Fmt2,file =paste(dir_tab$data_dir,"/NDNDData_Fmt2.Rdata",sep=""))

# simulations
for (s in 1:2000){
  print(s)
  
  Fmt2_sim<-SimNDND(NDNDData_Fmt2)
  
  bio_ser<-cbind(Fmt2_sim$Output$BiomassSeries)
  flow_ser<-cbind(Fmt2_sim$Output$FlowSeries)
  phi_ser<-cbind(Fmt2_sim$Output$PhiSeries)
  
  bio_Fmt2<-abind(bio_Fmt2,bio_ser,along=3)
  flow_Fmt2<-abind(flow_Fmt2,flow_ser,along=3)
  phi_Fmt2<-abind(phi_Fmt2,phi_ser,along = 3)
}

simFmt2<-list(bio_Fmt2,flow_Fmt2,phi_Fmt2)
save(simFmt2,file=paste(dir_tab$outputs_dir,"/simFmt2",Sys.Date(),".Rdata",sep=""))
# 7. Run for reference temperature and more fishing mortality -------------

# BSECO project set a scenario to 50% more fishing mortality
# Implement variation in fishing mortality
NDNDData_Fp<-NDNDData
NDNDData_Fp$Fmp<-NDNDData_Fp$Fmp*1.5
save(NDNDData_Fp,file =paste(dir_tab$data_dir,"/NDNDData_Fp.Rdata",sep=""))

# create the files where the times series are stocked 
bio_Fp<-NULL
flow_Fp<-NULL
phi_Fp<-NULL

# simulations
for (s in 1:2000){
  print(s)
  
  Fp_sim<-SimNDND(NDNDData_Fp)
  
  bio_ser<-cbind(Fp_sim$Output$BiomassSeries)
  flow_ser<-cbind(Fp_sim$Output$FlowSeries)
  phi_ser<-cbind(Fp_sim$Output$PhiSeries)
  
  bio_Fp<-abind(bio_Fp,bio_ser,along=3)
  flow_Fp<-abind(flow_Fp,flow_ser,along=3)
  phi_Fp<-abind(phi_Fp,phi_ser,along = 3)
}

simFp<-list(bio_Fp,flow_Fp,phi_Fp)
save(simFp,file=paste(dir_tab$outputs_dir,"/simFp",Sys.Date(),".Rdata",sep=""))

# 8. Run for more fishing mortality and temperature increase by 1 degree  --------

# create the files where the times series are stocked
bio_Fpt1<-NULL
flow_Fpt1<-NULL
phi_Fpt1<-NULL

# Transform NDNDData
NDNDData_Fpt1<-TEOMR(NDNDData_Fp,c(1,1,1,1,1,1,0,0),1)
save(NDNDData_Fpt1,file =paste(dir_tab$data_dir,"/NDNDData_Fpt1.Rdata",sep=""))

# simulations
for (s in 1:2000){
  print(s)
  
  Fpt1_sim<-SimNDND(NDNDData_Fpt1)
  
  bio_ser<-cbind(Fpt1_sim$Output$BiomassSeries)
  flow_ser<-cbind(Fpt1_sim$Output$FlowSeries)
  phi_ser<-cbind(Fpt1_sim$Output$PhiSeries)
  
  bio_Fpt1<-abind(bio_Fpt1,bio_ser,along=3)
  flow_Fpt1<-abind(flow_Fpt1,flow_ser,along=3)
  phi_Fpt1<-abind(phi_Fpt1,phi_ser,along = 3)
}

simFpt1<-list(bio_Fpt1,flow_Fpt1,phi_Fpt1)
save(simFpt1,file=paste(dir_tab$outputs_dir,"/simFpt1",Sys.Date(),".Rdata",sep=""))

# 9. Run for more fishing mortality and temperature increase by 2 degrees --------

# create the files where the times series are stocked
bio_Fpt2<-NULL
flow_Fpt2<-NULL
phi_Fpt2<-NULL

# Transform NDNDData
NDNDData_Fpt2<-TEOMR(NDNDData_Fp,c(1,1,1,1,1,1,0,0),2)
save(NDNDData_Fpt2,file =paste(dir_tab$data_dir,"/NDNDData_Fpt2.Rdata",sep=""))

# simulations
for (s in 1:2000){
  print(s)
  
  Fpt2_sim<-SimNDND(NDNDData_Fpt2)
  
  bio_ser<-cbind(Fpt2_sim$Output$BiomassSeries)
  flow_ser<-cbind(Fpt2_sim$Output$FlowSeries)
  phi_ser<-cbind(Fpt2_sim$Output$PhiSeries)
  
  bio_Fpt2<-abind(bio_Fpt2,bio_ser,along=3)
  flow_Fpt2<-abind(flow_Fpt2,flow_ser,along=3)
  phi_Fpt2<-abind(phi_Fpt2,phi_ser,along = 3)
}

simFpt2<-list(bio_Fpt2,flow_Fpt2,phi_Fpt2)
save(simFpt2,file=paste(dir_tab$outputs_dir,"/simFpt2",Sys.Date(),".Rdata",sep=""))

T2<-Sys.time()

Tdif<-difftime(T2,T1)
Tdif







#####
# # Avec plot
# par(mfrow=c(3,3))
# for (s in 1:length(NDNDData$Species)) {
#   bio_ref[,,s]
#   plot(1:NDNDData$Tmax,bio_ref[,s,1],type="l",ylim = c(0,120),col="gray")
#   xlab=NDNDData$Species[s]
#   for (p in 2:2000){
#     lines(1:NDNDData$Tmax,bio_ref[,s,p],col="gray")
#   }
# }

# Avec ggplot
library(plyr)

f_data<-adply(bio_ref,c(1,3))
f_data<-gather(f_data,key="Species",value="Biomass",3:10)
colnames(f_data)<-c("Year","Simulation","Species","Biomass")

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gg<-ggplot(f_data,aes(x=Year,y=Biomass,group=Simulation))+
  geom_line(aes(color=Species))+
  facet_wrap(.~factor(Species,levels = ref_sim$Data$Species),ncol=3,scale="free")+
  scale_x_discrete(breaks = seq(0,200,25))+
  scale_color_brewer(palette = "Paired")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "gray"))

x11()
gg

ggsave("C:/Users/a22073/Documents/Model_NDND/NDND/ndnd_outputs/TM_biomass_ref_S1.pdf")
