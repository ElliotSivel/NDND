NDNDSimulation$Output$BiomassSeries
delta_bio<-NULL
C<-(1-exp(-(Mu+Phi))/(Mu+Phi))
for(i in 2:nrow(NDNDSimulation$Output$BiomassSeries)){
  db<-NDNDSimulation$Output$BiomassSeries[i,]-NDNDSimulation$Output$BiomassSeries[i-1,]
  delta_bio<-rbind(delta_bio,db)
}
