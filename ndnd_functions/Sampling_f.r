###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### sampling function
#### Version v1.1.5
#### 07.05.19
#### Author : Elliot Sivel & Co
#### License ???
###############################################################################################

Sampling<-function(pA,pb,samp_type){
  if (samp_type=="LimSolve"){

    G<-as.matrix(-pA)                       # Defining A (constraints on fluxes) for the sampling
    h<-as.matrix(-pb)                       # Defining b (constraints on biomasses) for th sampling

    Fsample<-xsample(G=G,H=h,iter=10,burninlength=25,outputlength=10,type="mirror")

    F0<-Fsample[[1]][sample(1:nrow(Fsample[[1]]),1),]
  }
  else if (samp_type=="cpgs"){

    G<-as.matrix(pA)                       # Defining A (constraints on fluxes) for the sampling
    h<-as.matrix(pb)

    Fsample<-cpgs(100,G,h)

    F0<-Fsample[sample(1:nrow(Fsample),1),]
  }
  else if (samp_type=="cpgs2"){

    G<-as.matrix(pA)                       # Defining A (constraints on fluxes) for the sampling
    h<-as.matrix(pb)

    x0<-chebycenter(G,h)

    Fsample<-cpgs2(100,G,h,x0)

    F0<-Fsample[sample(1:nrow(Fsample),1),]
  }
  return(F0)
}