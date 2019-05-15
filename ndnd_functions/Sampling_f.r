###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### sampling function
#### Version v1.1.5
#### 07.05.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

##### Sampling function
# Three option of samplng implemented 
# LIM algorithm for sampling based on mirror algorithm given by Soetaert and van Oevelen : slow
# cpgs : Gibbs algorithm on R, coded by Hilaire Drouineau on C++ and translated in R by Benjamin Planque
# cpgs2 : Gibbs algorithm on C++, coded by Hilaire Droouineau : Fastest option

Sampling<-function(pA,pb,samp_type){
  if (samp_type=="LimSolve"){

    G<-as.matrix(-pA)                       # Defining A for the sampling
    h<-as.matrix(-pb)                       # Defining b for th sampling
    
    # Sampling
    # iter = number of sampled items, 
    # burninlength = number of samples before saving sampling, 
    # outputlength= length of the object containg the sampled flows 
    Fsample<-xsample(G=G,H=h,iter=10,burninlength=25,outputlength=10,type="mirror")      

    F0<-Fsample[[1]][sample(1:nrow(Fsample[[1]]),1),]          # Choose randomly one vector of flows in the output
  }
  else if (samp_type=="cpgs"){

    G<-as.matrix(pA)                       # Defining A for the sampling
    h<-as.matrix(pb)                       # Defining b for the sampling
    
    x0<-chebycenter(G,h)                   # Need for a starting point x0 : chebycenter method to define it
    
    Fsample<-cpgs(100,G,h,x0)                 # Sample with Gibbs algorithm 100 vectors of flows 

    F0<-Fsample[sample(1:nrow(Fsample),1),]           # Select one vector of flow for the among the output vectors
  }
  else if (samp_type=="cpgs2"){

    G<-as.matrix(pA)                       # Defining A for the sampling
    h<-as.matrix(pb)                       # Defining b for the sampling

    x0<-chebycenter(G,h)                   # Need for a starting point x0 : chebycenter method to define it

    Fsample<-cpgs2(100,G,h,x0)             # Sample with Gibbs algorithm 100 vectors of flows

    F0<-Fsample[sample(1:nrow(Fsample),1),]           # Select one vector of flow for the among the output vectors
  }
  return(F0)
}