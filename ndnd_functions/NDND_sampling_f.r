###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### sampling function
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel
#### License ???
###############################################################################################

# Needed libraries
###
library(LIM)                                # Loading the LIM package -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)
library(linprog)                            # Loading the linprog package -- Linear programming and optimization
require(symengine)
require(ggplot2)
require(ggraph)
require(coda)
require(parallel)
require(RCaN)
###

##### Sampling function

Sampling<-function(pA,pb){
  G<-as.matrix(pA)                       # Defining A for the sampling
  h<-as.matrix(pb)                       # Defining b for the sampling
    
  x0<-chebycenter(G,h)                   # Need for a starting point x0 : chebycenter method to define it
    
  Fsample<-RCaN::cpgs(1000,G,h,x0)             # Sample with Gibbs algorithm 100 vectors of flows
    
  F0<-Fsample[sample(100:nrow(Fsample),1),]           # Select one vector of flow for the among the output vectors
  return(F0)
}