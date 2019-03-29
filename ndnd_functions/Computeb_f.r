###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Computeb function
#### Version v1.0
#### 28.03.19
#### Author : Elliot Sivel
###############################################################################################

##### Computeb function
# This function computes the right part of the inequalty
# Ax =< b
# Computes b, the vector of constraints applied on biomasses
# Need 3 vector of biomasses : Initial biomasses, Import values and Export values
# Need 4 vector pf parameters : Rho (Inertia), Sigma (Satiation), Beta (Refuge Biomass) and Mu (Metabolic losses)
# Need 1 initialization values : nn (Number of possible flows)

Computeb<-function(Biomass,Import,Export,Rho,Sigma,Bta,Mu,nn){
  
  # Set identification flag
  
  ID<-NULL
  ID[1]<-as.character(Sys.time())                           # Print Date and Time
  ID[2]<-R.Version()$version.string                         # Print R version
  ID[3]<-Sys.info()[1]                                      # Print the Computer software version (Windows, Mac or Linux)
  ID[4]<-Sys.info()[2]
  ID[5]<-Sys.info()[3]
  ID[6]<-Sys.info()[5]
  
  # Compute the C and D vectors : Simplification of elements of the master equation
  
  C=as.matrix((1-exp(-Mu))/Mu)                              # Computes C
  
  D=as.matrix(exp(-Mu))                                     # Computes D
  
  # Implementing constraints on biomasses
  # There are 5 constraints -- b1,b2,b3,b4,b5
  
  # First constraint : Biomasses are bounded below - Refuge Biomass
  
  b1<-as.matrix((1/C)*(D*Biomass-Bta)+Import-Export)                  # Computes b1
  
  # Second constraint : Biomass increases are bounded above - inertia (Rho)
 
  b2<-as.matrix(((1/C)*(exp(Rho)-D)*Biomass)+Export)                  # Computes b2

  # Third constraint : Flows are bounded above - satiation (Sigma)
  
  b3<-as.matrix(Biomass*Sigma)                                        # Computes b3

  # Fourth constraint, flows are positive
  
  b4<-matrix(0,nn,1)                                                  # Computes b4
  
  # Fifth constraint, Biomass decreases are bounded below - Inertia (-Rho)
  
  b5<-as.matrix(((1/C)*(D-exp(-Rho))*Biomass)+Import)                 # Computes b5
  
  # Create an object where all constraints on biomasses are present --> b
  
  b<-rbind(b1,b2,b3,b4,b5)
  
  return(b)
}
