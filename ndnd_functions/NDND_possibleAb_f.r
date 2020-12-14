###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### possibleAb function
#### Version v1.0
#### 24.06.19
#### Author : Elliot Sivel
###############################################################################################

##### possibleAb function
# Matrices built in ComputeA and Computeb are taking in account all possible flows
# Not all possible flows are real
# Flow values needed ony for existing flows
# possibleAb delete the elements that are not existing from the constraints matrices
# Need matrices A and b computed with ComputeA and Computeb
# Need the vector of possible flows PFv

possibleAb<-function(A,b,PFv){
  
  # Step 1 : Delete the not existing flows from the matrix of constraints on flows -- A, then the constraints that are equal to 0  
  
  Ap<-A[,which(PFv==1)]                               # Keep the columns for which PFv = 1 in a matrix -- Ap
  if (ncol(Ap)>1) {                                   # If ncol(Ap) > 1
    Alines<-apply(abs(Ap), 1,sum)>0                   # Sum up the matrix to keep the rows
                                                      # Identify the values that are bigger than 0
    Ap<-Ap[which(Alines==TRUE),]                      # Delete the lines that are = 0 
  } else if (ncol(Ap)==1) {                           # If ncol(Ap) = 1 
    Alines<-Ap>0                                      # Test if the values in the remaining vector are = 0 
    Ap<-Ap[which(Alines==TRUE)]                       # Delete the values of the vector that are = 0 
  } else {stop("no topology given")}                  # If ncol(Ap) = 0, no topology
  
  # The simulation stops and an error message appears
  # Step 2 : Select the constraints on biomasses for which we have constraints on flows

  bp<-as.matrix(b[which(Alines==T)])                  # Select elements in the vector b that are left in matrix Ap
  
  return(list(Ap,bp))
}