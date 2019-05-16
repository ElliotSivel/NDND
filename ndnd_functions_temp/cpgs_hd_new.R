library(RcppEigen)
library(Rcpp)
library(inline)

library(tidyr)
library(ggplot2)







setwd("~/Bureau/cgps/")

  
#creation de la fonction CPP    
cpgsCPP='#include <Eigen/Cholesky>
#include <Rcpp.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<MatrixXd>A(as<Map<MatrixXd> >(AA));
const Map<VectorXd>b(as<Map<VectorXd> >(bb));
const Map<VectorXd>x0(as<Map<VectorXd> >(xx0));

int N=Rcpp::as<int>(NN); 
int p=A.cols();
int m=A.rows();
int runup=100;
int discard=100;

// Check input arguments
if (m < (p+1)){
  throw std::range_error("dimensions mismatch");
}
// Initialisation

MatrixXd X(N+runup+discard,p);
int n=0;
MatrixXd x(p,1);
MatrixXd y(p,1);
x=x0;

// Initialize variables for keeping track of sample mean, covariance
// and isotropic transform.
MatrixXd M(p,1);
M.setZero();
MatrixXd S2(p,p);
S2.setZero();

// outer products.
MatrixXd S(p,p);
MatrixXd S0(p,p);
S.setIdentity();

MatrixXd T1(p,p);
T1.setIdentity();
MatrixXd W(m,p);

W = A;
MatrixXd d(m,1);
MatrixXd delta0(p,1);
MatrixXd delta1(p,1);
MatrixXd z(m,1);
while (n < (N+runup+discard)){               //sampling loop
  y=x;
  NumericVector alea2=runif(p);
  // compute approximate stochastic transformation
  if (n == runup){
    T1=S.transpose().llt().matrixU().transpose();
    W = A*T1;
  }
  y=T1.inverse()*y;
  
  // choose p new components
  VectorXd e(p);
  for (int i=0;i<p;++i){
    //Find points where the line with the (p-1) components x_i
    //fixed intersects the bounding polytope.
    e.setZero();
    e(i)= 1;
    z = (W*e); //prevent any divisions by 0
    
    d=(b - W*y);
    d=d.cwiseQuotient(z); 
    double tmin=-9e9;
    double tmax=9e9;
    for (int j=0;j<m;++j){
      if (z(j)<0 && tmin<d(j)) tmin=d(j);
      if (z(j)>0 && tmax>d(j)) tmax=d(j);
    }
    y(i)+=(tmin+(tmax-tmin)*alea2(i));
  }
  x=T1*y;
  X.row(n)= x.col(0);
  ++n;
  // Incremental mean and covariance updates
  delta0 = x - M; // delta new point wrt old mean
  M+= delta0/(double)n;     // sample mean
  delta1= x - M;      // delta new point wrt new mean
  
  
  if (n > 1){
    S2 +=(n-1)/(double)(n*n)*(delta0*delta0.transpose())+(delta1*delta1.transpose());
    S0 = S;
    S = S2/(double)(n-1);           // sample covariance
  } else {
    S.setIdentity();
  }
}

return wrap(X.bottomRows(N));'


chebycenter <- function(A,b){
   n <- dim(A)[1]
   p <- dim(A)[2];
   an <- rowSums(A^2)^0.5
   A1 <- matrix(data=0,nrow = n,ncol = p+1)
   A1[,1:p] <- A;
   A1[,p+1] <- an;
   f <- matrix(data=0,nrow=p+1,ncol=1)
   f[p+1] <- -1;
   
   #d <- linprog(cc=f,A=A1,b=b);
   #x <- d$x[1:p];
   d <- lp(direction="min",objective.in = f,const.mat = A1,const.rhs = as.numeric(b),const.dir = rep("<=",n))
   x <- d$solution
   return(x[-p-1])
 }



#compilation (à ne faire qu'une seule fois)
cpgs2<-cxxfunction(signature(NN="integer",AA="matrix",bb="vector",xx0="vector"),cpgsCPP,plugin="RcppEigen",verbose=TRUE)
  
#jeu de test
A=read.table('pA1.txt')
b=read.table('pb1.txt')
f=read.table('f1.txt') # one solution (=vector of fluxes) provided by Matlab
# the problem to solve is A*f<=b
# In LIM, this is rephrased in G*x>=h
G=as.matrix(-A)
h=as.matrix(-b)
x0=as.matrix(f)


N=10000
ndiscard=nrunup=100

#test de la fonction 
debut=Sys.time()
res=cpgs2(N=N,-G,-h,x0)
end=Sys.time()
duree=end-debut
print(duree)
colnames(res)=c('F01','F02','F03','F04','F05','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16','F17','F18')
X2=gather(as.data.frame(res),key='trophic.flow')
ggplot(data=X2,aes(x=value))+facet_wrap(~trophic.flow,nrow=4,scales="free")+geom_histogram(aes(y=..density..), colour="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")
  
#test de l'ancienne fonction
debut=Sys.time()
res2=cpgs(N=N,A=-G,b=-h,x0=x0[,1])
end=Sys.time()
duree=end-debut
print(duree)
colnames(res2)=c('F01','F02','F03','F04','F05','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16','F17','F18')
X3=gather(as.data.frame(res2),key='trophic.flow')
x11()
ggplot(data=X3,aes(x=value))+facet_wrap(~trophic.flow,nrow=4,scales="free")+geom_histogram(aes(y=..density..), colour="black", fill="white")+geom_density(alpha=.2, fill="#FF6666")

  
  
