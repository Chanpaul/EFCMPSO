//#include <Rcpp.h>
//#include <stdlib.h>
#include <RcppArmadillo.h>
//#include <map>

//[[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace arma;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
arma::mat test(NumericMatrix x){
  
  return as<arma::mat>(x);
}

//[[Rcpp::export]]
NumericMatrix rcppdist(NumericMatrix dt) {
  int nrow=dt.nrow();
  //int ncol=dt.ncol();
  int h=0;
  NumericMatrix tdist(ceil(nrow*(nrow+1)/2),3);
  arma::mat matdt=as<arma::mat>(dt);
  for (int i=0;i<nrow-1;i++){
    for(int j=i+1;j<nrow ;j++){
      double temp=sum(pow(matdt.row(i)-matdt.row(j),2));
      
      tdist(h,0)=i+1;
      tdist(h,1)=j+1;
      tdist(h,2)=sqrt(temp);
      h=h+1;
    }
  }
  return(tdist);
}


// [[Rcpp::export]]
arma::vec cppComRho(arma::mat distMx,int dc){
  int nd= max(distMx.col(1));
  arma::vec distvec=sort(distMx.col(2));
  double dc_th=distvec(dc*0.01*distMx.n_rows)+0.00001;
  arma::vec rho=zeros(nd)+0.00000001;
  for (int i=0;i<distMx.n_rows;i++){
    int k=distMx(i,0)-1;
    int j=distMx(i,1)-1;
    double kj=exp(-1*pow(distMx(i,2)/dc_th,2));
    rho(k)=rho(k)+kj;
    rho(j)=rho(j)+kj;
  }
  return rho;
}
// [[Rcpp::export]]
arma::uvec testsort(arma::vec x){
  return sort_index(x);  
}

// [[Rcpp::export]]
arma::mat dist2mat(arma::mat x){
  int nr=max(x.col(0));
  arma::mat ydist(nr,nr,fill::zeros);
  for (int i=0;i<x.n_rows;i++){
    int k=x(i,0);
    int h=x(i,1);
    double v=x(i,2);
    ydist(k,h)=v;
    ydist(h,k)=v;
  }
  return ydist;
}

// [[Rcpp::export]]
arma::mat cppComDelta(arma::mat inDistMx,arma::vec rho){
  int nd=inDistMx.n_rows;
  arma::mat delta(nd,2,fill::zeros); //1st column: nearest neighbor; 2nd column distance;
  arma::uvec sort_rho_id=sort_index(rho,"descend");   //sort the density in descending order;
  //nneig(sort_rho_id(0))=sort_rho_id(0);
  delta(sort_rho_id(0),0)=sort_rho_id(0);
  
  arma::mat distMx=dist2mat(inDistMx);
  for (int i=1;i<nd;i++){
    int k=sort_rho_id(i);
    double nearDist=delta(k,1);
    int nearNeig=delta(k,0);
    for (int j=0;j<i;j++){
      int t=sort_rho_id(j);
      double tdist=distMx(k,t);
      if (tdist<nearDist){
        nearNeig=t;
        nearDist=tdist;
      }
    }
    delta(k,0)=nearNeig;
    delta(k,1)=nearDist;
  }
  delta(sort_rho_id(0),1)=max(delta.col(1));
  return delta;
}
