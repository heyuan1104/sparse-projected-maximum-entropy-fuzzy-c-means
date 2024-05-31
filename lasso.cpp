// lasso.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;



// [[Rcpp::export]]
arma::vec Lasso(arma::vec Bv,arma::mat Z, int lambda, int Bvd, int j) {
  
  if(Bvd > 0 && abs(Bvd) > (lambda/2)) {
    Bv(j-1) = (Bvd-lambda/2) / sum(Z.col(j-1).t() * Z.col(j-1));
  }
  if(Bvd < 0 && abs(Bvd) > (lambda/2)) {
    Bv(j-1) = (Bvd+lambda/2) / sum(Z.col(j-1).t() * Z.col(j-1));
  }
  if(abs(Bvd)<=(lambda/2)) {
    Bv(j-1) = 0;
  }
  
  return(Bv);
}
// END