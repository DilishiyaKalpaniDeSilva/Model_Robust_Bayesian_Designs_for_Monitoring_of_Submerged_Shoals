#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat expitCpp(arma::mat x) {
  return 1-(1/(1+arma::exp(x)));
}

// [[Rcpp::export]]
arma::mat dbinomCpp(arma::vec x,  arma::mat prob) {
  int n = x.size();
  arma::vec com(n);
  for(int i = 0; i < n; ++i) {
    com(i) =  Rf_choose(20, x[i]);
  }
  arma::mat probmi = 1-prob;
  arma::mat pp = (pow(prob.each_col(),x)) % (pow(probmi.each_col(),20-x)) ;
  return sum((log( pp.each_col() % com)),0) ;
}
