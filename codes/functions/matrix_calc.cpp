// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B,
                      int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult3(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B,
                      Eigen::Map<Eigen::MatrixXd> C,
                      int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd D = A * B * C;
  return Rcpp::wrap(D);
}

