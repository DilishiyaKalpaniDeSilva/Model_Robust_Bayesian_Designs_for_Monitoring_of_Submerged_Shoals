RcppFunctions <- function(){
  library(Rcpp)
  library(RcppArmadillo)
  library(RcppEigen)
  library(RcppParallel)
  sourceCpp("codes/functions/cpp_functions.cpp")
  sourceCpp("codes/functions/matrix_calc.cpp")
}


