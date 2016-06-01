//temp.cpp
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
arma::mat multiply_test(NumericMatrix a, double b){
  arma::mat c(a.nrow(),a.ncol());
  c = as<arma::mat>(a)*b;
  return(c);
}

//[[Rcpp::export]]
arma::mat cov_test(NumericMatrix a){
  return(cov(as<arma::mat>(a)));
}

//[[Rcpp::export]]
NumericVector subset_test(NumericVector a, NumericVector b){
  return(a[b]);
}

//[[Rcpp::export]]
arma::mat subset_mat(arma::mat a, arma::uvec b){
  return(a.submat(b,b));
}

//[[Rcpp::export]]
arma::cube cube_test(){
arma::cube a(5,5,3);
return(a);
}


//[[Rcpp::export]]
NumericMatrix create_blocks(NumericVector a){
  int no_infections = (a.size()-2)/6;
  NumericMatrix blocks(no_infections+1,2);
  blocks(0,0) = 0;
  blocks(0,1) = 2;

  for(int i = 1; i < blocks.nrow();++i){
    blocks(i,0) = 3 + 6*(i-1);
    blocks(i,1) = 8 + 6*(i-1);
  }
  return(blocks);
}
  

//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}
//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}
//[[Rcpp::export]]
arma::cube test_func(NumericMatrix start_covmatrix,NumericVector scales, NumericVector startvalue){
  int no_infections = (startvalue.size()-2)/6;
  
  // Create matrix that stores the indices of the blocks (one per infection and one for stats parameters)
  NumericMatrix blocks(no_infections+1,2);
  arma::mat tmp(6,6);
  blocks = create_blocks(startvalue);
  
  for(int i = 0; i < blocks.nrow();++i){
    cout << blocks(i,0) << " " << blocks(i, 1) << endl;
  }

  arma::cube cov_matrix(6,6,no_infections);
  arma::cube cov_matrix_scaled(6,6,no_infections);

  for(int i = 0; i < no_infections;++i){
    cov_matrix.slice(i) = (as<arma::mat>(start_covmatrix)).submat(blocks(i+1,0),blocks(i+1,0),blocks(i+1,1),blocks(i+1,1));
    cov_matrix_scaled.slice(i) = cov_matrix.slice(i)*scales(i+1);
  }
  
  return(cov_matrix_scaled);
}

//[[Rcpp::export]]
arma::mat subset_test2(NumericMatrix a, int b, int c){
  return((as<arma::mat>(a)).submat(b,b,c,c));
}

//[[Rcpp::export]]
arma::uvec testingagain(NumericVector a, int b){
  NumericVector c(a.size());
  c = a-(b-1)*6 - 3;
  return(as<arma::uvec>(c));
}


