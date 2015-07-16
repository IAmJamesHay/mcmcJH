//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <iostream>
#include <fstream>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace std;
using namespace Rcpp;


//[[Rcpp::export]]
int randNumber(const int n){return floor(unif_rand()*n);}

//[[Rcpp::export]]
int test(NumericVector myfixed){
  vector<int> fixed;
  
  for(int i=0;i<myfixed.length();++i){
    if(myfixed[i] == 0){
      fixed.push_back(i);
    }    
  }

  for(int i = 0; i <fixed.size();++i){
    cout << fixed[i] << endl;
  }
  int jj = randNumber(fixed.size());
  return(fixed[jj]);
}
  
