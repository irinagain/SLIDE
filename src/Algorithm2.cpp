// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List slide_givenS_C(const arma::mat& X, const arma::vec& pvec, const arma::mat& S, arma::mat& U, double eps = 1e-06, int k_max = 1000) {
  
  arma::colvec a = arma::zeros(1);
  arma::colvec pcum = arma::join_cols(a, cumsum(pvec));
  int r = S.n_cols;
  int p = X.n_cols;
  int n = X.n_rows;
  int d = pvec.n_elem;
  arma::mat UVold(n, p);
  
  arma::mat Q, R;
  arma::vec s;
  
  arma::colvec error(k_max);
  
  // Initialize V
  arma::mat V = arma::zeros(p, r);
  for (int i=0; i<d; i++) {
    // Index the measurements corresponding to ith dataset
    arma::uvec index = arma::regspace<arma::uvec>(pcum[i], pcum[i+1] - 1);
    // Identify columns in V that are present for the dataset i
    arma::uvec nonzero = arma::find(S.row(i) == 1);
    if (nonzero.n_elem > 0) {
      V(index, nonzero) = X.cols(index).t() * U.cols(nonzero);
    }
  }

  int k = 0;
  error[k] = 1000;
  while ((k < k_max) & (error[k] > eps)) {
    k++;
    UVold = U * V.t();
     
    // Refit U
    arma::svd(Q, s, R, X * V);
    U = Q.cols(0, r-1) * R.t();
       
    // Refit V
    for (int i=0; i<d; i++) {
      // Index the measurements corresponding to ith dataset
      arma::uvec index = arma::regspace<arma::uvec>(pcum[i], pcum[i+1] - 1);
      // Identify columns in V that are present for the dataset i
      arma::uvec nonzero = arma::find(S.row(i) == 1);
      if (nonzero.n_elem > 0) {
       V(index, nonzero) = X.cols(index).t() * U.cols(nonzero);
      }
    }
    
     
    // Calculate the difference due to refitting
    error[k] = accu(square(U * V.t() - UVold));
  }
  
  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V, Rcpp::Named("k") = k, Rcpp::Named("error") = error);
}
