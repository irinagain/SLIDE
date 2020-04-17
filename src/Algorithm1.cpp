// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void updateVoptimC(const arma::mat &XU, double lambda, const arma::colvec &pvec, arma::mat& V) {
  int d = pvec.n_elem;
  int r = XU.n_cols;
  //arma::mat V = arma::zeros(XU.n_rows, r);
  //arma::colvec pcum = arma::zeros(d+1);
  //pcum(arma::span(1, pcum.end())) = cumsum(pvec);
  V.zeros();
  arma::colvec a = arma::zeros(1);
  arma::colvec pcum = arma::join_cols(a, cumsum(pvec));
  arma::rowvec norms(r);
  
  for (int i=0; i<d; i++ ){
    // Index
    // arma::colvec index = arma::span(pcum[i] + 1, pcum[i+1]);
    arma::uvec index = arma::regspace<arma::uvec>(pcum[i], pcum[i+1] - 1);
    
    // Norms of each column
    norms = arma::sqrt(arma::sum(arma::square(XU(arma::span(pcum[i], pcum[i+1] - 1), arma::span::all)), 0));

    // Soft-thresholding
    arma::uvec indexr = arma::find(norms > lambda);
    if (indexr.n_elem > 1) {
     V(index, indexr) = XU(index, indexr) * arma::diagmat(1 - lambda/norms(indexr));
    } else if (indexr.n_elem == 1) {
      V(index, indexr) = XU(index, indexr) * (1 - lambda/norms(indexr));
    }
  }	                     
  //return V;
}


// [[Rcpp::export]]
double evaluatef_optimC(const arma::mat& XU, const arma::mat& V, double lambda, const arma::colvec& pvec) {
  double f = accu(square(V))/2 - sum(diagvec(V.t() * XU));
  int d = pvec.n_elem;
  arma::colvec a = arma::zeros(1);
  arma::colvec pcum = arma::join_cols(a, cumsum(pvec));
  for (int i=0; i<d; i++){
    // Index
    arma::uvec index = arma::regspace<arma::uvec>(pcum[i], pcum[i+1] - 1);
    // Adjust objective
    f += lambda * as_scalar(arma::sum(arma::sqrt(arma::sum(arma::square(V.rows(index)), 0)),1));
  }
  return(f);
}

// [[Rcpp::export]]
Rcpp::List solve_optimC(const arma::mat& X, arma::mat& U, double lambda, const arma::colvec& pvec, int k_max = 1000, double eps = 1e-06){
  arma::mat V = X.t() * U;
  // Current function value
  arma::colvec f(k_max);
  arma::colvec error(k_max);
  int k = 0;
  
  arma::mat Q, R;
  arma::vec s;
  
  f[k] = evaluatef_optimC(V, V, lambda, pvec);
  error[k] = 100;
  while ((k < k_max) & (error[k] > eps)) {
    k++;
    // Update V
    updateVoptimC(X.t() * U, lambda, pvec, V);
    // Update U
    arma::uvec nonzero = find(sum(abs(V),0) > 0);
    if (nonzero.n_elem == 0) {
      // V is exactly zero, terminate
      f[k] = 0;
      error[k] = f[k - 1] - f[k];
      // Return the list of values
      return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V, Rcpp::Named("k") = k, Rcpp::Named("error") = error, Rcpp::Named("f") = f, Rcpp::Named("fmin") = f[k]);
    } else {
      arma::svd_econ(Q, s, R, X * V);
      U = Q * R.t();
    }
    
    // Current function value
    f[k] = evaluatef_optimC(X.t() * U, V, lambda, pvec);
      
    // Current difference in function values
    error[k] = f[k - 1] - f[k];
  }
  // Return the list of values
  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V, Rcpp::Named("k") = k, Rcpp::Named("error") = error, Rcpp::Named("f") = f, Rcpp::Named("fmin") = f[k]);
}

// [[Rcpp::export]]
Rcpp::List solve_optim_seqC(const arma::mat& X, const arma::colvec& pvec, const arma::colvec& lambda_seq, int k_max = 1000, double eps = 1e-06) {
 
  int n_lambda = lambda_seq.n_elem;
  // Generate Ustart
  arma::mat Ustart, V;
  arma::vec s;
  arma::svd_econ(Ustart, s, V, X, "left");
  
  // Solve for each lambda value
  Rcpp::List param(n_lambda);
  Rcpp::List out;
  for (int l = n_lambda - 1; l >= 0; l--){
    // Use neighboring U as a new starting point
    //if (l < n_lambda - 1) {
    //  Ustart = param[l + 1]["U"];
    //}
    // Solve group lasso problem
    out = solve_optimC(X, Ustart, lambda_seq[l], pvec, k_max, eps);
    double fmin = out["fmin"];
    param[l] = Rcpp::List::create(Rcpp::Named("U") = out["U"], Rcpp::Named("V") = out["V"], Rcpp::Named("fmin") = fmin);
  }
  return(Rcpp::List::create(Rcpp::Named("param") = param, Rcpp::Named("lambda") = lambda_seq));
}