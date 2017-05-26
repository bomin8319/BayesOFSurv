#include <RcppArmadillo.h>
#include <cmath>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
//[[Rcpp::depends(RcppArmadillo)]]
using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;

using namespace Rcpp; 

void R_init_markovchain(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);	
}

// **********************************************************//
//     	 Likelihood function for parametric OF - Weibull     //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull (arma::vec Y,
					arma::mat X, 
					arma::vec betas, 
					arma::mat Z, 
					arma::vec gammas,
					arma::vec C,
					double lambda) {
	arma::vec XB = X * betas;
	arma::vec eXB = exp(XB);
	arma::vec ZG = Z * gammas;
	arma::vec alpha = 1 / (1 + exp(-ZG));
	arma::vec llik = C % (log(alpha % exp(-pow(eXB % Y, lambda)) + lambda * (1 - alpha) % eXB % pow(eXB % Y, lambda - 1) % exp(-pow(eXB % Y, lambda)))) +
					(1 - C) % (log(alpha) - pow(eXB % Y, lambda));		
	return sum(llik);
}