#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using std::log;
using std::exp;
using std::max;
using std::abs;
using std::sqrt;
using std::pow;

using namespace Rcpp; 

// **********************************************************//
//     	            Likelihood function for beta             //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull_betas (arma::vec Y,
					arma::mat X, 
					arma::vec betas, 
					arma::vec alpha,
					arma::vec C,
					double lambda) {
	arma::vec XB = X * betas;
	arma::vec eXB = exp(XB);
	arma::vec llik = C % (log(alpha % exp(-pow(eXB % Y, lambda)) + lambda * (1 - alpha) % eXB % pow(eXB % Y, lambda - 1) % exp(-pow(eXB % Y, lambda)))) +
					(1 - C) % (log(alpha) - pow(eXB % Y, lambda));		
	return sum(llik);
}

// **********************************************************//
//     	            Likelihood function for gamma            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull_gammas (arma::vec Y,
					arma::vec eXB, 
					arma::mat Z,
					arma::vec gammas,
					arma::vec C,
					double lambda) {
	arma::vec ZG = Z * gammas;
	arma::vec alpha = 1 / (1 + exp(-ZG))	;				
	arma::vec llik = C % (log(alpha % exp(-pow(eXB % Y, lambda)) + lambda * (1 - alpha) % eXB % pow(eXB % Y, lambda - 1) % exp(-pow(eXB % Y, lambda)))) +
					(1 - C) % (log(alpha) - pow(eXB % Y, lambda));		
	return sum(llik);
}

// **********************************************************//
//     	           Likelihood function for lambda            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull_lambda (arma::vec Y,
					arma::vec eXB, 
					arma::vec alpha,
					arma::vec C,
					double lambda) {			
	arma::vec llik = C % (log(alpha % exp(-pow(eXB % Y, lambda)) + lambda * (1 - alpha) % eXB % pow(eXB % Y, lambda - 1) % exp(-pow(eXB % Y, lambda)))) +
					(1 - C) % (log(alpha) - pow(eXB % Y, lambda));		
	return sum(llik);
}