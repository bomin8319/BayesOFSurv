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
                          arma::vec eXB,
                          arma::vec alpha,
                          arma::vec C,
                          double lambda) {
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
  arma::vec dexp2 = pow(eXB % Y, lambda - 1);
  arma::vec dexp3 = pow(eXB % Y, lambda);
  arma::vec llik1 = log((1 - alpha) % dexp1 + lambda * alpha % eXB % dexp2 % dexp1);
  arma::uvec ids1 = find(llik1 == arma::datum::inf);
  llik1.elem(ids1).fill(exp(700));
  arma::uvec ids2 = find(dexp3 == arma::datum::inf);
  dexp3.elem(ids2).fill(exp(700));
  arma::vec llik = C % llik1 + (1 - C) % (- dexp3);
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
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
  arma::vec dexp2 = pow(eXB % Y, lambda - 1);
  arma::vec dexp3 = pow(eXB % Y, lambda);
  arma::vec llik1 = log((1 - alpha) % dexp1 + lambda * alpha % eXB % dexp2 % dexp1);
  arma::uvec ids1 = find(llik1 == arma::datum::inf);
  llik1.elem(ids1).fill(exp(700));
  arma::uvec ids2 = find(dexp3 == arma::datum::inf);
  dexp3.elem(ids2).fill(exp(700));
  arma::vec llik = C % llik1 + (1 - C) % (- dexp3);
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
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
  arma::vec dexp2 = pow(eXB % Y, lambda - 1);
  arma::vec dexp3 = pow(eXB % Y, lambda);
  arma::vec llik1 = log((1 - alpha) % dexp1 + lambda * alpha % eXB % dexp2 % dexp1);
  arma::uvec ids1 = find(llik1 == arma::datum::inf);
  llik1.elem(ids1).fill(exp(700));
  arma::uvec ids2 = find(dexp3 == arma::datum::inf);
  dexp3.elem(ids2).fill(exp(700));
  arma::vec llik = C % llik1 + (1 - C) % (- dexp3);
  return sum(llik);
}


// **********************************************************//
//     	            Likelihood function for beta             //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull_betas2 (arma::vec Y,
                          arma::vec eXB,
                          arma::vec alpha,
                          arma::vec C,
                          double lambda) {
	arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
	arma::vec dexp2 = pow(eXB % Y, lambda - 1);
	arma::vec dexp3 = pow(eXB % Y, lambda);
	arma::vec llik1 = log((1 - alpha) + lambda * alpha % eXB % dexp2 % dexp1);
	arma::uvec ids1 = find(llik1 == arma::datum::inf);
	llik1.elem(ids1).fill(exp(700));
	arma::uvec ids2 = find(dexp3 == arma::datum::inf);
	dexp3.elem(ids2).fill(exp(700));
	arma::vec llik = C % llik1 + (1 - C) % (log(alpha) - dexp3);
	return sum(llik);
}

// **********************************************************//
//     	            Likelihood function for gamma            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull_gammas2 (arma::vec Y,
                           arma::vec eXB,
                           arma::mat Z,
                           arma::vec gammas,
                           arma::vec C,
                           double lambda) {
  arma::vec ZG = Z * gammas;
	arma::vec alpha = 1 / (1 + exp(-ZG));
	arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
	arma::vec dexp2 = pow(eXB % Y, lambda - 1);
	arma::vec dexp3 = pow(eXB % Y, lambda);
	arma::vec llik1 = log((1 - alpha) + lambda * alpha % eXB % dexp2 % dexp1);
	arma::uvec ids1 = find(llik1 == arma::datum::inf);
	llik1.elem(ids1).fill(exp(700));
	arma::uvec ids2 = find(dexp3 == arma::datum::inf);
	dexp3.elem(ids2).fill(exp(700));
	arma::vec llik = C % llik1 + (1 - C) % (log(alpha) - dexp3);
	return sum(llik);
}

// **********************************************************//
//     	           Likelihood function for lambda            //
// **********************************************************//
// [[Rcpp::export]]
double llikWeibull_lambda2 (arma::vec Y,
					arma::vec eXB, 
					arma::vec alpha,
					arma::vec C,
					double lambda) {
  arma::vec dexp1 = exp(-pow(eXB % Y, lambda));
	arma::vec dexp2 = pow(eXB % Y, lambda - 1);
	arma::vec dexp3 = pow(eXB % Y, lambda);
	arma::vec llik1 = log((1 - alpha) + lambda * alpha % eXB % dexp2 % dexp1);
	arma::uvec ids1 = find(llik1 == arma::datum::inf);
	llik1.elem(ids1).fill(exp(700));
	arma::uvec ids2 = find(dexp3 == arma::datum::inf);
	dexp3.elem(ids2).fill(exp(700));
	arma::vec llik = C % llik1 + (1 - C) % (log(alpha) - dexp3);
	return sum(llik);
}

