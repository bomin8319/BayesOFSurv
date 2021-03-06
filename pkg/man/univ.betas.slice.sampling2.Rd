% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesOFsurv.R
\name{univ.betas.slice.sampling2}
\alias{univ.betas.slice.sampling2}
\title{univ.betas.slice.sampling2}
\usage{
univ.betas.slice.sampling2(betas.p, p, Sigma.b, Y, X, betas, alpha, C, lambda,
  w, m, lower = -Inf, upper = +Inf, form)
}
\arguments{
\item{betas.p}{current value of the pth element of betas}

\item{p}{pth element}

\item{Sigma.b}{variance estimate of betas}

\item{Y}{response variable}

\item{X}{covariates for betas}

\item{betas}{current value of betas}

\item{alpha}{probability of true censoring}

\item{C}{censoring indicator}

\item{lambda}{current value of lambda}

\item{w}{size of the slice in the slice sampling}

\item{m}{limit on steps in the slice sampling}

\item{lower}{lower bound on support of the distribution}

\item{upper}{upper bound on support of the distribution}

\item{form}{type of parametric model (Exponential or Weibull)}
}
\value{
One sample update using slice sampling
}
\description{
univariate slice sampling for betas.p
}
