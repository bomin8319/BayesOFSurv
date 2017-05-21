#' @useDynLib BayesOFsurv
#' @importFrom stats
#' @import grDevices
#' @importFrom graphics 
#' @importFrom MCMCpack riwish
#' @importFrom mvtnorm rmvnorm 
NULL

#' @title beta.slice.sampling
#' @description slice sampling for beta
#'
#' @param beta previous value of beta used in the slice sampling
#' @param Sigma.b variance estimate of beta
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for beta
#' @param phi current value of phi
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#'
#' @return One sample update using slice sampling
#'
#' @export
beta.slice.sampling = function(beta, Sigma.b, Y, C, X, phi, lambda, w, m = 100) {
  p1 = length(beta)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    beta[p] = univ.beta.slice.sampling(beta[p], p, beta, Sigma.b, Y, C, X, phi, lambda, w, m)
  }
  return(beta)
}

#' @title univ.beta.slice.sampling
#' @description univariate slice sampling for beta.p
#'
#' @param beta.p current value of the pth element of beta
#' @param p pth element
#' @param beta current value of beta
#' @param Sigma.b variance estimate of beta
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for beta
#' @param phi current value of phi
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.beta.slice.sampling = function(beta.p, p, beta, Sigma.b, Y, C, X, phi, lambda, w, m = 100, lower = -Inf, upper = +Inf) {
  b0 = beta.p
  b.post0 = beta.post(b0, p, beta, Sigma.b, Y, C, X, phi, lambda)
  
  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat 
    { if (L <= lower) break
      if (beta.post(L, p, beta, Sigma.b, Y, C, X, phi, lambda) <= b.post0) break
      L = L - w
    }
    repeat 
    {
      if (R >= upper) break
      if (beta.post(R, p, beta, Sigma.b, Y, C, X, phi, lambda) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (beta.post(L, p, beta, Sigma.b, Y, C, X, phi, lambda) <= b.post0) break
      L = L - w
      J = J - 1
    }
  
    while (K > 0) {
      if (R >= upper) break
      if (beta.post(R, p, beta, Sigma.b, Y, C, X, phi, lambda) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat 
  {
    b1 = runif(1, L, R)
    b.post1 = beta.post(b1, p, beta, Sigma.b, Y, C, X, phi, lambda)
    
    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
  }

#' @title gamma.slice.sampling
#' @description slice sampling for gamma
#'
#' @param gamma previous value of gamma used in the slice sampling
#' @param Sigma.g variance estimate of gamma
#' @param Y response variable
#' @param C censoring indicator
#' @param XB current value of XB
#' @param Z covariates for gamma
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#'
#' @return One sample update using slice sampling
#'
#' @export
gamma.slice.sampling = function(gamma, Sigma.g, Y, C, XB, Z, lambda, w, m = 100) {
  p2 = length(gamma)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gamma[p] = univ.gamma.slice.sampling(gamma[p], p, gamma, Sigma.g, Y, C, XB, Z, lambda, w, m)
  }
  return(gamma)
}

#' @title univ.gamma.slice.sampling
#' @description univariate slice sampling for gamma.p
#'
#' @param gamma.p current value of the pth element of gamma
#' @param p pth element
#' @param gamma current value of gamma
#' @param Sigma.g variance estimate of gamma
#' @param Y response variable
#' @param C censoring indicator
#' @param XB current value of XB
#' @param Z covariates for gamma
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#'
#' @return One sample update using slice sampling
#'
#' @export
univ.gamma.slice.sampling = function(gamma.p, p, gamma, Sigma.g, Y, C, XB, Z, lambda, w, m = 100, lower = -Inf, upper = +Inf) {
  g0 = gamma.p
  g.post0 = gamma.post(g0, p, gamma, Sigma.g, Y, C, XB, Z, lambda)
  
  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat 
    { if (L <= lower) break
      if (gamma.post(L, p, gamma, Sigma.g, Y, C, XB, Z, lambda) <= g.post0) break
      L = L - w
    }
    repeat 
    {
      if (R >= upper) break
      if (gamma.post(R, p, gamma, Sigma.g, Y, C, XB, Z, lambda) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (gamma.post(L, p, gamma, Sigma.g, Y, C, XB, Z, lambda) <= g.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (gamma.post(R, p, gamma, Sigma.g, Y, C, XB, Z, lambda) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat 
  {
    g1 = runif(1, L, R)
    g.post1 = gamma.post(g1, p, gamma, Sigma.g, Y, C, XB, Z, lambda)
    
    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

#' @title lambda.slice.sampling
#' @description univariate slice sampling for lambda
#'
#' @param lambda current value of lambda
#' @param Y response variable
#' @param C censoring indicator
#' @param XB current value of XB
#' @param phi current value of phi
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#'
#' @return One sample update using slice sampling
#'
#' @export
lambda.slice.sampling = function(lambda, Y, C, XB, phi, w, m = 100, lower = 0 + 10^(-10), upper = +Inf) {
  l0 = lambda
  l.post0 = lambda.post(l0, Y, C, XB, phi)
  
  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat 
    { if (L <= lower) break
      if (lambda.post(L, Y, C, XB, phi) <= l.post0) break
      L = L - w
    }
    repeat 
    {
      if (R >= upper) break
      if (lambda.post(R, Y, C, XB, phi) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J
    
    while (J > 0) {
      if (L <= lower) break
      if (lambda.post(L, Y, C, XB, phi) <= l.post0) break
      L = L - w
      J = J - 1
    }
    
    while (K > 0) {
      if (R >= upper) break
      if (lambda.post(R, Y, C, XB, phi) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }
  
  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }
  
  repeat 
  {
    l1 = runif(1, L, R)
    l.post1 = lambda.post(l1, Y, C, XB, phi)
    
    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  names(l1) = l.post1
  return(l1)
}


#' @title beta.post
#' @description log-posterior distribution of beta with pth element fixed as beta.p 
#'
#' @param beta.p current value of the pth element of beta
#' @param p pth element
#' @param beta current value of beta
#' @param Sigma.b variance estimate of beta
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for beta
#' @param phi current value of phi
#' @param lambda current value of lambda
#'
#' @return log- posterior density of beta
#'
#' @export
beta.post = function(beta.p, p, beta, Sigma.b, Y, C, X, phi, lambda) {
  beta[p] = beta.p
  XB = X %*% beta
  llik = C * (log((1 - phi) + phi * exp(XB) * lambda * 
         ((exp(XB) * Y)^(lambda - 1)) * exp(-(exp(XB) * Y)^lambda))) + 
         (1 - C) * (log(phi) -(exp(XB) * Y)^lambda)
  lprior = dmvnorm(beta, rep(0, length(beta)), Sigma.b, log = TRUE)
  lpost = sum(llik) + lprior
  return(lpost)
}

#' @title gamma.post
#' @description log-posterior distribution of gamma with pth element fixed as gamma.p 
#'
#' @param gamma.p current value of the pth element of gamma
#' @param p pth element
#' @param gamma current value of gamma
#' @param Sigma.g variance estimate of gamma
#' @param Y response variable
#' @param C censoring indicator
#' @param XB current value of XB
#' @param Z covariates for gamma
#' @param lambda current value of lambda
#'
#' @return log- posterior density of beta
#'
#' @export
gamma.post = function(gamma.p, p, gamma, Sigma.g, Y, C, XB, Z, lambda) {
  gamma[p] = gamma.p
  ZG = Z %*% gamma
  phi = 1 / (1 + exp(-ZG))
  llik = C * (log((1 - phi) + phi * exp(XB) * lambda * 
                    ((exp(XB) * Y)^(lambda - 1)) * exp(-(exp(XB) * Y)^lambda))) + 
    (1 - C) * (log(phi) -(exp(XB) * Y)^lambda)
  lprior = dmvnorm(gamma, rep(0, length(gamma)), Sigma.g, log = TRUE)
  lpost = sum(llik) + lprior
  return(lpost)
}

#' @title lambda.post
#' @description log-posterior distribution of lambda 
#'
#' @param lambda current value of lambda
#' @param Y response variable
#' @param C censoring indicator
#' @param XB current value of XB
#' @param phi current value of phi
#' @param a shape parameter of gamma prior
#' @param b scale parameter of gamma prior
#'
#' @return log- posterior density of beta
#'
#' @export
lambda.post = function(lambda, Y, C, XB, phi, a = 1, b = 1) {
  llik = C * (log((1 - phi) + phi * exp(XB) * lambda * 
         ((exp(XB) * Y)^(lambda - 1)) * exp(-(exp(XB) * Y)^lambda))) + 
         (1 - C) * (log(phi) -(exp(XB) * Y)^lambda)
  lprior = dgamma(lambda, a, b, log = TRUE)
  lpost = sum(llik) + lprior
  return(lpost)
}

#' @title mcmcOF
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian parametric OF model
#'
#' @param Y response variable
#' @param C censoring indicator
#' @param X covariates for beta
#' @param Z covariates for gamma
#' @param N number of MCMC iterations
#' @param burn burn-in to be discarded
#' @param thin thinning to prevent from autocorrelation
#' @param form type of parametric model (Exponential or Weibull)
#' @param w size of the slice in the slice sampling for (beta, gamma, lambda)
#'
#' @return chain of the variables of interest
#'
#' @export
mcmcOF<- function(Y, C, X, Z, N, burn, thin, w = c(1, 1, 1), form, seed = 1) {
  set.seed(seed)
  p1 = dim(X)[2]
  p2 = dim(Z)[2]

  # initial values
  Sigma.b = 10 * p1 * diag(p1)  # multiply 10 to ensure large enough variance in the early stages
  Sigma.g = 10 * p2 * diag(p2)  # multiply 10 to ensure large enough variance in the early stages
  beta = rep(0, p1)
  gamma = rep(0, p2)
  lambda = 1
  XB = X %*% beta
  ZG = Z %*% gamma
  phi = 1 / (1 + exp(-ZG))
  
  beta.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
  gamma.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p2)
  lambda.samp = rep(NA, (N - burn) / thin)
  loglike.samp = rep(NA, (N - burn) / thin)
  for (iter in 1:N) {
    if (iter %% 500 == 0) print(iter)
    beta = beta.slice.sampling(beta, Sigma.b, Y, C, X, phi, lambda, w[1])
    XB = X %*% beta
    gamma = gamma.slice.sampling(gamma, Sigma.g, Y, C, XB, Z, lambda, w[2])
    ZG = Z %*% gamma
    phi = 1 / (1 + exp(-ZG))
    if (form %in% "Weibull") {
    lambdasamp = lambda.slice.sampling(lambda, Y, C, XB, phi, w[3])
    lambda = lambdasamp
    loglike = names(lambdasamp)
    } else {
      loglike = sum(C * (log((1 - phi) + phi * exp(XB) * exp(-(exp(XB) * Y)))) + 
        (1 - C) * (log(phi) -(exp(XB) * Y)))
    }
    if (iter > burn) {
      Sigma.b = riwish(1 + p1, beta %*% t(beta) + p1 * diag(p1))
      Sigma.g = riwish(1 + p2, gamma %*% t(gamma) + p2 * diag(p2))
      if ((iter - burn) %% thin == 0) {
      beta.samp[(iter - burn) / thin, ] = beta
      gamma.samp[(iter - burn) / thin, ] = gamma
      lambda.samp[(iter - burn) / thin] = lambda
      loglike.samp[(iter - burn) / thin] = loglike
      }
    }
  }
  return(list(loglike = loglike.samp, beta = beta.samp, gamma = gamma.samp, lambda = lambda.samp))
}
