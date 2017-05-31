load('Thyne.rdata')
summary(Thyne)

library(mvtnorm)
library(MCMCpack)
library(BayesOFsurv)

Z=as.matrix(cbind(1,Thyne$lforest,
                  Thyne$lngdppc_ipolate,Thyne$lnbdead))
X=as.matrix(cbind(1,Thyne$lforest,
                  Thyne$lngdppc_ipolate,Thyne$lnbdead))
Y=as.matrix(Thyne$`_t`)
C=as.matrix(Thyne$`_d`)

######### Try fit the data using Bayesian OF model #############  
Weibull = mcmcOF(Y, C, X, Z, N = 10000, burn = 5000, thin = 20,  w = c(1, 1, 1), m = 100, form = "Weibull")
summary(mcmc(Weibull$beta))
summary(mcmc(Weibull$gamma))
summary(mcmc(Weibull$lambda))
par(mfrow = c(2,3))
for (p in 1:2) {
  plot(Weibull$beta[,p], type = 'l')
}
for (p in 1:3) {
  plot(Weibull$gamma[,p], type = 'l')
}
plot(Weibull$lambda, type = 'l')

Exponential = mcmcOF(Y, C, X, Z, N = 10000, burn = 5000, thin = 20,  w = c(1, 1, 1), m = 100, form = "Exponential")
summary(mcmc(Exponential$beta))
summary(mcmc(Exponential$gamma))
summary(mcmc(Exponential$lambda))
par(mfrow = c(2,3))
for (p in 1:2) {
  plot(Exponential$beta[,p], type = 'l')
}
for (p in 1:3) {
  plot(Exponential$gamma[,p], type = 'l')
}
plot(Exponential$lambda, type = 'l')
