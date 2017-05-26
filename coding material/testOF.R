library(mvtnorm)
library(MCMCpack)
library(BayesOFsurv)

n = 1000
nsims = 1
#create covariates
x<-runif(n, min=-2.5, max=12)
z<-log(runif(n, min=1, max=100))
tru.est<-matrix(NA,nrow=nsims,ncol=8)
i = 1
  #Assign parameter values
  tru.est[i,1]<-1
  tru.est[i,2]<-3.5
  tru.est[i,3]<--2
  tru.est[i,4]<-2
  tru.est[i,5]<-3
  tru.est[i,6]<-1
  
  myrates <- exp(tru.est[i,1]+(tru.est[i,2]*x)) 
  y <- rexp(n, rate = myrates) # generates the r.v.
  cen <- rexp(n, rate = 1 )
  ycen <- pmin(y, cen)
  di <- as.numeric(y <= cen)
  tru.est[i,7]<-table(di)[1]
  
  #create parameters for ZG
  alpha<-1/(1+exp(-(tru.est[i,3]+tru.est[i,4]*z+tru.est[i,5]*x)))
  print(mean(alpha))
  yzero<-matrix(1,n,1)
  error<--1*rlogis(n)
  flag<-error<qlogis(alpha)
  yzero[flag]<-error[flag]
  flag<-yzero==1
  di[flag]<-ifelse(di[flag]==0,yzero[flag],di[flag])
  tru.est[i,8]<-table(di)[1]
  
  data<-cbind(ycen,di,x,z)
  
  data<-data
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  Z<-cbind(1,z,x)
  
######### Try fit the data using Bayesian OF model #############  
Weibull = mcmcOF(Y, C, X, Z, N = 10000, burn = 5000, thin = 20,  w = c(1, 1, 1), m = 10, form = "Weibull", seed = 100)
summary(mcmc(Weibull$beta))
summary(mcmc(Weibull$gamma))
summary(mcmc(Weibull$lambda))
par(mfrow = c(2,4))
plot(Weibull$loglike, type = 'l')
for (p in 1:2) {
  plot(Weibull$beta[,p], type = 'l')
}
for (p in 1:3) {
  plot(Weibull$gamma[,p], type = 'l')
}
plot(Weibull$lambda, type = 'l')


Exponential = mcmcOF(Y, C, X, Z, N = 10000, burn = 5000, thin = 10,  w = c(1, 2, 0.1), form = "Exponential", seed = 100)
summary(mcmc(Exponential$beta))
summary(mcmc(Exponential$gamma))
par(mfrow = c(2,3))
plot(Exponential $loglike, type = 'l')
for (p in 1:2) {
  plot(Exponential$beta[,p], type = 'l')
}
for (p in 1:3) {
  plot(Exponential $gamma[,p], type = 'l')
}

