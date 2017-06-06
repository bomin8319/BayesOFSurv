###
###1. Store True values for X0, X1, Z0, Z1, Z2, P
###2. Store proportion censored pre & post
###3. Simmulate n of 1000, do this 1000 times
###4. Estimate cox, weibull, store all relevant coefficient estimates (exponentiate p's where applicable)
###5. For each value in 4, calculate CPs and RMSEs, store.
###6. Estimate zombie exp and zombie weibull-->store all relevant coefficient estimates (exponentiate p's where applicable).
###7. For each value in 6, calculate CPs and RMSEs, store.
###8. Estimate Bayesian zombie exp and Bayesian zombie weibull-->store all relevant coefficient estimates (exponentiate p's where applicable).
###9. For each value in 8, calculate CPs and RMSEs, store.
##############
####Set Up####
##############

#clear memory
rm( list=ls() )

#load necessary libraries 						                                 
library(foreign)
library(Zelig)
library(car)
library(MASS)
library(VGAM)
library(plotrix)
library(pscl)
library(survival)
library(msm)
library(verification)
library(corpcor)
library(Design)
library(mvtnorm)
library(MCMCpack)
#library(devtools)
#install_github('bomin8319/BayesOFsurv/pkg')
library(BayesOFsurv)

#set working directory
setwd("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Mixture DGP/")

##########################################################################
##########################################################################
############################Monte Carlo###################################
##########################################################################

#set seed
set.seed(3)   

#set the number of observations
n<-100

#set the number of simulations, and create matrices to store the results
nsims<-1000

#history matrix for true estimates
tru.est<-matrix(NA,nrow=nsims,ncol=8)
#history matrix for cox estimates
cox.est<-matrix(NA,nrow=nsims,ncol=2)
#history matrix for exp estimates
exp.est<-matrix(NA,nrow=nsims,ncol=24)
#history matrix for weibull estimates
weib.est<-matrix(NA,nrow=nsims,ncol=30)
#history matrix for cox RMSE
cox.rmse<-matrix(NA,nrow=nsims,ncol=1)
#history matrix for exp RMSE
exp.rmse<-matrix(NA,nrow=nsims,ncol=12)
#history matrix for exp RMSE
weib.rmse<-matrix(NA,nrow=nsims,ncol=15)
#history matrix for cox CP
cox.cp<-matrix(NA,nrow=nsims,ncol=1)
#history matrix for exp CP
exp.cp<-matrix(NA,nrow=nsims,ncol=12)
#history matrix for exp CP
weib.cp<-matrix(NA,nrow=nsims,ncol=15)

#create covariates
x<-runif(n, min=-2.5, max=12)
z<-log(runif(n, min=1, max=100))


#create a dependent variable, begin the simmulations
for(i in 1:nsims){

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
phi<-1/(1+exp(-(tru.est[i,3]+tru.est[i,4]*z+tru.est[i,5]*x)))
print(mean(phi))
yzero<-matrix(1,n,1)
error<--1*rlogis(n)
flag<-error<qlogis(phi)
yzero[flag]<-error[flag]
flag<-yzero==1
di[flag]<-ifelse(di[flag]==0,yzero[flag],di[flag])
tru.est[i,8]<-table(di)[1]

data<-cbind(ycen,di,x,z)


#####################################################################################
###################################COX Model#########################################
#####################################################################################

#store estimate and se
cox.est[i,1]<-summary(coxph(Surv(ycen, di)~x,coxph.control(iter.max = 10000)))$coef[1]
cox.est[i,2]<-summary(coxph(Surv(ycen, di)~x,coxph.control(iter.max = 10000)))$coef[3]

#store rmse
cox.rmse[i,1]<-sqrt((tru.est[i,2]-cox.est[i,1])^2)

#calculate upper and lower 95% CI's
b1.lower<-cox.est[i,1]-(1.959964*cox.est[i,2])
b1.upper<-cox.est[i,1]+(1.959964*cox.est[i,2])

#store coverage parameters
cox.cp[i,1]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)




#############################################################################
########################Simple Exponential Model#############################
#############################################################################


Exponential<- function(est,Y,C,X,data) {					      
	n=nrow(data)							      					  
	llik <- matrix(0, nrow=n, ncol = 1)
	beta<-est[1:length(est)]
	XB<-X%*%beta
	llik<-C*(XB-exp(XB)*Y)+(1-C)*(-exp(XB)*Y)
	llik<--1*sum(llik)
	return(llik)
	
	}



#set starting parameters
est<-rbind(.01,.01)

#set data, Y and X
data<-data
Y<-ycen
C<-di
X<-cbind(1,x)


#optimize
output.Exponential<-try(optim(f=Exponential,  p=est, X=X,Y=Y,C=C, method="BFGS", control=list(maxit=10000),  data=data, hessian=TRUE), TRUE)

if(class(output.Exponential)=="list"){
	ifelse(is.positive.definite(output.Exponential$hessian)==TRUE,vcv<-solve(output.Exponential$hessian),vcv<-matrix(data=NA,nrow=2,ncol=2))

#store betas and ses
exp.est[i,1]<-output.Exponential$par[1]
exp.est[i,2]<-sqrt(vcv[1,1])
exp.est[i,3]<-output.Exponential$par[2]
exp.est[i,4]<-sqrt(vcv[2,2])

#store rmse
exp.rmse[i,1]<-sqrt((tru.est[i,1]-exp.est[i,1])^2)
exp.rmse[i,2]<-sqrt((tru.est[i,2]-exp.est[i,3])^2)

#calculate upper and lower 95% CI's
b0.lower<-exp.est[i,1]-(1.959964*exp.est[i,2])
b0.upper<-exp.est[i,1]+(1.959964*exp.est[i,2])
b1.lower<-exp.est[i,3]-(1.959964*exp.est[i,4])
b1.upper<-exp.est[i,3]+(1.959964*exp.est[i,4])


#store coverage parameters
exp.cp[i,1]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
exp.cp[i,2]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)

}

################################################################################
#########################Simple Weibull Model ##################################
################################################################################

#Note this estiamtes the model via hazard rates, a la Stata

test<-survreg(Surv(ycen, di)~x, dist="weibull")
summary(test)


Weibull<- function(est,Y,C,X,data) {					      
	n=nrow(data)							      					  
	llik <- matrix(0, nrow=n, ncol = 1)
	beta<-est[1:length(est)-1]
	p<-est[length(est)]
	p<-exp(p)
	XB<-X%*%beta
	llik<-C*(log(exp(XB+1/p)*p*((exp(XB+1/p)*Y)^(p-1))*exp(-(exp(XB+1/p)*Y)^p)))+(1-C)*log(exp(-(exp(XB+1/p)*Y)^p))
	llik<--1*sum(llik)
	return(llik)
	
	}



#set starting parameters
est<-rbind(exp.est[i,1],exp.est[i,3],.01)

#set data, Y and X
data<-data
Y<-ycen
C<-di
X<-cbind(1,x)


#optimize
output.Weibull<-try(optim(f=Weibull,  p=est, X=X,Y=Y,C=C, method="BFGS", control=list(maxit=10000),  data=data, hessian=TRUE), TRUE)

if(class(output.Weibull)=="list"){
	ifelse(is.positive.definite(output.Weibull$hessian)==TRUE,vcv<-solve(output.Weibull$hessian),vcv<-matrix(data=NA,nrow=3,ncol=3))

#store betas and ses
weib.est[i,1]<-output.Weibull$par[1]+1/exp(output.Weibull$par[3])
coeff<-c(weib.est[i,1],output.Weibull$par[3])
varcov<-matrix(NA,2,2)
varcov[1,1]<-vcv[1,1]
varcov[1,2]<-vcv[1,3]
varcov[2,1]<-vcv[3,1]
varcov[2,2]<-vcv[3,3]
weib.est[i,2]<-deltamethod(~(x1+1/exp(x2)), coeff, varcov, ses=TRUE)
weib.est[i,3]<-output.Weibull$par[2]
weib.est[i,4]<-sqrt(vcv[2,2])
weib.est[i,5]<-exp(output.Weibull$par[3])
coeff<-c(weib.est[i,5])
varcov<-matrix(NA,1,1)
varcov[1,1]<-vcv[3,3]
weib.est[i,6]<-deltamethod(~(exp(x1)), coeff, varcov, ses=TRUE)
 

#store rmse
weib.rmse[i,1]<-sqrt((tru.est[i,1]-weib.est[i,1])^2)
weib.rmse[i,2]<-sqrt((tru.est[i,2]-weib.est[i,3])^2)
weib.rmse[i,3]<-sqrt((tru.est[i,6]-weib.est[i,5])^2)

#calculate upper and lower 95% CI's
b0.lower<-weib.est[i,1]-(1.959964*weib.est[i,2])
b0.upper<-weib.est[i,1]+(1.959964*weib.est[i,2])
b1.lower<-weib.est[i,3]-(1.959964*weib.est[i,4])
b1.upper<-weib.est[i,3]+(1.959964*weib.est[i,4])
p.lower<-weib.est[i,5]-(1.959964*weib.est[i,6])
p.upper<-weib.est[i,5]+(1.959964*weib.est[i,6])

#store coverage parameters
weib.cp[i,1]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
weib.cp[i,2]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
weib.cp[i,3]<-ifelse(tru.est[i,6]>p.lower & tru.est[i,6]<p.upper, 1,0)

}

###logit estimates###
dataset<-as.data.frame(data)
logitcoef1<-glm(di~ z+x, data = dataset, family = "binomial")$coef[1]
logitcoef2<-glm(di~ z+x, data = dataset, family = "binomial")$coef[2]
logitcoef3<-glm(di~ z+x, data = dataset, family = "binomial")$coef[3]

###############################################################################
##########################Zombie Exponential Model#############################
###############################################################################


#This program estimates the Exponential loglikelihood function returning hazard rate form coefficients

ZExponential<- function(est,Y,C,X,Z,data) {					      
  n=nrow(data)							      					  
  llik <- matrix(0, nrow=n, ncol = 1)
  gamma<-est[1:ncol(Z)]
  beta<-est[(ncol(Z)+1):length(est)]
  XB<-X%*%beta
  ZG<-Z%*%gamma
  phi<-1/(1+exp(-ZG))
  llik<-C*(log((1-phi)+phi*exp(XB)*exp(-exp(XB)*Y)))+(1-C)*(log(phi)+-exp(XB)*Y)
  llik<--1*sum(llik)
  return(llik)
  
}



#set starting parameters
est<-rbind(.01,.01,.01,exp.est[i,1],exp.est[i,3])

#set data, Y and X
data<-data
Y<-ycen
C<-di
X<-cbind(1,x)
Z<-cbind(1,z,x)


#optimize
output.ZExponential<-try(optim(f=ZExponential,  p=est, X=X,Y=Y,C=C,Z=Z, method="BFGS", control=list(maxit=10000),  data=data, hessian=TRUE), TRUE)

if(class(output.ZExponential)=="list"){
  ifelse(is.positive.definite(output.ZExponential$hessian)==TRUE,vcv<-solve(output.ZExponential$hessian),vcv<-matrix(data=NA,nrow=5,ncol=5))
  
  #store betas and ses
  exp.est[i,5]<-output.ZExponential$par[1]
  exp.est[i,6]<-sqrt(vcv[1,1])
  exp.est[i,7]<-output.ZExponential$par[2]
  exp.est[i,8]<-sqrt(vcv[2,2])
  exp.est[i,9]<-output.ZExponential$par[3]
  exp.est[i,10]<-sqrt(vcv[3,3])
  exp.est[i,11]<-output.ZExponential$par[4]
  exp.est[i,12]<-sqrt(vcv[4,4])
  exp.est[i,13]<-output.ZExponential$par[5]
  exp.est[i,14]<-sqrt(vcv[5,5])
  
  #store rmse
  exp.rmse[i,3]<-sqrt((tru.est[i,3]-exp.est[i,5])^2)
  exp.rmse[i,4]<-sqrt((tru.est[i,4]-exp.est[i,7])^2)
  exp.rmse[i,5]<-sqrt((tru.est[i,5]-exp.est[i,9])^2)
  exp.rmse[i,6]<-sqrt((tru.est[i,1]-exp.est[i,11])^2)
  exp.rmse[i,7]<-sqrt((tru.est[i,2]-exp.est[i,13])^2)
  
  #calculate upper and lower 95% CI's
  g0.lower<-exp.est[i,5]-(1.959964*exp.est[i,6])
  g0.upper<-exp.est[i,5]+(1.959964*exp.est[i,6])
  g1.lower<-exp.est[i,7]-(1.959964*exp.est[i,8])
  g1.upper<-exp.est[i,7]+(1.959964*exp.est[i,8])
  g2.lower<-exp.est[i,9]-(1.959964*exp.est[i,10])
  g2.upper<-exp.est[i,9]+(1.959964*exp.est[i,10])
  b0.lower<-exp.est[i,11]-(1.959964*exp.est[i,12])
  b0.upper<-exp.est[i,11]+(1.959964*exp.est[i,12])
  b1.lower<-exp.est[i,13]-(1.959964*exp.est[i,14])
  b1.upper<-exp.est[i,13]+(1.959964*exp.est[i,14])
  
  
  #store coverage parameters
  exp.cp[i,3]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
  exp.cp[i,4]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
  exp.cp[i,5]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
  exp.cp[i,6]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  exp.cp[i,7]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  
}

#####################################################################################
##########################Zombie Weibull Model ######################################
#####################################################################################


#This program estimates the Exponential loglikelihood function returning hazard rate form coefficients

ZWeibull<- function(est,Y,C,X,Z,data) {					      
  n=nrow(data)							      					  
  llik <- matrix(0, nrow=n, ncol = 1)
  gamma<-est[1:ncol(Z)]
  beta<-est[(ncol(Z)+1):(length(est)-1)]
  p<-est[length(est)]
  p<-exp(p)
  XB<-X%*%beta
  ZG<-Z%*%gamma
  phi<-1/(1+exp(-(ZG+1/p)))
  llik<-C*(log((1-phi)+phi*exp(XB+1/p)*p*((exp(XB+1/p)*Y)^(p-1))*exp(-(exp(XB+1/p)*Y)^p)))+(1-C)*(log(phi)+-(exp(XB+1/p)*Y)^p)
  llik<--1*sum(llik)
  return(llik)
  
}



#set starting parameters
est<-rbind(.01,.01,.01,output.Weibull$par[1],output.Weibull$par[2],output.Weibull$par[3])

#set data, Y and X
data<-data
Y<-ycen
C<-di
X<-cbind(1,x)
Z<-cbind(1,z,x)


#optimize
output.ZWeibull<-try(optim(f=ZWeibull,  p=est, X=X,Y=Y,C=C,Z=Z, method="BFGS", control=list(maxit=10000),  data=data, hessian=TRUE), TRUE)

if(class(output.ZWeibull)=="list"){
  ifelse(is.positive.definite(output.ZWeibull$hessian)==TRUE,vcv<-solve(output.ZWeibull$hessian),vcv<-matrix(data=NA,nrow=6,ncol=6))
  
  #store betas and ses
  weib.est[i,7]<-output.ZWeibull$par[1]+1/exp(output.ZWeibull$par[6])
  coeff<-c(weib.est[i,7],output.ZWeibull$par[6])
  varcov<-matrix(NA,2,2)
  varcov[1,1]<-vcv[1,1]
  varcov[1,2]<-vcv[1,6]
  varcov[2,1]<-vcv[6,1]
  varcov[2,2]<-vcv[6,6]
  weib.est[i,8]<-deltamethod(~(x1+1/exp(x2)), coeff, varcov, ses=TRUE)
  weib.est[i,9]<-output.ZWeibull$par[2]
  weib.est[i,10]<-sqrt(vcv[2,2])
  weib.est[i,11]<-output.ZWeibull$par[3]
  weib.est[i,12]<-sqrt(vcv[3,3])
  weib.est[i,13]<-output.ZWeibull$par[4]+1/exp(output.ZWeibull$par[6])
  coeff<-c(weib.est[i,13],output.ZWeibull$par[6])
  varcov<-matrix(NA,2,2)
  varcov[1,1]<-vcv[4,4]
  varcov[1,2]<-vcv[4,6]
  varcov[2,1]<-vcv[6,4]
  varcov[2,2]<-vcv[6,6]
  weib.est[i,14]<-deltamethod(~(x1+1/exp(x2)), coeff, varcov, ses=TRUE)
  weib.est[i,15]<-output.ZWeibull$par[5]
  weib.est[i,16]<-sqrt(vcv[5,5])
  weib.est[i,17]<-exp(output.ZWeibull$par[6])
  coeff<-c(weib.est[i,17])
  varcov<-matrix(NA,1,1)
  varcov[1,1]<-vcv[6,6]
  weib.est[i,18]<-deltamethod(~(exp(x1)), coeff, varcov, ses=TRUE)
  
  #store rmse
  weib.rmse[i,4]<-sqrt((tru.est[i,3]-weib.est[i,7])^2)
  weib.rmse[i,5]<-sqrt((tru.est[i,4]-weib.est[i,9])^2)
  weib.rmse[i,6]<-sqrt((tru.est[i,5]-weib.est[i,11])^2)
  weib.rmse[i,7]<-sqrt((tru.est[i,1]-weib.est[i,13])^2)
  weib.rmse[i,8]<-sqrt((tru.est[i,2]-weib.est[i,15])^2)
  weib.rmse[i,9]<-sqrt((tru.est[i,6]-weib.est[i,17])^2)
  
  #calculate upper and lower 95% CI's
  g0.lower<-weib.est[i,7]-(1.959964*weib.est[i,8])
  g0.upper<-weib.est[i,7]+(1.959964*weib.est[i,8])
  g1.lower<-weib.est[i,9]-(1.959964*weib.est[i,10])
  g1.upper<-weib.est[i,9]+(1.959964*weib.est[i,10])
  g2.lower<-weib.est[i,11]-(1.959964*weib.est[i,12])
  g2.upper<-weib.est[i,11]+(1.959964*weib.est[i,12])
  b0.lower<-weib.est[i,13]-(1.959964*weib.est[i,14])
  b0.upper<-weib.est[i,13]+(1.959964*weib.est[i,14])
  b1.lower<-weib.est[i,15]-(1.959964*weib.est[i,16])
  b1.upper<-weib.est[i,15]+(1.959964*weib.est[i,16])
  p.lower<-weib.est[i,17]-(1.959964*weib.est[i,18])
  p.upper<-weib.est[i,17]+(1.959964*weib.est[i,18])
  
  #store coverage parameters
  weib.cp[i,4]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
  weib.cp[i,5]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
  weib.cp[i,6]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
  weib.cp[i,7]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  weib.cp[i,8]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  weib.cp[i,9]<-ifelse(tru.est[i,6]>p.lower & tru.est[i,6]<p.upper, 1,0)
}


###############################################################################
######################Bayesian Zombie Exponential Model########################
###############################################################################
#set data, Y and X
data<-data
Y<-ycen
C<-di
X<-cbind(1,x)
Z<-cbind(1,z,x)
BayesZExponential = mcmcOF2(Y, C, X, Z, N = 3000, burn = 1000, thin = 20,  w = c(1, 1, 1), m = 10, form = "Exponential")
output.BayesZExponential = list(par = c(summary(mcmc(BayesZExponential$beta))[[1]][,1], summary(mcmc(BayesZExponential$gamma))[[1]][,1]), 
								se = c(summary(mcmc(BayesZExponential$beta))[[1]][,2], summary(mcmc(BayesZExponential$gamma))[[1]][,2]),
								CI = rbind(summary(mcmc(BayesZExponential$beta))[[2]], summary(mcmc(BayesZExponential$gamma))[[2]]))
exp.est[i,15]<-output.BayesZExponential$par[1]
exp.est[i,16]<-output.BayesZExponential$se[1]
exp.est[i,17]<-output.BayesZExponential$par[2]
exp.est[i,18]<-output.BayesZExponential$se[2]
exp.est[i,19]<-output.BayesZExponential$par[3]
exp.est[i,20]<-output.BayesZExponential$se[3]
exp.est[i,21]<-output.BayesZExponential$par[4]
exp.est[i,22]<-output.BayesZExponential$se[4]
exp.est[i,23]<-output.BayesZExponential$par[5]
exp.est[i,24]<-output.BayesZExponential$se[5]

#store rmse
exp.rmse[i,8]<-sqrt((tru.est[i,3]-exp.est[i,15])^2)
exp.rmse[i,9]<-sqrt((tru.est[i,4]-exp.est[i,17])^2)
exp.rmse[i,10]<-sqrt((tru.est[i,5]-exp.est[i,19])^2)
exp.rmse[i,11]<-sqrt((tru.est[i,1]-exp.est[i,21])^2)
exp.rmse[i,12]<-sqrt((tru.est[i,2]-exp.est[i,23])^2)

#calculate upper and lower 95% CI's
# b0.lower<-output.BayesZExponential$CI[1,1]
# b0.upper<-output.BayesZExponential$CI[1,5]
# b1.lower<-output.BayesZExponential$CI[2,1]
# b1.upper<-output.BayesZExponential$CI[2,5]
# g0.lower<-output.BayesZExponential$CI[3,1]
# g0.upper<-output.BayesZExponential$CI[3,5]
# g1.lower<-output.BayesZExponential$CI[4,1]
# g1.upper<-output.BayesZExponential$CI[4,5]
# g2.lower<-output.BayesZExponential$CI[5,1]
# g2.upper<-output.BayesZExponential$CI[5,5]
b0.lower<-exp.est[i,15]-(1.959964*exp.est[i,16])
b0.upper<-exp.est[i,15]+(1.959964*exp.est[i,16])
b1.lower<-exp.est[i,17]-(1.959964*exp.est[i,18])
b1.upper<-exp.est[i,17]+(1.959964*exp.est[i,18])
g0.lower<-exp.est[i,19]-(1.959964*exp.est[i,20])
g0.upper<-exp.est[i,19]+(1.959964*exp.est[i,20])
g1.lower<-exp.est[i,21]-(1.959964*exp.est[i,22])
g1.upper<-exp.est[i,21]+(1.959964*exp.est[i,22])
g2.lower<-exp.est[i,23]-(1.959964*exp.est[i,24])
g2.upper<-exp.est[i,23]+(1.959964*exp.est[i,24])


#store coverage parameters
exp.cp[i,8]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
exp.cp[i,9]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
exp.cp[i,10]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
exp.cp[i,11]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
exp.cp[i,12]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)

###############################################################################
########################Bayesian Zombie Weibull Model##########################
###############################################################################
#set data, Y and X
data<-data
Y<-ycen
C<-di
X<-cbind(1,x)
Z<-cbind(1,z,x)
BayesZWeibull = mcmcOF2(Y, C, X, Z, N = 3000, burn = 1000, thin = 20,  w = c(1, 1, 1), m = 10, form = "Weibull")
output.BayesZWeibull = list(par = c(summary(mcmc(BayesZWeibull$beta))[[1]][,1], summary(mcmc(BayesZWeibull$gamma))[[1]][,1], 
									summary(mcmc(BayesZWeibull$lambda))[[1]][1]), 
								se = c(summary(mcmc(BayesZWeibull$beta))[[1]][,2], summary(mcmc(BayesZWeibull$gamma))[[1]][,2], 
									   summary(mcmc(BayesZWeibull$lambda))[[1]][2]),
								CI = rbind(summary(mcmc(BayesZWeibull$beta))[[2]], summary(mcmc(BayesZWeibull$gamma))[[2]], 
										summary(mcmc(BayesZWeibull$lambda))[[2]]))

weib.est[i,19]<-output.BayesZWeibull$par[1]
weib.est[i,20]<-output.BayesZWeibull$se[1]
weib.est[i,21]<-output.BayesZWeibull$par[2]
weib.est[i,22]<-output.BayesZWeibull$se[2]
weib.est[i,23]<-output.BayesZWeibull$par[3]
weib.est[i,24]<-output.BayesZWeibull$se[3]
weib.est[i,25]<-output.BayesZWeibull$par[4]
weib.est[i,26]<-output.BayesZWeibull$se[4]
weib.est[i,27]<-output.BayesZWeibull$par[5]
weib.est[i,28]<-output.BayesZWeibull$se[5]
weib.est[i,29]<-output.BayesZWeibull$par[6]
weib.est[i,30]<-output.BayesZWeibull$se[6]

#store rmse
weib.rmse[i,10]<-sqrt((tru.est[i,3]-weib.est[i,19])^2)
weib.rmse[i,11]<-sqrt((tru.est[i,4]-weib.est[i,21])^2)
weib.rmse[i,12]<-sqrt((tru.est[i,5]-weib.est[i,23])^2)
weib.rmse[i,13]<-sqrt((tru.est[i,1]-weib.est[i,25])^2)
weib.rmse[i,14]<-sqrt((tru.est[i,2]-weib.est[i,27])^2)
weib.rmse[i,15]<-sqrt((tru.est[i,6]-weib.est[i,29])^2)

#calculate upper and lower 95% CI's
# b0.lower<-output.BayesZWeibull$CI[1,1]
# b0.upper<-output.BayesZWeibull$CI[1,5]
# b1.lower<-output.BayesZWeibull$CI[2,1]
# b1.upper<-output.BayesZWeibull$CI[2,5]
# g0.lower<-output.BayesZWeibull$CI[3,1]
# g0.upper<-output.BayesZWeibull$CI[3,5]
# g1.lower<-output.BayesZWeibull$CI[4,1]
# g1.upper<-output.BayesZWeibull$CI[4,5]
# g2.lower<-output.BayesZWeibull$CI[5,1]
# g2.upper<-output.BayesZWeibull$CI[5,5]
# p.lower<-output.BayesZWeibull$CI[6,1]
# p.upper<-output.BayesZWeibull$CI[6,2]
g0.lower<-weib.est[i,19]-(1.959964*weib.est[i,20])
g0.upper<-weib.est[i,19]+(1.959964*weib.est[i,20])
g1.lower<-weib.est[i,21]-(1.959964*weib.est[i,22])
g1.upper<-weib.est[i,21]+(1.959964*weib.est[i,22])
g2.lower<-weib.est[i,23]-(1.959964*weib.est[i,24])
g2.upper<-weib.est[i,23]+(1.959964*weib.est[i,24])
b0.lower<-weib.est[i,25]-(1.959964*weib.est[i,26])
b0.upper<-weib.est[i,25]+(1.959964*weib.est[i,26])
b1.lower<-weib.est[i,27]-(1.959964*weib.est[i,28])
b1.upper<-weib.est[i,27]+(1.959964*weib.est[i,28])
p.lower<-weib.est[i,29]-(1.959964*weib.est[i,30])
p.upper<-weib.est[i,29]+(1.959964*weib.est[i,30])


#store coverage parameters
weib.cp[i,10]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
weib.cp[i,11]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
weib.cp[i,12]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
weib.cp[i,13]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
weib.cp[i,14]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
weib.cp[i,15]<-ifelse(tru.est[i,6]>p.lower & tru.est[i,6]<p.upper, 1,0)

}
#combine matrices and label variables
main.data<-cbind(tru.est, cox.est, exp.est, weib.est, cox.rmse, exp.rmse, weib.rmse, cox.cp, exp.cp, weib.cp)
colnames(main.data)<-c("true.x0","true.x1","true.z0","true.z1","true.z2","true.p","cen.lat","cen.obs",
	"cox.x1","cox.x1.se",
	"exp.x0","exp.x0.se","exp.x1","exp.x1.se",
	"zexp.z0","zexp.z0.se","zexp.z1","zexp.z1.se","zexp.z2","zexp.z2.se","zexp.x0","zexp.x0.se","zexp.x1","zexp.x1.se",
	"bzexp.x0","zexp.x0.se","bzexp.x1","bzexp.x1.se","bzexp.z0","bzexp.z0.se","bzexp.z1","bzexp.z1.se","bzexp.z2","bzexp.z2.se",
	"wei.x0","wei.x0.se","wei.x1","wei.x1.se","wei.p","wei.p.se",
	"zwei.z0","zwei.z0.se","zwei.z1","zwei.z1.se","zwei.z2","zwei.z2.se","zwei.x0","zwei.x0.se","zwei.x1","zwei.x1.se","zwei.p","zwei.p.se",
	"bzwei.x0","bzwei.x0.se","bzwei.x1","bzwei.x1.se","bzwei.z0","bzwei.z0.se","bzwei.z1","bzwei.z1.se","bzwei.z2","bzwei.z2.se","bzwei.p","bzwei.p.se",
	"cox.x1.rmse",
	"exp.x0.rmse","exp.x1.rmse","zexp.z0.rmse","zexp.z1.rmse","zexp.z2.rmse","zexp.x0.rmse","zexp.x1.rmse","bzexp.x0.rmse","bzexp.x1.rmse","bzexp.z0.rmse","bzexp.z1.rmse","bzexp.z2.rmse",
	"wei.x0.rmse","wei.x1.rmse","wei.p.rmse","zwei.z0.rmse","zwei.z1.rmse","zwei.z2.rmse",
	"zwei.x0.rmse","zwei.x1.rmse","zwei.p.rmse", "bzwei.x0.rmse","bzwei.x1.rmse","bzwei.z0.rmse","bzwei.z1.rmse","bzwei.z2.rmse","bzwei.p.rmse",
	"cox.x1.cp","exp.x0.cp","exp.x1.cp","zexp.z0.cp","zexp.z1.cp","zexp.z2.cp","zexp.x0.cp","zexp.x1.cp","bzexp.x0.cp","bzexp.x1.cp","bzexp.z0.cp","bzexp.z1.cp","bzexp.z2.cp",
	"wei.x0.cp","wei.x1.cp","wei.p.cp",
	"zwei.z0.cp","zwei.z1.cp","zwei.z2.cp","zwei.x0.cp","zwei.x1.cp","zwei.p.cp", "bzwei.x0.cp","bzwei.x1.cp","bzwei.z0.cp","bzwei.z1.cp","bzwei.z2.cp","bzwei.p.cp")

#save dataset
main.data<-as.data.frame(main.data)
write.dta(main.data2,"main.data2.dta", )

#the end