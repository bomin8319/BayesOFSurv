
llFun <- function(est,Y,C,X,Z,data) {
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
  llik<- -1*sum(llik)
  return(llik)
}

calculateDIC = function(Y,C,X,Z,data, theta, llFun) {
  #Calculate L
  theta_post = cbind(theta$gammas, theta$betas, theta$lambda)
  theta_hat = apply(theta_post, 2, mean)
  L = llFun(theta_hat,Y,C,X,Z,data)
  
  #Calculate P
  S = nrow(theta_post) #S = number of iterations
  #Add up the log likelihoods of each iteration
  llSum = 0
  for (s in 1:S) {
    theta_s = theta_post[s,]
    llSum = llSum + llFun(theta_s,Y,C,X,Z,data)
  }
  
  P = 2 * (L - (1 / S * llSum))
  
  #Calculate DIC
  DIC = -2 * (L - P)
  
  #Return the results
  list(DIC=DIC, P=P, L=L)
}

#---------------------------------------------------------------------;


X <- cbind(1, rbs$lnlevel, rbs$necon, rbs$presi, rbs$tag, rbs$rel, rbs$ethn, rbs$prevdem, (rbs$openc)/10, rbs$calinv, rbs$calileve)
X <- as.matrix(X)
Z <- cbind(rep(1,nrow(rbs)))
Z <- as.matrix(Z)
Y <- as.matrix(rbs$'_t')
C <- as.matrix(rbs$'_d')


