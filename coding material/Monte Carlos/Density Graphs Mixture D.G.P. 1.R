#graph B2

##Ben Bagozzi
##
##June 8th 2011
##
##Combine MC results
##
##
##
##

#clear memory  
rm(list=ls())                                           
library(car)
library(Hmisc)                                                     
library(mvtnorm)
library(foreign)
library(graphics)
library(MASS)
library(lattice)
library(tseries)
library(Matrix)
library(Design)
library(msm)
library(corpcor)
library(Zelig)


Results<-read.dta("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Mixture DGP/main.data1000.dta")
Results<-as.matrix(Results, )

resize.win <- function(Width=6, Height=6)
{
        # works for windows
    dev.off(); # dev.new(width=6, height=6)
    windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(2,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(-0.3,1.6)
local.ylim<-c(0,6.6)
plot(density(Results[,11],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.7, y =6.4, "BFGS.Exp",col="red")
text (x=.7, y =6, "BFGS.Wei",col="forestgreen")
text (x=1.33, y =6.25, "Exp",col="blue")
text (x=1.5, y =6.25, "Wei",col="purple")
text (x=0.7, y =5.5, "Bayes.Exp",col="pink")
text (x=0, y =6.25, "Bayes.Wei",col="orange")
par(new=FALSE)
 

#B1
local.xlim<-c(3.40,3.6)
local.ylim<-c(0,46.5)
plot(density(Results[,13],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,23],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,49],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="forestgreen")

text (x=3.55, y =43, "BFGS.Exp",col="red")
text (x=3.425, y =40, "BFGS.Wei",col="purple")
text (x=3.425, y =43, "Exp",col="blue")
text (x=3.55, y =34, "Wei",col="forestgreen")
text (x=3.55, y =40, "Bayes.Exp",col="pink")
text (x=3.55, y =37, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,20)
local.xlim<-c(0.75,1.1)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=1.06, y =15.5, "BFGS.Wei",col="purple")
text (x=1.06, y =14.2, "Bayes.Wei",col="orange")
text (x=.80, y =18, "Wei",col="forestgreen")
par(new=FALSE)

#B1
local.ylim<-c(0,.5)
local.xlim<-c(-8,2.5)
plot(density(Results[,15],na.rm=TRUE), main = "", ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=-2,lty=3)
par(new=TRUE)
plot(density(Results[,29],na.rm=TRUE), main = "", xlab = "Gamma 0",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,41],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,57],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=-4, y =.38, "BFGS.Exp",col="red")
text (x=1, y =.38, "Bayes.Exp",col="pink")
text (x=1, y =.35, "Bayes.Wei",col="orange")
text (x=-4, y =0.35, "BFGS.Wei",col="purple")
par(new=FALSE)


#B1
local.ylim<-c(0,1.25)
local.xlim<-c(0.5,6)
plot(density(Results[,17],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=2,lty=3)
par(new=TRUE)
plot(density(Results[,31],na.rm=TRUE), main = "", xlab = "Gamma 1",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,43],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,59],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=2.75, y =.9, "BFGS.Exp",col="red")
text (x=0.9, y =.82, "Bayes.Exp",col="pink")
text (x=2.75, y =.80, "BFGS.Wei",col="purple")
text (x=0.9, y =.92, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,0.8)
local.xlim<-c(1,8.5)
plot(density(Results[,19],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=3,lty=3)
par(new=TRUE)
plot(density(Results[,33],na.rm=TRUE), main = "", xlab = "Gamma 2",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,45],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,61],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=4.5, y =.6, "BFGS.Exp",col="red")
text (x=4.5, y =.65, "Bayes.Exp",col="pink")
text (x=4.5, y =.55, "BFGS.Wei",col="purple")
text (x=4.5, y =.70, "Bayes.Wei",col="orange")
par(new=FALSE)





Results<-read.dta("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Regular DGP/main.data.dta")
Results<-as.matrix(Results, )

resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(1,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(-0.2,1.6)
local.ylim<-c(0,6.5)
plot(density(Results[,11],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.7, y =6.5, "BFGS.Exp",col="red")
text (x=.7, y =6.3, "BFGS.Wei",col="forestgreen")
text (x=1.33, y =6.25, "Exp",col="blue")
text (x=1.33, y =6, "Wei",col="purple")
text (x=0.7, y =6, "Bayes.Exp",col="pink")
text (x=0, y =6.25, "Bayes.Wei",col="orange")
par(new=FALSE)


#B1
local.xlim<-c(3.40,3.6)
local.ylim<-c(0,46.5)
plot(density(Results[,13],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,23],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,49],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="forestgreen")

text (x=3.55, y =42, "BFGS.Exp",col="red")
text (x=3.46, y =40, "BFGS.Wei",col="purple")
text (x=3.46, y =42, "Exp",col="blue")
text (x=3.55, y =36, "Wei",col="forestgreen")
text (x=3.55, y =40, "Bayes.Exp",col="pink")
text (x=3.55, y =38, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,19)
local.xlim<-c(0.90,1.12)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=1.06, y =15.2, "BFGS.Wei",col="purple")
text (x=1.06, y =14.2, "Bayes.Wei",col="orange")
text (x=1.06, y =13.2, "Wei",col="forestgreen")
par(new=FALSE)





#1500

Results<-read.dta("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Mixture DGP/main.data1500.dta")
Results<-as.matrix(Results, )

resize.win <- function(Width=6, Height=6)
{
        # works for windows
    dev.off(); # dev.new(width=6, height=6)
    windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(2,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(-0.4,1.6)
local.ylim<-c(0,8)
plot(density(Results[,11],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.7, y =6.4, "BFGS.Exp",col="red")
text (x=.7, y =6, "BFGS.Wei",col="forestgreen")
text (x=1.5, y =6.75, "Exp",col="blue")
text (x=1.5, y =6.25, "Wei",col="purple")
text (x=0.7, y =5.5, "Bayes.Exp",col="pink")
text (x=-0.25, y =6.25, "Bayes.Wei",col="orange")
par(new=FALSE)
 

#B1
local.xlim<-c(3.40,3.6)
local.ylim<-c(0,55)
plot(density(Results[,13],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,23],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,49],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="forestgreen")

text (x=3.55, y =50, "BFGS.Exp",col="red")
text (x=3.425, y =38, "BFGS.Wei",col="purple")
text (x=3.425, y =42, "Exp",col="blue")
text (x=3.55, y =38, "Wei",col="forestgreen")
text (x=3.55, y =46, "Bayes.Exp",col="pink")
text (x=3.55, y =42, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,25)
local.xlim<-c(0.75,1.2)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=1.08, y =16.2, "BFGS.Wei",col="purple")
text (x=1.08, y =14.2, "Bayes.Wei",col="orange")
text (x=.78, y =18, "Wei",col="forestgreen")
par(new=FALSE)

#B1
local.ylim<-c(0,.6)
local.xlim<-c(-6,3)
plot(density(Results[,15],na.rm=TRUE), main = "", ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=-2,lty=3)
par(new=TRUE)
plot(density(Results[,29],na.rm=TRUE), main = "", xlab = "Gamma 0",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,41],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,57],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=-4, y =.38, "BFGS.Exp",col="red")
text (x=1, y =.38, "Bayes.Exp",col="pink")
text (x=1, y =.35, "Bayes.Wei",col="orange")
text (x=-4, y =0.35, "BFGS.Wei",col="purple")
par(new=FALSE)


#B1
local.ylim<-c(0,1.3)
local.xlim<-c(.5,4.5)
plot(density(Results[,17],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=2,lty=3)
par(new=TRUE)
plot(density(Results[,31],na.rm=TRUE), main = "", xlab = "Gamma 1",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,43],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,59],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=2.75, y =.86, "BFGS.Exp",col="red")
text (x=1, y =.86, "Bayes.Exp",col="pink")
text (x=2.75, y =.80, "BFGS.Wei",col="purple")
text (x=1, y =.92, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,0.95)
local.xlim<-c(1,7)
plot(density(Results[,19],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=3,lty=3)
par(new=TRUE)
plot(density(Results[,33],na.rm=TRUE), main = "", xlab = "Gamma 2",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,45],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,61],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=4.4, y =.6, "BFGS.Exp",col="red")
text (x=4.4, y =.65, "Bayes.Exp",col="pink")
text (x=4.4, y =.55, "BFGS.Wei",col="purple")
text (x=4.4, y =.70, "Bayes.Wei",col="orange")
par(new=FALSE)





Results<-read.dta("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Regular DGP/main.data1500_2.dta")
Results<-as.matrix(Results, )

resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(1,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(-0.2,1.6)
local.ylim<-c(0,7.6)
plot(density(Results[,11],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.7, y =6.5, "BFGS.Exp",col="red")
text (x=.7, y =6.3, "BFGS.Wei",col="forestgreen")
text (x=1.33, y =6.25, "Exp",col="blue")
text (x=1.33, y =6, "Wei",col="purple")
text (x=0.7, y =6, "Bayes.Exp",col="pink")
text (x=0.2, y =6.25, "Bayes.Wei",col="orange")
par(new=FALSE)


#B1
local.xlim<-c(3.40,3.6)
local.ylim<-c(0,55)
plot(density(Results[,13],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,23],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,49],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="forestgreen")

text (x=3.55, y =42, "BFGS.Exp",col="red")
text (x=3.46, y =40, "BFGS.Wei",col="purple")
text (x=3.46, y =42, "Exp",col="blue")
text (x=3.55, y =36, "Wei",col="forestgreen")
text (x=3.55, y =40, "Bayes.Exp",col="pink")
text (x=3.55, y =38, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,19)
local.xlim<-c(0.90,1.12)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=1.06, y =15.2, "BFGS.Wei",col="purple")
text (x=1.06, y =14.2, "Bayes.Wei",col="orange")
text (x=1.06, y =13.2, "Wei",col="forestgreen")
par(new=FALSE)





#2000

Results<-read.dta("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Mixture DGP/main.data2000.dta")
Results<-as.matrix(Results, )

resize.win <- function(Width=6, Height=6)
{
        # works for windows
    dev.off(); # dev.new(width=6, height=6)
    windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(2,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(-0.45,1.6)
local.ylim<-c(0,9.5)
plot(density(Results[,11],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.7, y =6.4, "BFGS.Exp",col="red")
text (x=.7, y =6, "BFGS.Wei",col="forestgreen")
text (x=1.5, y =6.75, "Exp",col="blue")
text (x=1.5, y =6.25, "Wei",col="purple")
text (x=0.7, y =5.5, "Bayes.Exp",col="pink")
text (x=-0.25, y =6.25, "Bayes.Wei",col="orange")
par(new=FALSE)
 

#B1
local.xlim<-c(3.40,3.6)
local.ylim<-c(0,65)
plot(density(Results[,13],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,23],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,49],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="forestgreen")

text (x=3.55, y =50, "BFGS.Exp",col="red")
text (x=3.425, y =38, "BFGS.Wei",col="purple")
text (x=3.425, y =42, "Exp",col="blue")
text (x=3.55, y =38, "Wei",col="forestgreen")
text (x=3.55, y =46, "Bayes.Exp",col="pink")
text (x=3.55, y =42, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,25)
local.xlim<-c(0.75,1.2)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=1.09, y =16.2, "BFGS.Wei",col="purple")
text (x=1.09, y =14.2, "Bayes.Wei",col="orange")
text (x=.78, y =18, "Wei",col="forestgreen")
par(new=FALSE)

#B1
local.ylim<-c(0,.6)
local.xlim<-c(-6,3)
plot(density(Results[,15],na.rm=TRUE), main = "", ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=-2,lty=3)
par(new=TRUE)
plot(density(Results[,29],na.rm=TRUE), main = "", xlab = "Gamma 0",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,41],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,57],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=-4, y =.38, "BFGS.Exp",col="red")
text (x=1, y =.38, "Bayes.Exp",col="pink")
text (x=1, y =.35, "Bayes.Wei",col="orange")
text (x=-4, y =0.35, "BFGS.Wei",col="purple")
par(new=FALSE)


#B1
local.ylim<-c(0,1.4)
local.xlim<-c(.5,4.5)
plot(density(Results[,17],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=2,lty=3)
par(new=TRUE)
plot(density(Results[,31],na.rm=TRUE), main = "", xlab = "Gamma 1",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,43],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,59],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=2.75, y =.9, "BFGS.Exp",col="red")
text (x=1, y =.82, "Bayes.Exp",col="pink")
text (x=2.75, y =.80, "BFGS.Wei",col="purple")
text (x=1, y =.92, "Bayes.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,1)
local.xlim<-c(1,7)
plot(density(Results[,19],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=3,lty=3)
par(new=TRUE)
plot(density(Results[,33],na.rm=TRUE), main = "", xlab = "Gamma 2",xlim=local.xlim, ylim=local.ylim,col="pink")
par(new=TRUE)
plot(density(Results[,45],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple")
par(new=TRUE)
plot(density(Results[,61],na.rm=TRUE), main = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange")
text (x=4.4, y =.6, "BFGS.Exp",col="red")
text (x=4.4, y =.65, "Bayes.Exp",col="pink")
text (x=4.4, y =.55, "BFGS.Wei",col="purple")
text (x=4.4, y =.70, "Bayes.Wei",col="orange")
par(new=FALSE)





Results<-read.dta("/Users/bomin8319/Desktop/BayesOFsurv/coding material/Monte Carlos/Regular DGP/main.data2000_2.dta")
Results<-as.matrix(Results, )

resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(1,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(-0.2,1.6)
local.ylim<-c(0,9)
plot(density(Results[,11],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.7, y =6.5, "Z.Exp",col="red")
text (x=.7, y =6.3, "Z.Wei",col="forestgreen")
text (x=1.33, y =6.25, "Exp",col="blue")
text (x=1.33, y =6, "Wei",col="purple")
text (x=0.7, y =6, "BZ.Exp",col="pink")
text (x=0.2, y =6.25, "BZ.Wei",col="orange")
par(new=FALSE)


#B1
local.xlim<-c(3.40,3.6)
local.ylim<-c(0,60)
plot(density(Results[,13],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,23],na.rm=TRUE), main = "", ylab = "",  xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="purple", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,49],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="forestgreen")

text (x=3.55, y =42, "Z.Exp",col="red")
text (x=3.46, y =40, "Z.Wei",col="purple")
text (x=3.46, y =42, "Exp",col="blue")
text (x=3.55, y =36, "Wei",col="forestgreen")
text (x=3.55, y =40, "BZ.Exp",col="pink")
text (x=3.55, y =38, "BZ.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,22)
local.xlim<-c(0.90,1.12)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=1.06, y =15.2, "Z.Wei",col="purple")
text (x=1.06, y =14.2, "BZ.Wei",col="orange")
text (x=1.06, y =13.2, "Wei",col="forestgreen")
par(new=FALSE)
