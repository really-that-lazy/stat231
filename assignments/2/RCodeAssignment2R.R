###################################################################################
# run this code only once 
library(moments)
library(MASS)     # truehist is in the library MASS
#
###################################################################################
cat("\n\nProblem 1\n")
###################################################################################
# Problem 1: R code for Binomial relative likelihood function
# The following R code will plot the Binomial relative likelihood function and the line y = 0.15
# which can be used to determine a 15% likelihood interval.
##id<-20456484
id<-20466075
set.seed(id)
#generate a random value of theta from Beta distribution last 2 digits of ID avoiding zero values
theta<-rbeta(1, max(1,id-10*trunc(id/10)), max(1,trunc(id/10)-10*trunc(id/100)))  
n<-40
y<-rbinom(1,n,theta)    # observation from Binomial(n,theta) distribution
thetahat<-y/n    #maximum likelihood estimate of theta
# display values
cat("n = ", n, ", theta = ", theta, ", y = ", y, ", thetahat = ",thetahat)
s<-(thetahat*(1-thetahat)/n)^0.5
# interval of values for plotting relative likelihood function
th<-seq(max(0,thetahat-4*s), min(1,thetahat+4*s),0.001)
# create function to calculate Binomial relative likelihood function
BinRLF <- function(x) {exp(y*log(x/thetahat)+(n-y)*log((1-x)/(1-thetahat)))}

png("problem1-1.png", width=680, height=480, res=120)
plot(th,BinRLF(th),xlab="theta",type="l",lwd=2)    # plot relative likelihood function
# draw a horizontal line at 0.15
abline(a=0.15,b=0,col="red",lwd=2)
title(main="Binomial Relative Likelihood Function")
dev.off()
###################################################################################

###################################################################################
# R code for finding the endpoints of a 15% likelihood interval for Binomial example
cat("Root Finder")
uniroot(function(x) BinRLF(x)-0.15,lower=0.1,upper=0.20)
uniroot(function(x) BinRLF(x)-0.15,lower=0.35,upper=0.62)
###################################################################################
cat("\n\nProblem 2\n")
###################################################################################
# Problem 2: R code for Exponential relative likelihood function
# The following R code will plot the Exponential relative likelihood function and the line y = 0.15
# which can be used to determine a 15% likelihood interval.
set.seed(id)
theta<-rexp(1, 1/(max(1,id-10*trunc(id/10))))    #generate a random value of theta
n<-40
# generate random sample of n observation from Exponential(theta) distribution 
y<-sort(rexp(n,1/theta))    
y[1:3]                 # display first 3 numbers in the data set
thetahat<-mean(y)    #maximum likelihood estimate of theta
cat("n = ",n,", theta = ",theta,", thetahat = ",thetahat)   # display values
s<-thetahat/(n^0.5)
# interval of values for plotting relative likelihood function
th<-seq(max(0,thetahat-4*s), thetahat+4*s,0.001)
# create function to calculate Exponential relative likelihood function 
ExpRLF<-function(x) {(thetahat/x)^n*exp(n*(1-thetahat/x))}

png("problem2-1.png", width=680, height=480, res=120)
plot(th,ExpRLF(th),xlab="theta",type="l",lwd=2)    # plot relative likelihood function
# draw a horizontal line at 0.15
abline(a=0.15,b=0,col="red",lwd=2)
title(main="Exponential Relative Likelihood Function")
dev.off()
###################################################################################

###################################################################################
# R code for finding the endpoints of a 15% likelihood interval for Exponential example
uniroot(function(x) ExpRLF(x)-0.15,lower=2.8,upper=4)
uniroot(function(x) ExpRLF(x)-0.15,lower=5,upper=8)
###################################################################################
cat("\n\nProblem 3\n")
#############################################################################
# Problem 3: Sampling Distribution
set.seed(id)
# pop is a vector of 500 variate values for a finite population
pop<-rpois(500,max(1,id-10*trunc(id/10)))
mu<-mean(pop)             # population mean
stdeviation<-sd(pop)
cat('mu')
mu
cat('stdev')
stdeviation
(499*var(pop)/500)^0.5    # population standard deviation
png("problem3-1.png", width=680, height=480, res=120)
truehist(pop,h=1,main="Population Histogram",xlab="Variate Value",)
dev.off()
k<-10000           # number of simulations
n<-15              # sample size
sim<-rep(0,k)      # vector to store sample means
# Calculate k sample means for samples of size n drawn from population pop
for (i in 1:k)
sim[i]=mean(sample(pop,n,replace=F)) 

png("problem3-2.png", width=680, height=480, res=120)
truehist(sim,h=0.25,xlab="Sample Mean",main="Sampling Distribution of Sample Mean")
# percentage of times sample mean is within 0.5 of true mean mu
mean(abs(sim-mu)<0.5)
dev.off()
#############################################################################

