rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1)
integrate(function(y) 2*exp(-2*(y-1)), lower=0, upper=4)

#Problem 2)
#Find P(X=1) for f(x)=  2^x/x! e^(-2)
dpois(1,2)


#find P(-2 < X < 4)
rm(list=ls())  #clears environment
x<-0:4
y<-dpois(0:4,2)
data.frame("fx"=y,row.names=x)
plot(0:4, dpois(0:4,2), type='h', xlab="x", ylab="f(x)" )

#Problem 3)
#no r code needed

#Problem 4)
rm(list=ls())  #clears environment
#P(Y ??? 2) 
sum(dbinom(0:2,3,0.25))

#E(Y)
y.range<-(0:3)
Ey<-sum(y.range*dbinom(y.range,size=3,p=0.25))
Ey

#Var(Y)
y.range<-(0:3)
Vary<-sum((y.range-Ey)^2*dbinom(y.range,size=3,p=0.25))
Vary

#Problem 5
rm(list=ls())  #clears environment
#For X following a Chi-square distribution with degree of freedom m = 3, 
#compute the following: 
#P(1 < X < 4)
#note, if:
#X<y pchisq(y, m)
#X>y 1-pchisq(y, m)
#y<x<z pchisq(z, m) - pchisq(y, m)
pchisq(4,3)-pchisq(1,3)

#Monte Carlo rchisq(n, df, ncp = 0) where n =100,000 and df = 3
x<-rchisq(100000, 3, ncp = 0)
mean((1<x) & (x<4))

#Problem 6)
#no r code needed

#Problem 7)
rm(list=ls())  #clears environment
#a)
#What is the probability that a randomly chosen patient have the Zyxin gene 
#expression values between 1 and 1.6 when mean = 1.6 and sd = 0.4
#confirm my equation:
integrate(function(x) exp(-(x-1.6)^2/0.32)/0.4/sqrt(2*pi), lower=1, upper=1.6)
#samething as above but having R solve it
integrate(function(x) dnorm(x,mean=1.6,sd=0.4), lower=1, upper=1.6)

#b) Use a Monte Carlo simulation of sample size n=500,000 to estimate the 
#probability in part (a). Give your R code, and show the value of your estimate. 
#use rnorm(n, mean = 0, sd = 1)
x<-rnorm(500000, mean = 1.6, sd = 0.4)
mean((1<x) & (x<1.6))

#c)What is the probability that exactly 2 out of 5 patients have the Zyxin gene 
#expression values between 1 and 1.6? 
y<-dbinom(2,5,0.4331928)
y

#8) 
rm(list=ls())  #clears environment
#(a) Hand in a R script that calculates the mean and variance of two random variables 
#X~F(m=2,n=5) and Y~F(m=10,n=5) from their density functions. 

#for X~F(m=2,n=5) where o = number of observations = 1,000,000
#use rf(o, m, n)
x<-rf(1000000, 2, 5)

mean(x)
var(x)

#for Y~F(m=10,n=5) where o = number of observations = 1,000,000
#use rf(o, m, n)
y<-rf(1000000, 10, 5)

mean(y)
var(y)