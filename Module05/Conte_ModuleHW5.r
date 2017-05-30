rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#numerical value for MLE
obs<-c(1.636, 0.374, 0.534, 3.015, 0.932, 0.179)
nl.x=function(x) -sum(log(dexp(obs,x)))#negative-likelihood function
optim(1,nl.x) #minimize nl.x with starting parameter value=1

#analytic formula
1/c(mean(obs))

#Problem 2
rm(list=ls())  #clears environment
x<-100.8
n<-53
s<-12.4

qt(c(0.05),df=n-1)

ci <- x+qt(c(0.05),df=n-1)*s/sqrt(n) #90% t-interval
ci

#Problem 3
rm(list=ls())  #clears environment
data(golub, package = "multtest") #load "golub" data in package "multtest" 
gol.fac <- factor(golub.cl,levels=0:1,labels=c("ALL","AML")) #get a factor variable, label values 0/1 to ALL/AML
#define ALL and AML variables
Zyxin.ALL<-golub[2124,gol.fac=="ALL"]
Zyxin.AML<-golub[2124,gol.fac=="AML"]

#a) Find the bootstrap 95% CIs for the mean and for the variance of the gene expression in each group separately.
#ALL mean and for the variance
n<-length(Zyxin.ALL) 
nboot<-1000 
boot.xbar <- rep(NA, nboot) 
boot.xbar2 <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin.ALL[sample(1:n,replace=TRUE)] 
  boot.xbar[i]<-mean(data.star)
  boot.xbar2[i]<-var(data.star)
} 
#ALL Mean
quantile(boot.xbar, c(0.025,0.975))
#ALL variance
quantile(boot.xbar2, c(0.025,0.975))


#AML mean and for the variance
n<-length(Zyxin.AML) 
nboot<-1000 
boot.xbar <- rep(NA, nboot) 
boot.xbar2 <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin.AML[sample(1:n,replace=TRUE)] 
  boot.xbar[i]<-mean(data.star)
  boot.xbar2[i]<-var(data.star)
} 
#AML Mean
quantile(boot.xbar, c(0.025,0.975))
#AML variance
quantile(boot.xbar2, c(0.025,0.975))

#b)
#parametric 95% CIs for ALL
n<-length(Zyxin.ALL) #sample size n 
cim.all <- mean(Zyxin.ALL)+qt(c(0.025,0.975),df=n-1)*sd(Zyxin.ALL)/sqrt(n) #95% t-interval 
civ.all <- var(Zyxin.ALL)+qt(c(0.025,0.975),df=n-1)*sd(Zyxin.ALL)/sqrt(n) #95% t-interval 
print(cim.all) #print the 95% CI for ALL patients (stored in vector cim.all) 
print(civ.all)

#parametric 95% CIs for ALL
n<-length(Zyxin.AML) 
cim.aml <- mean(Zyxin.AML)+qt(c(0.025,0.975),df=n-1)*sd(Zyxin.AML)/sqrt(n) 
civ.aml <- var(Zyxin.AML)+qt(c(0.025,0.975),df=n-1)*sd(Zyxin.AML)/sqrt(n) 
print(cim.aml)
print(civ.aml) 

#c)
#Find the bootstrap 95% CI for the median gene expression for ALL
n<-length(Zyxin.ALL) 
nboot<-1000 
boot.xbar <- rep(NA, nboot) 
for (i in 1:nboot) {
  data.star <- Zyxin.ALL[sample(1:n,replace=TRUE)] 
  boot.xbar[i]<-median(data.star)
} 
#ALL Median
quantile(boot.xbar, c(0.025,0.975))


#the bootstrap 95% CI for the median gene expression for AML
n<-length(Zyxin.AML) 
nboot<-1000 
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin.AML[sample(1:n,replace=TRUE)] 
  boot.xbar[i]<-median(data.star)
} 
#AML Median
quantile(boot.xbar, c(0.025,0.975))

#Problem 4)
#Find the 90% CI for lambda = 0.1
rm(list=ls())  #clears environment
x<-0.1
n<-50
nsim<-1000 
#to calculate the mean
mean.x <- rep(NA, nsim) 
#to calculate the variance
var.x <- rep(NA, nsim)
for (i in 1:nsim) {
  data.star <- rpois(n,x) 
  mean.x[i]<-mean(data.star)
  var.x[i]<-var(data.star)
} 
#Mean
quantile(mean.x, c(0.05,0.95))
#Variance
quantile(var.x, c(0.05,0.95))

#This is the analytic calculation
#CI from sample mean
ci.m <- x+qt(c(0.05,0.95),df=n-1)*sqrt(x/n) #90% CI
ci.m
#cI from variance
ci.v <- ((n-1)*x)/ qchisq(c(0.95,0.05),df=n-1)
ci.v

#Find the 90% CI for lambda = 1
rm(list=ls())  #clears environment
x<-1
n<-50
nsim<-1000 
#to calculate the mean
mean.x <- rep(NA, nsim) 
#to calculate the variance
var.x <- rep(NA, nsim)
for (i in 1:nsim) {
  data.star <- rpois(n,x) 
  mean.x[i]<-mean(data.star)
  var.x[i]<-var(data.star)
} 
#Mean
quantile(mean.x, c(0.05,0.95))
#Variance
quantile(var.x, c(0.05,0.95))

#This is the analytic calculation
#CI from sample mean
ci.m <- x+qt(c(0.05,0.95),df=n-1)*sqrt(x/n) #90% CI
ci.m
#cI from variance
ci.v <- ((n-1)*x)/ qchisq(c(0.95,0.05),df=n-1)
ci.v

#Find the 90% CI for lambda = 10
rm(list=ls())  #clears environment
x<-10
n<-50
nsim<-1000 
#to calculate the mean
mean.x <- rep(NA, nsim) 
#to calculate the variance
var.x <- rep(NA, nsim)
for (i in 1:nsim) {
  data.star <- rpois(n,x) 
  mean.x[i]<-mean(data.star)
  var.x[i]<-var(data.star)
} 
#Mean
quantile(mean.x, c(0.05,0.95))
#Variance
quantile(var.x, c(0.05,0.95))

#This is the analytic calculation
#CI from sample mean
ci.m <- x+qt(c(0.05,0.95),df=n-1)*sqrt(x/n) #90% CI
ci.m
#cI from variance
ci.v <- ((n-1)*x)/ qchisq(c(0.95,0.05),df=n-1)
ci.v


