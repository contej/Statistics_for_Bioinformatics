rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio


#Problem 1

#a)
#Hypothesis Testing Process
#1) ALL patients should be greater than -0.9
#H0: mu > -0.9 HA: mu < -0.9

#2) Here we just use the Golub et al. (1999) data set, and take those 
#Gdf5 gene expression values for the ALL patients.
rm(list=ls())  #clears environment

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
H4j_ALL <- golub[2972,gol.fac=="ALL"]

#3) summarize the data and carry out a test
t.test(H4j_ALL, mu=-0.9, alternative = "greater")

#b)
#Hypothesis Testing Process
#1) "H4/j gene" gene expression value in ALL group differs from the mean "H4/j gene"
#H0: mu does not equal H4/j gene AML; HA: mu equals H4/j gene AML

#2) Here we just use the Golub et al. (1999) data set, and take those 
#Gdf5 gene expression values for the AML patients.
H4j_AML <- golub[2972,gol.fac=="AML"]

#3) summarize the data and carry out a test
t.test(H4j_AML,H4j_ALL)

#c)
#Hypothesis Testing Process
#1) "H4/j gene" gene is lower than the mean expression value for the 
#"APS Prostate specific antigen" gene
#H0: mu < APS gene ALL; HA: mu > APS gene ALL

#2) Here we just use the Golub et al. (1999) data set, and take those 
#APS gene expression values for the ALL patients.
APS_ALL <- golub[2989,gol.fac=="ALL"]

#3) summarize the data and carry out a test
t.test(H4j_ALL,APS_ALL, paired=T, alternative = "less" )

#d)
#Hypothesis Testing Process
#1) "H4/j gene" gene is lower than the mean expression value for the 
#"APS Prostate specific antigen" gene
#H0: mu < APS gene; HA: mu > APS gene

#2) Here we just use the Golub et al. (1999) data set, and take  
#APS and H4j gene expression values.
H4j <- golub[2972,]
APS <- golub[2989,]

#3) summarize the data and carry out a test
t.test(H4j,APS, paired=T, alternative = "less" )

#e)
#Hypothesis Testing Process
#1) "H4/j gene" gene is lower than the mean expression value for the 
#"APS Prostate specific antigen" gene
#H0: mu < APS gene; HA: mu > APS gene

#2) Here we just use the Golub et al. (1999) data set, and take those 
#H4j gene expression values.


#3) summarize the data and carry out a test
t.test(H4j, mu = -0.6, alternative = "greater" )

#Problem 2
rm(list=ls())  #clears environment
pbinom(90,2000,0.05)

#Problem 3
#Calculating Power Using Monte Carlo simulation
rm(list=ls())  #clears environment
x.sim<-matrix(rnorm(10000*10, mean=2, sd = 4), ncol=20)
tstat<-function(x) (mean(x)-3)/sd(x)*sqrt(length(x))
tstat.sim<-apply(x.sim,1,tstat) #Calculate t-test statistic for each data set
power.sim<-mean(tstat.sim>qt(c(0.3,0.4),df=19)) #Calculate the rejection rate
#Display rejection rate (power) with its 95% CI
power.sim+c(-1,0,1)*qnorm(0.975)*sqrt(power.sim*(1-power.sim)/10000)


#Problem 4
rm(list=ls())  #clears environment
data(golub, package = "multtest") 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
p.values <- apply(golub[,gol.fac=="ALL"], 1, function(x) t.test(x)$p.value) 
p.bon <-p.adjust(p=p.values, method="bonferroni") 
p.fdr <-p.adjust(p=p.values, method="fdr") 
#ALL Values
sum(p.values<0.05)
sum(p.bon<0.05)
sum(p.fdr<0.05)
o <- order(p.values,decreasing=FALSE)
golub.gnames[o[1:3],2]

p.values <- apply(golub[,gol.fac=="AML"], 1, function(x) t.test(x)$p.value) 
p.bon <-p.adjust(p=p.values, method="bonferroni") 
p.fdr <-p.adjust(p=p.values, method="fdr") 
#AML Values
sum(p.values<0.05)
sum(p.bon<0.05)
sum(p.fdr<0.05)
o <- order(p.values,decreasing=FALSE)
golub.gnames[o[1:3],2]

#b)
pt <- apply(golub, 1, function(x) t.test(x ~ gol.fac)$p.value)
o <- order(pt,decreasing=FALSE)
golub.gnames[o[1:3],2]


#Problem 5
rm(list=ls())  #clears environment

#a)
#declare variables
NCL=0.95 #conf.level
p=0.2
n=40
x=n*p

#Wald
Wald.ci<-function(x,n,conf.level=NCL){
  
  alpha = 1-NCL 
  z <- qnorm(1-alpha/2)
  x<-n*p
  
  p+c(-z,z)*sqrt((p*(1-p))/n)
  
}

Wald.ci(x,n,NCL)


#wilson
wilson.ci<-function(x,n,conf.level=NCL){
  
  alpha = 1-NCL 
  z <- qnorm(1-alpha/2)
  x<-n*p
  
  a<-p+((z^2)/(2*n))
  b<-z*sqrt(((p*(1-p))+((z^2)/(4*n)))/n)
  c<-1+((z^2)/n)
  
  (a+c(-b,b))/c
}

wilson.ci(x,n,NCL)


#Agresti-Coull
AC.ci<-function(x,n,conf.level=NCL){
  
  alpha = 1-NCL 
  z <- qnorm(1-alpha/2)
  x<-n*p
  
  n2=n+(z^2)
  p2=(1/n2)*(x+((z^2)/2))
  
  p2+c(-z,z)*sqrt((p2*(1-p2))/n2)
}

AC.ci(x,n,NCL)

#b)
#Wald Monte Carlo
n=40
p=0.2

x.sim=rbinom(10000,size=n,prob=p)
Wald.sim=matrix(Wald.ci(x.sim,n=n, conf.level=0.95),nrow=2)

mean(Wald.sim[1,]);mean(Wald.sim[2,])


#Wilson Monte Carlo
n=40
p=0.2

x.sim=rbinom(10000,size=n,prob=p)
wilson.sim=matrix(wilson.ci(x.sim,n=n, conf.level=0.95),nrow=2)

mean(wilson.sim[1,]);mean(wilson.sim[2,])


#agresti-coull Monte Carlo
n=40
p=0.2

x.sim=rbinom(10000,size=n,prob=p)
AC.sim=matrix(AC.ci(x.sim,n=n, conf.level=0.95),nrow=2)

mean(AC.sim[1,]);mean(AC.sim[2,])
