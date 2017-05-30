rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#n =5 mean = 5 sd=3 and 1000 samples

my.data.5<-matrix(rnorm(5*1000,mean = 5, sd = 3),nrow=5) 
my.stat.5<-sqrt(5)*(apply(my.data.5,2,mean)-2) 

par(mfrow=c(1,2))
hist(my.stat.5, freq = FALSE, main="histogram of 1000 S",xlab="S", ylim=c(0,0.15))
#simulated distribution
curve(dnorm(x,mean=5,sd=3),add=T) #overlay chi-square(4) density curve
qqplot(qnorm(ppoints(1000), mean=5,sd=3), my.stat.5, 
       main = "Q-Q plot when n=5",xlab="",ylab="S") 
qqline(my.stat.5, distribution = function(p) qnorm(p,mean=5,sd=3), datax=TRUE, col=2)

#problem 2
rm(list=ls())  #clears environment
1 - pnorm(16, mean=14,sd=sqrt(4.2))

#Problem 3
rm(list=ls())  #clears environment
#µ1 = 9; ??1^2 = 3; µ2 = 10; ??2^2 = 9; ??12 = 2
#this is the equation that I worked out in the home work section
require(mvtnorm)
mydata.mat<-as.matrix(rmvnorm(50*1000, mean=c(9,10), sigma=matrix(c(3,2,2,5), nrow=2) ))

#this will be the x column
X<-mydata.mat[,1]
#this will be the y column
Y<-mydata.mat[,2]
#this will make a 3rd column, x+0.5
mydata.mat<-cbind(mydata.mat, mydata.mat[,1]+0.5)
#this is now the 3rd column which is X+0.5
X.5<-mydata.mat[,3]

#this gets P(X + 0.5 < Y)
E1<-sum(as.numeric((Y - X) > 0.5))
PXY<-E1/NROW(mydata.mat)
PXY

#the 95% confidence interval 
mean(Y > X.5)+c(-1,1)*1.96*sqrt(var(Y > X.5)/50)

#Problem 4
rm(list=ls())  #clears environment
#chi where df=10
X1<-rchisq(10000, df=10)
#gamma where alpha=1 and beta=2
X2<-rgamma(10000, shape=1, scale = 2)
#t-distribution where m=3
X3<-rt(10000, df=3)
#calculate mean
y<-sqrt(X1)*X2+4*(X3^2)
mean(y)

#problem 5
rm(list=ls())  #clears environment

n <- 1000
an <- sqrt(2*log(n)) - 0.5*(log(log(n))+log(4*pi))*(2*log(n))^(-1/2)
bn <- (2*log(n))^(-1/2)
e <- double()
for (i in 1:1000) e[i] <- (max(rnorm(n))-an)/bn
par(mfrow=c(1,1))
plot(density(e), ylim=c(0,0.5), lwd=2,main="Question 5")
f<-function(x){exp(-x)*exp(-exp(-x))}
curve(f,range(density(e)$x),add=TRUE,col = "blue", lwd=2)
curve(dnorm,add=TRUE,col = "red", lwd=2)
legend("top",title="Legend", horiz=TRUE, c("Extreme value","Density of generated extreme data", "Normal distribution"), 
       col=c("blue","black","red"),bg="grey96", inset=.02, lty=1, lwd=2, text.width=3, cex=.55)