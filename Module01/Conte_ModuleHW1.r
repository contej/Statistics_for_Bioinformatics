# Module 1 Homework

# Problem 1
# (a) What is the class of the object defined be vec <-c(5,TRUE) ? 
rm(list = ls())
vec <-c(5,TRUE) 
# This returns a numeric value

# (b) Suppose I have vectors x <- 1:4 and y <- 1:2. What is the result of the expression x + y?
rm(list = ls())
x <- 1:4 
y <- 1:2
x + y

#(c) Suppose I define the following function in R: 
#fsin<-function(x) sin(pi*x) 
#What will be returned by fsin(1) ? 
rm(list = ls())
fsin<-function(x) sin(pi*x)
fsin(1)

#(d) What is returned by the R command c(1,2) %*% t(c(1,2)) ? 
rm(list = ls())
c(1,2) %*% t(c(1,2))

#(e) Suppose I define the following function in R: 
#Consider the following function: 
#  f <- function(x) { 
#    g <- function(y) { 
#      y + z 
#    } 
#    z <- 4 
#    x + g(x) 
#  } 
#If I then run in R the following statements 
#z <- 15 
#f(3) 
#What value is returned? 
rm(list = ls())
f <- function(x) { 
  g <- function(y) { 
    y + z 
  } 
  z <- 4 
  x + g(x) 
} 

z <- 15 
f(3) 

# Problem 2
rm(list = ls())
x <- 1:1000
sum(x^2)

#Question 3
#a)Create a vector X of length 20, with the kth element in X = 2k, 
#for k=1...20. Print out the values of X.
rm(list = ls())
k<- 1:20
X=2*k
X[1:20]

#b)Create a vector Y of length 20, with all elements in Y equal to 0. Print
#out the values of Y
rm(list = ls())
Y<- rep(0,20)
Y

# c) Using a for loop, reassigns the value of the k-th element in Y, for k
# = 1...20. When k < 12, the kth element of Y is reassigned as the cosine
# of (3k). When the k greater or equal to 12, the kth element of Y is reassigned as the
# value an integral
rm(list = ls())
Y <- c()

for(k in 1:20) {
  
  if(k<12){
    Y[k] <- cos(3*k)
    
  }else{
    Y[k] <- integrate(sqrt, 0, k)$value
  }
  
}

Y