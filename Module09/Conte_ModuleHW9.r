rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
data(golub,package="multtest") #Load golub data from multtest package 
#The row number of the GRO2 gene. 
grep("GRO2 GRO2 oncogene",golub.gnames[,2])
#The row number of the GRO3 gene. 
grep("GRO3 GRO3 oncogene",golub.gnames[,2])

#(a) Find the correlation between the expression values of these two genes.
x <- golub[2714,] #GRO2 gene.
y <- golub[2715,] #GRO3 gene.
cor(x,y) #calculate correlation

plot (x,y)

#(b) Find the parametric 90% confident interval for the correlation with cor.test().
#use ?cor.test to learn how to set the confidence level different from the default 
#value of 95%.
?cor.test 
#set conf.level=0.9 to make it 90 percent confidence interval
cor.test(x,y, conf.level=0.9) #test if true correlation=0

#(c) Find the bootstrap 90% confident interval for the correlation.
nboot <- 2000 # We will resample 2000 times
boot.cor <- matrix(0,nrow=nboot, ncol = 1) #A vector to save the resampled statistics
data <- cbind(x,y) #Data set with x and y in two columns.
for (i in 1:nboot){
  dat.star <- data[sample(1:nrow(data),replace=TRUE), ] #Resample the
  pairs
  boot.cor[i,] <- cor(dat.star[,1], dat.star[,2]) #Correlation on resampled
  data
}
quantile(boot.cor[,1],c(0.05,0.95)) #Find quantiles for resampled statistics

#(d) Test the null hypothesis that correlation = 0.64
n<-length(x) #sample size n = number of pairs
T.obs<- 0.64 #correlation = 0.64
n.perm=2000 # We will permute 2000 times
T.perm = rep(NA, n.perm) #A vector to save the permuted statistics
for(i in 1:n.perm) {
  x.perm = sample(x, n, replace=F) #permute data (x only)
  T.perm[i] = cor(x.perm, y) #Permuted statistic is the correlation
}
mean(abs(T.perm)>abs(T.obs)) #p-value 

#problem 2
#On the Golub et al. (1999) data set, we consider the correlation between the Zyxin 
#gene expression values and each of the gene in the data set. 
rm(list=ls())  #clears environment
data(golub, package="multtest")
Zyxin <- (golub[2124,])
n <- 3051
t.perm <- rep(NA, n) #A vector to save the permuted statistics
for (i in 1:n){
  y.perm = golub[i,]
  t.perm[i] = cor(Zyxin, y.perm) #Permuted statistic is the correlation
}
#(a) How many of the genes have correlation values less than negative 0.5? 
#(Those genes are highly negatively correlated with Zyxin gene). 
sum(t.perm<(-0.5))

#(b) Find the gene names for the top five genes that are most negatively 
#correlated with Zyxin gene.
o <- order(t.perm,decreasing=FALSE)
golub.gnames[o[1:5],2]

#(c) Using the t-test, how many genes are negatively correlated with the Zyxin gene?
p.values <- apply(golub, 1, function(x) cor.test(Zyxin,x)$p.value) 
p.fdr <-p.adjust(p=p.values, method="fdr") 
sum(p.fdr<0.05)

#problem 3
#On the Golub et al. (1999) data set, regress the expression values for the GRO3 
#GRO3 oncogene on the expression values of the GRO2 GRO2 oncogene.
rm(list=ls())  #clears environment
data(golub,package="multtest") #Load golub data from multtest package 
GRO2_d3 <- golub[2714,] #GRO2 gene.
GRO3_x6 <- golub[2715,] #GRO3 gene.

#(a) Is there a statistically significant linear relationship between the two genes' 
#expression?
cor.test(GRO2_d3,GRO3_x6) #test if true correlation=0

#(b)Test if the slope parameter is less than 0.5 at the ?? = 0.05 level. 
reg.fit<-lm(GRO3_x6 ~ GRO2_d3) #Regression GRO2_x6 = b0+ b1*GRO3_d3
reg.fit #Results of the regression fit

#b0 = -0.8426 and b1 = 0.3582
summary(reg.fit) #summary of regression results

#To get the 90% confidence intervals for b0 and b1
confint(reg.fit, level=0.9) #Show 90% 2-sided CIs from regression fit

#(c)Find an 80% prediction interval for the GRO3 GRO3 oncogene expression when 
#GRO2 GRO2 oncogene is not expressed (zero expression value).
predict(reg.fit, newdata=data.frame(GRO2_d3=0), interval="prediction", level=0.8)

#(d) Check the regression model assumptions. 
shapiro.test(residuals(reg.fit)) #normality test on residuals
plot(reg.fit,which=1)
plot(reg.fit,which=2)

#problem 4
#For this problem, work with the data set stackloss that comes with R. 
#You can get help on the data set with ?stackloss command.
rm(list=ls())  #clears environment
?stackloss
#(a) Regress stack.loss on the other three variables. 
#What is the fitted regression equation?
lin.reg<-lm(stack.loss~Air.Flow+Water.Temp+Acid.Conc., data=stackloss) #multiple regression of stack.loss on 3 variables
summary(lin.reg) #summary of regression results

#(c) Find a 90% confidence interval and 90% prediction interval for stack.loss when Air.Flow=60, Water.Temp=20 and Acid.Conc.=90.
#confidence interval
predict(lin.reg, newdata=data.frame(Air.Flow=60, Water.Temp=20, Acid.Conc.=90), interval="confidence", level=0.9)
#prediction interval
predict(lin.reg, newdata=data.frame(Air.Flow=60, Water.Temp=20, Acid.Conc.=90), interval="prediction", level=0.9)
