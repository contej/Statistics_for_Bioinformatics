rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#a)
#ALL Values
data(golub, package = "multtest") 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
p.values <- apply(golub, 1, function(x) wilcox.test(x~gol.fac, alternative="greater")$p.value)
p.fdr <-p.adjust(p=p.values, method="fdr")
sum(p.fdr<0.05)


#AML Values
p.values.AML <- apply(golub[,gol.fac=="AML"], 1, function(x) wilcox.test(x)$p.value) 
p.fdr.AML <-p.adjust(p=p.values.AML, method="fdr") 
sum(p.fdr.AML<0.05)


#b)
#ALL Names
o.ALL <- order(p.values.ALL,decreasing=FALSE)
golub.gnames[o.ALL[1:3],2]

#AML Names
o.AML <- order(p.values.AML,decreasing=FALSE)
golub.gnames[o.AML[1:3],2]

#Problem 2
#AML Values
p.values.AML <- apply(golub[,gol.fac=="AML"], 1, function(x) shapiro.test(x)$p.value) 
p.fdr.AML <-p.adjust(p=p.values.AML, method="fdr") 
#pass
sum(p.fdr.AML<0.05)
#do not pass
sum(p.fdr.AML>0.05)

#Problem 3
HOXA9<- golub[1391,gol.fac=="ALL"]
CD33<- golub[808,gol.fac=="ALL"]
data(golub, package='multtest') 
wilcox.test (x= HOXA9, y= CD33, paired=T, alternative="two.sided")

#Problem 4
rm(list=ls())  #clears environment

data(UCBAdmissions)
#where [x,,] = Admit [,x,] = gender [,,x] = dept
for ( i in 1:6 ) { # for each department
  info <- UCBAdmissions[,,i] # that department's data
  m <- info[1,1] / (info[1,1]+info[2,1]) # Men's admission rate
  w <- info[1,2] / (info[1,2]+info[2,2]) # Women's admission rate
  
  print ( c ( m, w ) ) # print them
}


 #Problem 5
 rm(list=ls())  #clears environment
 data(golub, package = "multtest") 
 gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
 data<-golub[808,]
 n<-length(data)
 T.obs<-abs(var(data[gol.fac=="ALL"]) / var(data[gol.fac=="AML"]))
 #Observed statistic
 n.perm=2000
 T.perm = rep(NA, n.perm)
 for(i in 1:n.perm) {
   data.perm = sample(data, n, replace=F) #permute data
   T.perm[i] = abs(mean(data.perm[gol.fac=="ALL"])-
                     mean(data.perm[gol.fac=="AML"])) #Permuted statistic
 }
 mean(T.perm>=T.obs) #p-value
