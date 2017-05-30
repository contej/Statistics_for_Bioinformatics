rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#(a) Conduct the one-way ANOVA. Do the disease stages affect the mean gene expression value?
library(ALL);data(ALL) 
ALLB <- ALL[,ALL$BT %in% c("B","B1","B2","B3","B4")] #patients in five stages 
y <- exprs(ALLB)["109_at",] #exprs function gives gene expression values. We only take gene 1241_at values. 
anova(lm(y ~ ALLB$BT)) #anova. Group indicator in ALLB$BT


#(b) From the linear model fits, find the mean gene expression value among B3 patients.

summary(lm(y ~ ALLB$BT)) #summary of anova fit
#the mean gene expression for B3 patients is 6.8102-0.1249=6.6853

#(c) Which group's mean gene expression value is different from that of group B?
pairwise.t.test(y, ALLB$BT) #Do all pairwise comparison of group means

#(d) Use the pairwise comparisons at FDR=0.05 to find which group means are different. What is your conclusion?
pairwise.t.test(y,ALLB$BT,p.adjust.method='fdr') #FDR-adjusted pairwise tests

#(e) Check the ANOVA model assumptions with diagnostic tests? Do we need to apply robust ANOVA tests here? If yes, apply the appropriate tests and state your conclusion.
shapiro.test(residuals(lm(y ~ ALLB$BT))) #normality test on residuals

library(lmtest)
bptest(lm(y ~ ALLB$BT), studentize = FALSE) #test equal variances

#problem 2
rm(list=ls())  #clears environment

#(a) Use FDR adjustments at 0.05 level. How many genes are expressed differently in some of the groups? 
library(ALL); data(ALL)
library(lmtest)
ALLB <- ALL[,ALL$BT %in% c("B","B1","B2","B3","B4")]

pkw <- apply(exprs(ALLB), 1, function(x) kruskal.test(x ~ ALLB$BT)$p.value)
p.fdr.ALL <-p.adjust(p=pkw, method="fdr") 
sum(p.fdr.ALL<0.05)

#(b) Find the probe names for the top five genes with smallest p-values. 
o.ALL <- order(pkw,decreasing=FALSE)
featureNames(ALLB)[o.ALL[1:5]]

#Problem 3
#(a) Conduct the appropriate ANOVA analysis. Does any of the two factors affects 
#the gene expression values? Are there interaction between the two factors?
rm(list=ls())  #clears environment
library("ALL"); data(ALL) 
library(lmtest)
ALLBm <- ALL[,which(ALL$BT %in% c("B1","B2","B3","B4") & ALL$sex %in% c("F","M"))] #select patients 
y<-exprs(ALLBm)["38555_at",] #gene 32069_at expression values 
Bcell<-ALLBm$BT # B-cell stages 
sex<-ALLBm$sex #molecular biology types 
anova(lm(y~ Bcell*sex)) #full two-way ANOVA
anova(lm(y~ Bcell+sex)) #Additive (no interaction) two-way ANOVA
summary(lm(y~ Bcell+sex)) #summary of anova fit

#(b) Check the ANOVA model assumption with diagnostic tests? Are any of the assumptions violated? 
shapiro.test(residuals(lm(y~ Bcell+sex)))

bptest(lm(y ~ (Bcell+sex)), studentize = FALSE) #test equal variances

#Problem 4
rm(list=ls())  #clears environment
data(ALL,package="ALL");library(ALL)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")] #patients in 3 stages
data<- exprs(ALLB123)["1242_at",] #gene 1242_g_at expression values
group<-ALLB123$BT[,drop=T] #drop unused levels, keep only B1,B2,B3
n<-length(data) #sample size n
n.group<-length(by(data,group,mean)) #this is g
group.means<-mean(by(data,group,mean)) #this is uj
T.obs<-(1/(n.group-1))*sum((group.means-mean(group.means))^2)
n.perm=2000 # we will do 2000 permutations
T.perm = rep(NA, n.perm) #A vector to save permutated statistic
for(i in 1:n.perm) {
  data.perm = sample(data, n, replace=F) #permute data
  group.means.perm<-mean(by(data.perm,group,mean)) #permute uj
  n.group.perm<-length(by(data.perm,group,mean)) #permute g
  T.perm[i] = (group.means.perm-mean(group.means.perm)^2)/(n.group.perm-1) #Permuted statistic
}
mean(T.perm>=T.obs) #p-value


