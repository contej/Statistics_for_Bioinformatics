rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#(a) Define an indicator variable IsB such that IsB=TRUE for B-cell patients 
#and IsB=FALSE for T-cell patients. 
#(b) Use two genes "39317_at" and "38018_g_at" to fit a classification tree for 
#IsB. Print out the confusion matrix. Plot ROC curve for the tree. 
library("ALL"); data(ALL); #load package and data 
library(ROCR) ##Load library
allB <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))] #select patients 
allBT <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4","T","T1","T2","T3","T4"))] #select patients
prob.name <- c("39317_at","38018_g_at") 
expr.data <-exprs(allBT)[prob.name,]
IsB <- (ALL$BT==allB$BT) #A boolean class indicator (in B)
data.lgr <- data.frame( IsB, t(expr.data)) #A data frame of class indicator and expression data. 
fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.lgr) #logistic regression of class indicator on other variables (1389_at expression values) 
pred.prob <- predict(fit.lgr, data=data.lgr$expr.data, type="response") #predict class probability (response) using the regression fit on expression data. 
pred.B<- factor(pred.prob> 0.5, levels=c(TRUE,FALSE), labels=c("B","not")) #classify according to prediction probability, and label TRUE to 'B' and FALSE to 'not' 
IsB<-factor(IsB, levels=c(TRUE,FALSE), labels=c("B","not")) #class indicator label to 'B' and 'not' 
table(pred.B, IsB) #confusion matrix that compares predicted versus true classes

#Classification Tree on ALL Data
library("hgu95av2.db");
ALLB <- ALL[,ALL$BT %in% c("B","B1","B2","B3","B4")] #patients in B1,B2,B3 stages
pano <- apply(exprs(ALLB), 1, function(x) anova(lm(x ~ ALLB$BT))$Pr[1]) #do ANOVA on every gene 
names <- featureNames(ALL)[pano<0.000001] #get probe names for genes passing the ANOVA gene filter 
symb <- mget(names, env = hgu95av2SYMBOL) #get the gene names 
ALLBTnames <- ALLB[names, ] #keep only the gene passing filter 
probedat <- as.matrix(exprs(ALLBTnames)) #save expression values in probedat 
row.names(probedat)<-unlist(symb) #change row names to gene names
require(rpart) #load library 
B.stage <- factor(ALLBTnames$BT) #The true classes: In all B stages 
c.tr <- rpart(B.stage ~ ., data = data.frame(t(probedat))) #Fit the tree on data set 'probedat'. The transpose is needed since we are classifying the patients (rows). 
plot(c.tr, branch=0,margin=0.1); #Plot the tree with V-shaped branch (=0), leave margins for later text additions. 
text(c.tr, digits=3) #Add text for the decision rules on the tree. use 3 digits for numbers 
rpartpred <- predict(c.tr, type="class") #predicted classes from the tree 
table(rpartpred, B.stage) #confusion matrix compare predicted versus true classes

#ROC curve for the tree
pred.prob <- predict(fit.lgr, data=data.lgr[,-1], type="response")
pred <- prediction(pred.prob, IsB) 
perf <- performance(pred, "tpr", "fpr" ) 
plot(perf) #add green ROC curve for this classifier


#(c) Find its empirical misclassification rate (mcr), false negative rate (fnr) 
#and specificity. Find the area under curve (AUC) for the ROC curve.
#empirical misclassification rate
#Training and Validation (testing) 
set.seed(131) #Set an arbitrary random seed. 
testID <- sample(1:128, 51, replace = FALSE) #randomly select 31 test (validation) cases out of 78 
data.tr<-data.lgr[-testID, ] #training data, remove test cases with "-" 
data.test<-data.lgr[testID, ] #test (validation) data 
fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.tr) #fit logistic regression on training data. 
pred.prob <- predict(fit.lgr, data=data.tr, type="response") #training data prediction 
pred.B<- factor(pred.prob> 0.5)
mcr.tr<- sum(pred.B!=data.tr$IsB)/length(data.tr$IsB) #mcr on training data. (# of misclassification)/(# total cases) 
pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #test data prediction 
pred.B<- factor(pred.prob> 0.5) 
mcr.test<- sum(pred.B!=data.test$IsB)/length(data.test$IsB) #mcr on test data 
data.frame(mcr.tr, mcr.test) #print out mcr

#area under curve (AUC) for the ROC curve
performance(pred,"auc")


#(d) Use 10-fold cross-validation to estimate its real false negative rate (fnr). What is your estimated fnr? 
## 10 fold cross-validation 
require(caret) 
n<-dim(data.lgr)[1] #size of the whole sample 
index<-1:n #index for data points 
K<-10 #number of folds 
flds <- createFolds(index, k=K) #create K folds 
mcr.cv.raw<-rep(NA, K) ## A vector to save raw mcr 
for (i in 1:K) { 
  testID<-flds[[i]] #the i-th fold as validation set 
  data.tr<-data.lgr[-testID,] #remove the test cases from training data 
  data.test<-data.lgr[testID,] #validation (test) data 
  fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.tr) #train model 
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #prediction probability 
  pred.B<- (pred.prob> 0.5) #prediction class 
  mcr.cv.raw[i]<- sum(pred.B!=data.test$IsB)/length(pred.B)#mcr on testing case 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over K folds.
mcr.cv

#(e) Do a logistic regression, using genes "39317_at" and "38018_g_at" to predict 
#IsB. Find an 80% confidence interval for the coefficient of gene "39317_at".
allB <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))] #select patients 
allBT <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4","T","T1","T2","T3","T4"))] #select patients
prob.name <- c("39317_at","38018_g_at") 
expr.data <-exprs(allBT)[prob.name,]
IsB <- (ALL$BT==allB$BT) #A boolean class indicator (in B)
data.lgr <- data.frame( IsB, t(expr.data)) #A data frame of class indicator and expression data. 
fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.lgr) #logistic regression of class indicator on other variables (1389_at expression values) 
pred.prob <- predict(fit.lgr, data=data.lgr$expr.data, type="response") #predict class probability (response) using the regression fit on expression data. 
pred.B<- factor(pred.prob> 0.5, levels=c(TRUE,FALSE), labels=c("B","not")) #classify according to prediction probability, and label TRUE to 'B' and FALSE to 'not' 
IsB<-factor(IsB, levels=c(TRUE,FALSE), labels=c("B","not")) #class indicator label to 'B' and 'not' 
table(pred.B, IsB) #confusion matrix that compares predicted versus true classes
confint(fit.lgr, level=0.8) #Find 80% 2-sided CIs for the parameters

#f) Use n-fold cross-validation to estimate misclassification rate (mcr) of the 
#logistic regression classifier. What is your estimated mcr?
n<-dim(data.lgr)[1] #size of the whole sample 
index<-1:n #index for data points 
K<-n #number of folds 
flds <- createFolds(index, k=K) #create K folds 
mcr.cv.raw<-rep(NA, K) ## A vector to save raw mcr 
for (i in 1:K) { 
  testID<-flds[[i]] #the i-th fold as validation set 
  data.tr<-data.lgr[-testID,] #remove the test cases from training data 
  data.test<-data.lgr[testID,] #validation (test) data 
  fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.tr) #train model 
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #prediction probability 
  pred.B<- (pred.prob> 0.5) #prediction class 
  mcr.cv.raw[i]<- sum(pred.B!=data.test$IsB)/length(pred.B)#mcr on testing case 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over K folds.
mcr.cv

#(g) Conduct a PCA on the scaled variables of the whole ALL data set (NOT just 
#the two genes used above). 
pca.ALL<-prcomp(exprs(ALL), scale=TRUE) #Apply PCA on ALL data, scale each variable first
summary(pca.ALL) #print out summary of the PCA fit

#(h) Do a SVM classifier of IsB using only the first five PCs. (The number K=5 
#is fixed so that we all use the same classifier. You do not need to choose this 
#number in the previous part (g).) What is the sensitivity of this classifier?

B.stage<-factor(ALLBTnames$BT[1:5]) #the classes 
data.pca<-pca.ALL$x[,1:5] #keep only the first five principal components
n<-length(B.stage)

pca.ALL.svm <- svm(data.pca, B.stage, type = "C-classification", kernel = "linear") #train SVM 
svmpred <- predict(pca.ALL.svm , data.pca) #get SVM prediction. 
mcr.svm<- mean(svmpred!=B.stage) #misclassification rate 
### leave-one-out cross validation 
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation 
for (i in 1:n) { 
  svmest <- svm(data.pca[-i,], B.stage[-i], type = "C-classification", kernel = "linear") #train SVM without i-th observation 
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here transpose t() is used to make the vector back into a 1 by ncol matrix 
  mcr.cv.raw[i]<- mean(svmpred!=B.stage[i]) #misclassification rate 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.svm, mcr.cv)

#(i) Use leave-one-out cross-validation to estimate misclassification rate (mcr) 
#of the SVM classifier. Report your estimate.

pca.ALL<-data.frame(B.stage, data.pca) #combine response variable with PCA data 
fit <- rpart(B.stage ~ ., data = pca.ALL, method = "class") #Fit tree pca.ALL data 
pred.tr<-predict(fit, pca.ALL, type = "class") #predict classes from the tree 
mcr.tr <- mean(pred.tr!=B.stage) #misclassification rate ### leave-one-out cross validation 
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation 
for (i in 1:n) { 
  fit.tr <- rpart(B.stage ~ ., data = pca.ALL[-i,], method = "class") #train the tree without i-th observation 
  pred <- predict(fit.tr, pca.ALL[i,], type = "class")#use tree to predict i-th observation class 
  mcr.cv.raw[i]<- mean(pred!=B.stage[i]) #check misclassifion 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.tr, mcr.cv)

#Problem 2
rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

## 10 fold cross-validation 
pca.iris<-prcomp(iris[,1:4], scale=TRUE) #Apply PCA on first four variables in iris data, scale each variable first
Species<-iris$Species #response variable with true classes 
data.pca<-pca.iris$x[,1:3] #keep only the first three principal components 
n<-length(Species) #sample size n
iris2<-data.frame(Species, data.pca) #combine response variable with PCA data
n<-dim(iris2)[1] #size of the whole sample 
index<-1:n #index for data points 

#k=1
K<-1 #number of folds 
flds <- createFolds(index, k=K) #create K folds 
mcr.cv.raw<-rep(NA, K) ## A vector to save raw mcr 
for (i in 1:K) { 
  testID<-flds[[i]] #the i-th fold as validation set 
  data.tr<-iris2[-testID,] #remove the test cases from training data 
  data.test<-iris2[testID,] #validation (test) data 
  fit.lgr <- glm(Species ~ ., family=binomial(link='logit'),data = data.tr) #train model 
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #prediction probability 
  pred.B1<- (pred.prob> 0.5) #prediction class 
  mcr.cv.raw[i]<- sum(pred.B1!=data.pca)/length(pred.B1)#mcr on testing case 
} 
mcr.k<-mean(mcr.cv.raw) #average the mcr over K folds.

  ### leave-one-out cross validation 
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation 
for (i in 1:n) { 
  svmest <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel = "linear") #train SVM without i-th observation 
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here transpose t() is used to make the vector back into a 1 by ncol matrix 
  mcr.cv.raw[i]<- mean(svmpred!=Species[i]) #misclassification rate 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.k, mcr.cv)

#k=2
K<-2 #number of folds 
flds <- createFolds(index, k=K) #create K folds 
mcr.cv.raw<-rep(NA, K) ## A vector to save raw mcr 
for (i in 1:K) { 
  testID<-flds[[i]] #the i-th fold as validation set 
  data.tr<-iris2[-testID,] #remove the test cases from training data 
  data.test<-iris2[testID,] #validation (test) data 
  fit.lgr <- glm(Species ~ ., family=binomial(link='logit'),data = data.tr) #train model 
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #prediction probability 
  pred.B1<- (pred.prob> 0.5) #prediction class 
  mcr.cv.raw[i]<- sum(pred.B1!=data.pca)/length(pred.B1)#mcr on testing case 
} 
mcr.k<-mean(mcr.cv.raw) #average the mcr over K folds.

### leave-one-out cross validation 
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation 
for (i in 1:n) { 
  svmest <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel = "linear") #train SVM without i-th observation 
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here transpose t() is used to make the vector back into a 1 by ncol matrix 
  mcr.cv.raw[i]<- mean(svmpred!=Species[i]) #misclassification rate 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.k, mcr.cv)

#k=3
K<-3 #number of folds 
flds <- createFolds(index, k=K) #create K folds 
mcr.cv.raw<-rep(NA, K) ## A vector to save raw mcr 
for (i in 1:K) { 
  testID<-flds[[i]] #the i-th fold as validation set 
  data.tr<-iris2[-testID,] #remove the test cases from training data 
  data.test<-iris2[testID,] #validation (test) data 
  fit.lgr <- glm(Species ~ ., family=binomial(link='logit'),data = data.tr) #train model 
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #prediction probability 
  pred.B1<- (pred.prob> 0.5) #prediction class 
  mcr.cv.raw[i]<- sum(pred.B1!=data.pca)/length(pred.B1)#mcr on testing case 
} 
mcr.k<-mean(mcr.cv.raw) #average the mcr over K folds.

### leave-one-out cross validation 
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation 
for (i in 1:n) { 
  svmest <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel = "linear") #train SVM without i-th observation 
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here transpose t() is used to make the vector back into a 1 by ncol matrix 
  mcr.cv.raw[i]<- mean(svmpred!=Species[i]) #misclassification rate 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.k, mcr.cv)

#k=4
K<-4 #number of folds 
flds <- createFolds(index, k=K) #create K folds 
mcr.cv.raw<-rep(NA, K) ## A vector to save raw mcr 
for (i in 1:K) { 
  testID<-flds[[i]] #the i-th fold as validation set 
  data.tr<-iris2[-testID,] #remove the test cases from training data 
  data.test<-iris2[testID,] #validation (test) data 
  fit.lgr <- glm(Species ~ ., family=binomial(link='logit'),data = data.tr) #train model 
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response") #prediction probability 
  pred.B1<- (pred.prob> 0.5) #prediction class 
  mcr.cv.raw[i]<- sum(pred.B1!=data.pca)/length(pred.B1)#mcr on testing case 
} 
mcr.k<-mean(mcr.cv.raw) #average the mcr over K folds.

### leave-one-out cross validation 
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation 
for (i in 1:n) { 
  svmest <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel = "linear") #train SVM without i-th observation 
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here transpose t() is used to make the vector back into a 1 by ncol matrix 
  mcr.cv.raw[i]<- mean(svmpred!=Species[i]) #misclassification rate 
} 
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.k, mcr.cv)
