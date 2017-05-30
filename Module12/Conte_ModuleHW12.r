rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#(a) Define an indicator variable ALL.fac such that ALL.fac=1 for 
#T-cell patients and ALL.fac=2 for B-cell patients.
library(ALL); data(ALL)
#remove numbers from the letters
B_T<-sub('\\d+', '', ALL$BT) 
#this replaces T with 1 and B with 2
ALL.fac <- factor(B_T,levels=c("T","B"), labels=c(1,2)) 
ALL.fac

#(b) Plot the histograms for the first three genes' expression values in one row.
par(mfrow=c(1,3))
for (i in 1:3) {
  hist(exprs(ALL)[,i], main="", nclass=15, freq=FALSE, xlab=featureNames(ALL)[i])
}

#(c) Plot the pairwise scatterplots for the first five genes.
pairs(exprs(ALL)[,1:5])

#(d) Do a 3D scatterplot for the genes "39317_at", "32649_at" and "481_at", 
#and color according to ALL.fac (give different colors for B-cell versus 
#T-cell patients). Can the two patient groups be distinguished using these three genes?

require(scatterplot3d) #load the library 'scatterplot3d'

x = (exprs(ALL)["39317_at",])
y = (exprs(ALL)["32649_at",])
z = (exprs(ALL)["481_at",])

par(mfrow=c(1,1))

scatterplot3d(x,y,z, xlab="39317_at", ylab="32649_at", zlab="481_at", 
              color=as.numeric(ALL.fac))

#(e)Do K-means clustering for K=2 and K=3 using the three genes in (d).
clusdata <- data.frame(exprs(ALL)["39317_at",], exprs(ALL)["32649_at",], exprs(ALL)["481_at",])
colnames(clusdata)<- c("39317_at", "32649_at","481_at")
cl.2mean <- kmeans(clusdata, centers=2, nstart = 10)
cl.2mean
cl.3mean <- kmeans(clusdata, centers=3, nstart = 10)
cl.3mean

table(cl.2mean$cluster,ALL.fac)
table(cl.3mean$cluster,ALL.fac)

plot(clusdata[c("39317_at", "32649_at","481_at")], col = cl.2mean$cluster, pch="+") 
points(clusdata[c("39317_at", "32649_at","481_at")], col = ALL.fac, pch="o") 
legend("topright",c("T","B"), pch= c("+","o"))

plot(clusdata[c("39317_at", "32649_at","481_at")], col = cl.3mean$cluster, pch="+") 
points(clusdata[c("39317_at", "32649_at","481_at")], col = ALL.fac, pch="o") 
legend("topright",c("T","B"), pch= c("+","o"))


#(f) Carry out the PCA on the ALL data set with scaled variables. 
#What proportion of variance is explained by the first principal 
#component? By the second principal component?
pr.ALL <- prcomp(exprs(ALL), scale=TRUE) #save (scaled) PCA fit on ALL
pr.ALL
summary(pr.ALL)

#print(t(pr.ALL$rotation[,1:2]), digits=3) #print out the loadings for 
#the first two PCs, transpose so that the loadings are shown in rows 
#instead of columns, display only 3 digits for each entry

#(g) Do a biplot of the first two principal components. Observe the pattern 
#for the loadings. What info is the first principal component summarizing?
biplot(pr.ALL, xlim=c(-0.06,0.03), ylim=c(-0.055,0.055), cex=0.5) #biplot 
#display of the PCA fit, set ranges of x and y, use smaller font (cex=0.5)


#(h) For the second principal component PC2, print out the three genes with 
#biggest PC2 values and the three genes with smallest PC2 values.

ob <- order(pr.ALL$x[,2], decreasing = FALSE)
ob3<-ob[1:3]
ob3

os <- order(pr.ALL$x[,2], decreasing = TRUE)
os3<-os[1:3]
os3

#(i) Find the gene names and chromosomes for the gene with biggest PC2 value 
#and the gene with smallest PC2 value. (Hint: review Module 10 on searching 
#the annotation.)

featureNames(ALL)

#Problem 2
#(a) Create a data set consisting of the first four numerical variables in the iris 
#data set (That is, to drop the last variable Species which is categorical). Then 
#make a scaled data set that centers each of the four variables (columns) to have 
#mean zero and variance one.
rm(list=ls())  #clears environment
library(datasets)
data(iris)

x<-iris[,1:4]
scaled.x <- scale(x)

# check that we get mean of 0 and sd of 1
colMeans(scaled.x)  # faster version of apply(scaled.dat, 2, mean)
apply(scaled.x, 2, sd)

#(b) Calculate the correlations between the columns of the data sets using the 
#cor() function. Show that these correlations are the same for scaled and the 
#unscaled data sets.
cor(x)
cor(scaled.x)

#(c) Calculate the Euclidean distances between the columns of the scaled data set 
#using dist() function. Show that the squares of these Euclidean distances are 
#proportional to the (1-correlation)s. What is the value of the proportional 
#factor here?
#calculation of 1-correlation
1-cor(scaled.x)
#Transpose Iris data set to get the columns into the rows for dist() to function properly
t.scaled.x <- t(scaled.x)
#Euclidean distances between the columns
e.d<-dist(t.scaled.x,method="euclidian")
#squares of these Euclidean distances
(e.d)^2
#sample calculations showing the proportional factor = 298
333.03580/1.1175698
54.25354/0.1820589

#(d) Show the outputs for doing PCA on the scaled data set and on the unscaled 
#data set. (Apply PCA on the two data sets with option "scale=FALSE". Do NOT 
#use option "scale=TRUE", which will scale data no matter which data set you 
#are using.) Are they the same?
pr.scaled.x <- prcomp(scaled.x, scale=FALSE) #save (scaled) PCA fit on ALL
pr.scaled.x

pr.x <- prcomp(x, scale=FALSE) #save (scaled) PCA fit on ALL
pr.x

#(e) What proportions of variance are explained by the first two principle 
#components in the scaled PCA and in the unscaled PCA?
summary(pr.scaled.x)
summary(pr.x)


#(f) Find a 90% confidence interval on the proportion of variance explained 
#by the second principal component, in the scaled PCA.
data <- pr.scaled.x$x; p <- ncol(data); n <- nrow(data) ; nboot<-1000 #define quantities 
sdevs <- array(dim=c(nboot,p)) #empty nboot by p matrix to save sdev for each PC in bootstrapped samples 
for (i in 1:nboot) {
  dat.star <- data[sample(1:n,replace=TRUE),] #bootstrap data 
  sdevs[i,] <- prcomp(dat.star)$sdev #sdev for PCs in bootrap data 
} 
print(names(quantile(sdevs[,1], c(0.05,0.95)))) #display "5%" "95%" 
for (j in 1:p) cat(j, as.numeric(quantile(sdevs[,j], c(0.05,0.95))),"\n") #for j-th PC, print "j" then the 90% CI limits in one line. "\n" starts a new line


a<-0.73+ 0.229+ 0.0367+ 0.00518
a
0.73+ 0.229
0.925+ 0.0531
