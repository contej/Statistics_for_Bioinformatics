rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#HW11
#Problem 1
#(a)	Conduct hierarchical clustering using single linkage and Ward linkage. 
#Plot the cluster dendrogram for both fit.
data(golub, package="multtest") 
clusdata <- data.frame(golub[1042,]) 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 
hcALL.sing<-hclust(dist(clusdata,method="euclidian"),method="single")
hcALL.ward<-hclust(dist(clusdata,method="euclidian"),method="ward.D2")
par(mfrow=c(1,2)) 
plot(hcALL.sing, labels=gol.fac)
plot(hcALL.ward, labels=gol.fac)


#Use function table() to compare the clusters with the two patient groups ALL/AML
sing.table = cutree(hcALL.sing,2)
table(sing.table, gol.fac)

ward.table = cutree(hcALL.ward,2)
table(ward.table, gol.fac)

#(b)Use k-means cluster analysis to get two clusters. 
#Use table() to compare the two clusters with the two patient groups ALL/AML.
colnames(clusdata)<- c("CCND3 Cyclin D3")
cl.2mean <- kmeans(clusdata, centers=2, nstart = 10)
table(cl.2mean$cluster,gol.fac)

#(d)Find the two cluster means from the k-means cluster analysis. 
#Perform a bootstrap on the cluster means. 
cl.2mean

initial <-cl.2mean$centers 
n <- dim(clusdata)[1]; nboot<-1000 
boot.cl <- matrix(NA,nrow=nboot,ncol = 2) 
for (i in 1:nboot){ 
  dat.star <- clusdata[sample(1:n,replace=TRUE),] 
  cl <- kmeans(dat.star, initial, nstart = 10) 
  boot.cl[i,] <- c(cl$centers) 
}
apply(boot.cl,2,mean)
quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))

#(e)Produce a plot of K versus SSE, for K=1, ., 30. 
K<-(1:30); sse<-rep(NA,length(K)) 
for (k in K) { 
  sse[k]<-kmeans(clusdata, centers=k,nstart = 10)$tot.withinss 
} 
par(mfrow=c(1,1)) 
plot(K, sse, type='o', xaxt='n'); axis(1, at = K, las=2)

#Problem 2
#(a)	Select the oncogenes and antigens from the Golub data
rm(list=ls())  #clears environment
library(multtest);data(golub);
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 
oncogene <- grep("oncogene",golub.gnames[,2])
antigen <- grep("antigen",golub.gnames[,2])

#(b) On the selected data, do clustering analysis for the genes 
#(not for the patients). Using K-means and K-medoids with K=2 
#to cluster the genes. Use table() to compare the resulting two 
#clusters with the two gene groups oncogenes and antigens for 
#each of the two clustering analysis.
#K-means
clusdata <- rbind((matrix(oncogene,ncol=2)), (matrix(antigen,ncol=2))) 
colnames(clusdata)<- c("oncogene","antigen")
cl.2mean <- kmeans(clusdata, 2)
table.2mean<-table(cl.2mean$cluster,cl.2mean$cluster)
table.2mean

#K-medoids
library(cluster)
cl.pam <- pam(clusdata, k=2)
table.pam<-table(cl.pam$clustering,cl.pam$clustering)
table.pam

#(c)	Use appropriate tests (from previous modules) to test the 
#marginal independence in the two by two tables in (b). Which 
#clustering method provides clusters related to the two gene groups? 
#To test the marginal independence, we use the chi-square and fisher tests.
#independence test for kmean
chisq.test(table.2mean)
fisher.test(table.2mean)

#independence test for pam
chisq.test(table.pam)
fisher.test(table.pam)

#(d)	Plot the cluster dendrograms for this part of golub data with single 
#linkage and complete linkage, using Euclidean distance.
#oncogene
par(mfrow=c(1,2)) 
plot(hclust(dist(golub[oncogene,],method="euclidian"),method="single"))
plot(hclust(dist(golub[oncogene,],method="euclidian"),method="complete"))

#antigen
par(mfrow=c(1,2)) 
plot(hclust(dist(golub[antigen,],method="euclidian"),method="single"))
plot(hclust(dist(golub[antigen,],method="euclidian"),method="complete"))


#problem 3
rm(list=ls())  #clears environment
#install.packages('ISLR')
library(ISLR)
ncidata<-NCI60$data
ncilabs<-NCI60$labs


#(a)	Using k-means clustering, produce a plot of K versus SSE, for 
#K=1,., 30. How many clusters appears to be there?
clusdata <- data.frame(ncidata)
cl.2mean <- kmeans(clusdata, centers=7, nstart = 10)
table(ncilabs, cl.2mean$cluster)

K<-(1:30); sse<-rep(NA,length(K)) 
for (k in K) { sse[k]<-kmeans(clusdata, centers=k,nstart = 10)$tot.withinss 
} 
plot(K, sse, type='o', xaxt='n'); axis(1, at = K, las=2)

#(b) Do K-medoids clustering (K=7) with 1-correlation as the 
#dissimilarity measure on the data.
table(ncilabs, pam(as.dist(1-cor(t(ncidata))), k=7)$cluster)

