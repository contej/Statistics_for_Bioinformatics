rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio

#Problem 1
#(a) Preprocess the raw data set into an expression data set using: 
#"mas" background correction method, the "quantiles" normalization method, 
#"pmonly" pm correction method and "medianpolish" summary method. 

library(ArrayExpress);library(affy)
#getAE('E-MEXP-1551', path = 'C:/yeast',  type = "full")

yeast.raw <-  ReadAffy(celfile.path= 'C:/yeast' )

yeast.raw
eset <- expresso(yeast.raw,bgcorrect.method="mas",
                 normalize.method="quantiles",pmcorrect.method="pmonly",
                 summary.method="medianpolish")
summary (eset)


#(b) Print out the mean expression values for the first five genes across all samples. 

firstFive<-eset[,eset$sample<=5] #Take first 5 samples
firstFive.copy<-firstFive #Make a copy data set
meanFive <- apply(exprs(firstFive), 2, mean) #mean for each of the 5 samples
meanFive

#(c) How many genes and how many samples are in the preprocessed expression data set?
yeast.raw
eset


#Problem 2
#(a) What is the annotation package for the yeast data set in question 1? 
#Install the annotation package from Bioconductor. 
annotation(yeast.raw)

#source("https://bioconductor.org/biocLite.R")
#biocLite("yeast2.db")

#(b) Search the 1769308_at gene GO numbers related to Molecular Function (MF). 
#How many GO numbers do you get?
library("GO.db")
library("yeast2.db")
library(annotate)

go1769308 <- get("1769308_at", env = yeast2GO)
gonr <- getOntology(go1769308,"MF")
gonr

#(c) Find the GO parents of the GO IDs in part (b). 
#How many GO parents are there? 

gP <- getGOParents(gonr) #parents for GO numbers 
pa <- sapply(gP,function(x) x$Parents) #extract GO numbers for parents
gonrpa <- unique(c(unlist(pa))) # GO numbers of parents)
gonrpa

#(d) Find the GO children of the GO IDs in part (b). 
#How many GO children are there?
gC <- getGOChildren(gonr) #children for GO numbers
ch <- sapply(gC,function(x) x$Children) #extract GO numbers for children
gonrpa <- unique(c(unlist(ch))) # GO numbers of children
gonrpa

#Problem 3
rm(list=ls())  #clears environment
library("genefilter")
data(ALL, package = "ALL")
ALLB <- ALL[,ALL$BT %in% c("B2","B3")]
f1 <- function(x) (wilcox.test(x, exact=F)$p.value < 0.001)
f2 <- function(x) (t.test(x)$p.value < 0.001)

sel2w <- genefilter(exprs(ALL[,ALLB$BT=="B2"]), filterfun(f1))
sel3w <- genefilter(exprs(ALL[,ALLB$BT=="B3"]), filterfun(f1))

sel2t <- genefilter(exprs(ALL[,ALLB$BT=="B2"]), filterfun(f2))
sel3t <- genefilter(exprs(ALL[,ALLB$BT=="B3"]), filterfun(f2))


library(limma)
x <- matrix(as.integer(c(sel2w,sel3w)),ncol = 2,byrow=FALSE)
colnames(x) <- c("sel2w","sel3w")
vcw <- vennCounts(x, include="both")
vennDiagram(vcw)

y <- matrix(as.integer(c(sel2t,sel3t)),ncol = 2,byrow=FALSE)
colnames(y) <- c("sel2t","sel3t")
vct <- vennCounts(y, include="both")
vennDiagram(vct)


#(d) What is the annotation package for the ALL data set? 
#Find the GO numbers for "oncogene".
annotation(ALL)

library("ALL"); data(ALL)
library("GO.db")
library("hgu95av2.db")
library(annotate)

library(help=hgu95av2.db)

ALLBt <- ALL[,ALL$BT %in% c("B","B1","B2","B3","B4")]
bnames <- featureNames(ALLBt)
goncogene <- mget(bnames, env = hgu95av2GENENAME)
length(grep("oncogene",goncogene))



#(e)
affynames <- featureNames(ALLB)
genenames <- mget(affynames, env = hgu95av2GENENAME)
ax<-grep("oncogene",genenames)

affytot <- unique(featureNames(ALLB))
genenamestot <- mget(affytot, env = hgu95av2GENENAME)
length(grep("oncogene",genenamestot))

#Problem 4
#(a) Select the persons with B-cell leukemia which are in 
#stage B1, B2, and B3. 
library("ALL")
library("limma");
allB <- ALL[,which(ALL$BT %in% c("B1","B2","B3"))]
facB123 <- factor(allB$BT)
design.ma <- model.matrix(~ 0 + facB123)
colnames(design.ma) <- c("B1","B2","B3")
fit <- lmFit(allB, design.ma)
fit <- eBayes(fit)
topTable(fit1, coef=2,5,adjust.method="fdr")


#(c) Select the persons with B-cell leukemia which are in 
#stage B1, B2, and B3. 
library("ALL")
library("limma");
allB <- ALL[,which(ALL$BT %in% c("B1","B2","B3"))]
facB123 <- factor(allB$BT)
cont.ma <- makeContrasts(B1-B2,B2-B3, levels=facB123)
design.ma <- model.matrix(~ 0 + facB123)
colnames(design.ma) <- c("B1","B2","B3")
fit <- lmFit(allB, design.ma)
fit1 <- contrasts.fit(fit, cont.ma)
fit1 <- eBayes(fit1)

dim( topTable(fit1, number=Inf, p.value=0.01,adjust.method="fdr"))

topTable(fit1, number=5,p.value=0.01,adjust.method="fdr")