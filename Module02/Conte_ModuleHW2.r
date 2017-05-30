rm(list=ls())  #clears environment
cat("\014") #clears the console in RStudio
library(multtest) # Load the 'multtest' package
data(golub) # Load the 'golub' data set

#Joshua Conte HW2
#Problem 1

#(a) Compute the mean expression values for every gene among "ALL" patients
#The data consist of gene expression values of 3051 genes (rows) from 38 
#leukemia patients (columns). The first Twenty seven patients are diagnosed as acute lymphoblastic 
#leukemia (ALL) and the remaining eleven as acute myeloid leukemia (AML).

#this line of code separates ALL and AML
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL", "AML"))

#The genes are in rows, so this will get us the mean expression for every genes for ALL patients
meanALL_row <- apply(golub[,gol.fac=="ALL"], 1, mean)
meanALL_row

#(b) Compute the mean expression values for every gene among "AML" patients.
#continuing from the information from problem 1 (a), we can take the information
#from meanALL_row and change golub[,gol.fac from ALL to AML to get mean for every AML gene.
meanAML_row <- apply(golub[,gol.fac=="AML"], 1, mean)
meanAML_row

#(c) Give the biological names of the three genes with the largest mean expression
#value among "ALL" patients.
#continuing from the information from problem 1 (a), we can take the information
#and arrange the order by row means in decreasing order
ALL_order <- order(meanALL_row,decreasing=TRUE) 

# Display the three genes with the largest mean expression
ALL_order[1:3]

#This gives us the gene numbers, which we can turn into a name by using the command
#golub.gnames[x,] where x is the gene number.
#the top 3 genes are 2784, 2478, and 748.

#Gene name for 2784
golub.gnames[2784,]

#Gene name for 2478
golub.gnames[2478,]

#Gene name for 748
golub.gnames[748,]

#(d) Give the biological names of the three genes with the largest mean expression 
#value among "AML" patients.
#We can take the same approach from 1 (c), we just need to change "ALL" to "AML.

# order by row means in decreasing order
AML_order <- order(meanAML_row,decreasing=TRUE) 

# Display the order of row means
AML_order[1:3]

#This gives us the gene numbers, which we can turn into a name by using the command
#golub.gnames[x,] where x is the gene number.
#the top 3 genes are 2586, 2440, and 2468.

#Gene name for 2586
golub.gnames[2586,]

#Gene name for 2440
golub.gnames[2440,]

#Gene name for 2468
golub.gnames[2468,]


#Problem 2
#(a) Save the expression values of the first five genes (in the first five rows) for the
#AML patients in a csv file "AML5.csv".

#this collects the first five genes for the AML patients and puts it into AML5
AML5 <- golub[1:5,gol.fac=="AML"]

#this writes the contents of AML5 into a CVS file.
write.csv(AML5,file="AML5.csv")

#(b) Save the expression values of the first five genes for the ALL patients in a plain
#text file "ALL5.txt".

#this collects the first five genes for the ALL patients and puts it into ALL5
ALL5 <- golub[1:5,gol.fac=="ALL"]

#this writes the contents of ALL5 into a CVS file.
write.csv(ALL5,file="ALL5.csv")

#(c) Compute the standard deviation of the expression values on the first patient,
#of the 100th to 200th genes (total 101 genes).

#this collects all the information for the 100th and 200th genes for 
#the first patient and puts it in SD1.
SD1 <- golub[100:200,1]

#this is the command to get the standard deviation of SD1
sd(SD1)

#(d) Compute the standard deviation of the expression values of every gene, across
#all patients. Find the number of genes with standard deviation greater than 1.

#this computes the standard deviation of all genes for each patient.
SD_every <- apply(golub[1:3051, 1:38],1,sd)

#this collects all standard deviations that are greater than 1
SDg1 <- SD_every[SD_every>1]

#the length lets us know how many elements are in SDg1, which is
#equal to the number of genes with standard deviation greater than 1
length(SDg1)

#(e) Do a scatter plot of the 101th gene expressions against the 102th gene
#expressions, label the x-axis and the y-axis with the genes' biological names using
#xlab= and ylab= control options.

#I will start by getting the names of the genes:
#Gene name for 101
golub.gnames[101,]

#Gene name for 102
golub.gnames[102,]

#101 = NUCLEAR PORE COMPLEX PROTEIN NUP214
#102 = PHOSPHATIDYLSERINE SYNTHASE I

#define x and y
x <- golub[101,]
y <- golub[102,]

#plot the information:
plot(x,y,xlab="NUCLEAR PORE COMPLEX PROTEIN NUP214",ylab="PHOSPHATIDYLSERINE SYNTHASE I", main="Scatter Plot")

#Problem 3.
#(a) Use exprs(ALL[,ALL$BT=="B1"] to extract the gene expressions from the patients in disease stage B1. Produce one histogram of these gene expressions in the this matrix.
#Since this is a new problem I clear the environment
rm(list=ls())  

#Load the library ALL
library(ALL); data(ALL)

#Extract B1 information
B1<-exprs(ALL[,ALL$BT=="B1"])

#draw histogram
hist(B1) 

#(b) Compute the mean gene expressions for every gene over these B1 patients.
meanB1 <- apply(exprs(ALL[,ALL$BT=="B1"]),1, mean)


#(c) Give the gene identifiers of the three genes with the largest mean.
#order by row means in decreasing order
o <- order(meanB1,decreasing=TRUE)

#the gene identifiers of the three genes with the largest mean.
meanB1[o[1:3]]

#Problem 4. (20 points)
#(a) Find the type of the trees data object.

# load trees dataset into workspace
data(trees)

#checking data object type
# is trees a vector?
is.vector(trees) 
# is trees a matrix?
is.matrix(trees)
# is trees a factor?
is.factor(trees) 
# is trees a array?
is.array(trees) 
# is trees a list?
is.list(trees) 
# is trees a data frame?
is.data.frame(trees) 
# are the elements in trees numbers?
is.numeric(trees) 
# are the elements in trees character strings?
is.character(trees) 

str(trees) # display info about trees. From below, we can see that trees has numbers in a two-dimens

# description (help) for the dataset:
?trees

#(b) Produce a figure with two overlaid scatterplots: 
#Height versus Girth, Volume versus Girth(The Girth is on the x-axis). 
#Do the Height plot with blue "+" symbols, and do the Volume plot with 
#red "o" symbols. You need to learn to set the ylim=control option so 
#that all points from the two plots can all show up on the merged figure.
#Hint: you should use plot() then points() to create the overlaid two scatterplots.	

#define Height, Girth, Volume
x<-trees$Girth
y<-trees$Height
z<-trees$Volume

#plot the data
plot(y~x,col="blue", pch=3, ylab="Height and Volume",xlab="Girth", ylim=c(min(z),max(y)), main="Height and Volume versus Girth of 31 Felled Black Cherry Trees")
points(x,z,col="red",pch=1 )
legend("bottomright",title="Legend", horiz=TRUE, c("Height","Volume"), pch=c(3,1),col=c("blue","red"),bg="grey96", inset=.02)