#Class 4
#Analysis of Microarray data

source("https://bioconductor.org/biocLite.R")
biocLite("SpikeInSubset") # biocLite is a function

#loading this package
library(SpikeInSubset) # library where we have datasets; SpikeInSubset is package
library(affy)

# Load data
data(spikein95) # spikein95 is dataset

#Chek the chips
image(spikein95)

#collect the gene ids from dataset and put it into an object
ids <- geneNames(spikein95) # we make an object, geneNames is fuction

#ids - show the ids
ids[1:10]


# perform Mas 5.0 normalization
Mas  <- mas5(spikein95)
#Mas <- mas5(spikein95, col=2:5) # put normalized data in object, mas5 is a function
mas_log <- log2(exprs(Mas)) # taking expression values with the function exprs

# Box plot for row data and boxplot for normalized data on log fold
boxplot(spikein95)
x11() # separate windows
boxplot(mas_log)

#density plot
density <- density(mas_log[,1]) #density function
plot(density, main="Mas mormalization")

density2 <- density(mas_log[,2]) # column 2 second sample
lines(density2, col="red")

density3 <- density(mas_log[,3])
lines(density3, col="blue")

# Making MA plots
# M: difference in average log intensities
# A: average log intensities

d <- rowMeans(mas_log[,1:3] - rowMeans(mas_log[,4:6]))
a <- rowMeans(mas_log)

# Plotting data
plot(a,d, ylim = c(-5,5), main= "Mas 5.0 normalized MA plot",
xlab = "A", ylab = "M", pch = ".") # 2 objects a nad b

abline(h = c(-1.5,1.5), col="red") 

rma.eset <- rma(spikein95)
rma.e <- exprs(rma.eset) # we collect expression values 

head(rma.e)
head(mas_log)

x11()
boxplot(mas_log, col = 2:5, main = "Mas 5 Norm")
x11()
boxplot(rma.e, col = 2:5, main = "RMA Norm")

# finding specific genes/probes in MA plot

spikedn <- colnames(pData(spikein95))
spikedIndex <- match(spikedn, featureNames(Mas))
points(a[spikedIndex], d[spikedIndex], pch=19, col="red")

# Volcano Plots
source("https://bioconductor.org/biocLite.R") #update - latest version of Bioconductor
biocLite("genefilter")  #??? Install specific packages
library("genefilter") 
pData(rma.eset) <- pData(Mas)
tt <- rowMeans(rma.e)
lod <- tt
plot(d, lod, cex = 0.25, main = "Volcano plot for MA", xlim = c(-2, 2), xlab = "M", ylab = "A", yaxt = "n")
axis(2, at = seq(0, 3, by = 1), labels = 10^(-seq(0, 3, by = 1)))
points(d[spikedIndex], lod[spikedIndex], pch = 19, col = "red")
abline(h = 2, v = c(-1, 1))


# Exercise to practice

https://rawgit.com/bioinformatics-core-shared-training/microarray-analysis/master/de-tutorial.html

https://rawgit.com/bioinformatics-core-shared-training/microarray-analysis/master/affymetrix.nb.html



library(affy)
setwd("C:/Users/Natalia/Documents/GitHub/Class_4/estrogen/estrogen")
targetsFile <- "estrogen/estrogen.txt"
targetsFile
pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)
pData(pd)

ER <- pData(pd)$estrogen
Time <- factor(pData(pd)$time.h)
design <- model.matrix(~ER+Time)
design

design2 <- model.matrix(~ER*Time)
design2

#raw <-ReadAffy(celfile.path = "C:/Users/Natalia/Documents/GitHub/Class_4/estrogen/estrogen", filenames=rownames(pData(pd)),phenoData = pd)
#raw
raw <-ReadAffy(celfile.path = "estrogen", filenames=rownames(pData(pd)),phenoData = pd) #!!!!
raw


boxplot(raw,col="red",las=2) #las=2 to make labels perpendicular to the axis

par(mfrow=c(2,1)) ## 2 gistograms in 1 windov

hist(log2(pm(raw[,1])),breaks=100,col="steelblue",main="PM",xlim=c(4,14))
hist(log2(mm(raw[,1])),breaks=100,col="steelblue",main="MM",xlim=c(4,14))



source("https://bioconductor.org/biocLite.R")
biocLite("affyPLM")  ##installation path not writeable, unable to update packages: foreign, survival???
library(affyPLM)

plmset <- fitPLM(raw) ##This function converts an AffyBatch into an PLMset by fitting a specified robust linear model to the probe level data. 
NUSE(plmset,las=2)

#bad <- ReadAffy(celfile.path = "C:/Users/Natalia/Documents/GitHub/Class_4/estrogen/estrogen",filenames="bad.cel")
bad <- ReadAffy(celfile.path = "estrogen/",filenames="bad.cel")
image(bad)


par(mfrow=c(2,4))
image(raw[,1])
image(raw[,2])

image(raw[,3])
image(raw[,4])

image(raw[,5])
image(raw[,6])

image(raw[,7])
image(raw[,8])
##If we are happy with the quality of the raw data we can proceed to the next step


## Normalization
eset <- rma(raw)

library(limma) ## ?????????? there is no package called ‘limma’
fit1 <- lmFit(eset, design)
fit1 <- eBayes(fit1)
topTable(fit1, coef=2)

fit2 <- lmFit(eset, design2)
fit2 <- eBayes(fit2)
topTable(fit2, coef=2)


#Annotation of samples
library(GEOquery)
library(limma)
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz"
filenm <- "GSE33126_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)
gse <- getGEO(filename=filenm)


