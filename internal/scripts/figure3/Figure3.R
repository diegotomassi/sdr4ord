# This script reproduces Figure 3 in the manuscript
# "Sufficient dimension reduction for ordinal predictors",
# by L. Forzani, R. Garcia arancibia, P. Llop and D. Tomassi.
#
###############################################################
library(MASS)
library(car)
source("auxfuns.R")

## SES index using PFCord
redDATA <- read.csv("./proj_PFCord.txt", header=FALSE)
ipcf = as.matrix(redDATA[,1])
SES_index = as.matrix(redDATA[,2])
SES_index = (SES_index-min(SES_index))/(max(SES_index)-min(SES_index))
plot(SES_index,ipcf)
m = lm(ipcf~SES_index)
mmp(m,SES_index,sd=TRUE)
library(MASS)
boxcox(m)
boxcox((ipcf+.00001)~SES_index)
boxcox((ipcf+.00001)^.33~SES_index)
ipcf =  ipcf+.00001
m1 = lm(ipcf^.33~SES_index+I(SES_index^2))
mmp(m1,SES_index,sd=TRUE,xlab="SES_PFCord")


## SES index derived from PCA
traindata_full_new <- read.csv("../data/EPH/traindata_full_new.txt", header=FALSE)
X = as.matrix(traindata_full_new[,-1])
#ipcf = as.matrix(traindata_full_new[,1])/1000
evec = polycorPCA(X,1)

SES_pca = -X%*%evec
SES_pca = (SES_pca-min(SES_pca))/(max(SES_pca)-min(SES_pca))
mpca = lm(ipcf~SES_pca)
mmp(mpca,SES_pca,sd=TRUE)
boxcox(mpca)
m1 = lm(ipcf^.33~SES_pca+I(SES_pca^2))
boxcox(m1)
mmp(m1,SES_pca,sd=TRUE,xlab="SES_PCApoly")

