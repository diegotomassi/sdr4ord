# prediction error for INCOME PER CAPITA using PCApoly
# Compare obtained results with entries for LASSOord in Table 3.
##########################################################################
require(polycor)
source("cv.PCApoly.R")

print("Running example for SES Index construction using PCA with polychoric correlations")
print("This can take some time...")

set.seed(160480)
oldw <- getOption("warn")
options(warn = -1)

# BUENOS AIRES
traindata_GBA <- read.csv("../data/EPH/traindata_GBA.txt", header=FALSE)
X = as.matrix(traindata_GBA[,-1])
Y = as.matrix(traindata_GBA[,1]/1000)
outGBA = cv.PCApoly(X,Y,model='linear',nfold=10)

# HUMID PAMPAS
traindata_Pampeana <- read.csv("../data/EPH/traindata_Pampeana.txt", header=FALSE)
X = as.matrix(traindata_Pampeana[,-1])
Y = as.matrix(traindata_Pampeana[,1]/1000)
outPampeana = cv.PCApoly(X,Y,model='linear',nfold=10)

# NORTHwEST
traindata_NOA <- read.csv("../data/EPH/traindata_NOA.txt", header=FALSE)
X = as.matrix(traindata_NOA[,-1])
Y = as.matrix(traindata_NOA[,1]/1000)
outNOA = cv.PCApoly(X,Y,model='linear',nfold=10)

# NORTHEAST
traindata_NEA <- read.csv("../data/EPH/traindata_NEA.txt", header=FALSE)
X = as.matrix(traindata_NEA[,-1])
Y = as.matrix(traindata_NEA[,1]/1000)
outNEA = cv.PCApoly(X,Y,model='linear',nfold=10)


# PATAGONIA
traindata_Patagonia <- read.csv("../data/EPH/traindata_Patagonia.txt", header=FALSE)
X = as.matrix(traindata_Patagonia[,-1])
Y = as.matrix(traindata_Patagonia[,1]/1000)
outPatagonia = cv.PCApoly(X,Y,model='linear',nfold=10)

output = c(outGBA$mse,outPampeana$mse,outNOA$mse,outNEA$mse,outPatagonia$mse)
print(output)

output = c(outGBA$std,outPampeana$std,outNOA$std,outNEA$std,outPatagonia$std)
print(output)



options(warn = oldw)




# 
# 
# ##
# 
# redDATA <- read.csv("C:/Users/diegoT/Dropbox/_proyectos/ordinal/private/codeMATLABv4/redDATA.txt", header=FALSE)
# ipcf = as.matrix(redDATA[,1])
# SES_index = -as.matrix(redDATA[,2])
# SES_index = (SES_index-min(SES_index))/(max(SES_index)-min(SES_index))
# plot(SES_index,ipcf)
# library(MASS)
# library(car)
# m = lm(ipcf~SES_index)
# mmp(m,SES_index,sd=TRUE)
# boxcox((ipcf+1)~SES_index)
# boxcox((ipcf+1)^.33~SES_index)
# 
# ipcf = ipcf+1
# m1 = lm(ipcf^.33~SES_index+I(SES_index^2))
# mmp(m1,SES_index,sd=TRUE)
# 
# 
# 
# 
# evec = polycorPCA(X,1)
# SES_pca = X%*%evec
# mpca = lm(Y~SES_pca)
# mmp(mpca,SES_pca,sd=TRUE)
# SES_pca = -X%*%evec
# SES_pca = (SES_pca-min(SES_pca))/(max(SES_pca)-min(SES_pca))
# m1 = lm(ipcf^.33~SES_pca+I(SES_pca^2))
# mmp(m1,SES_pca,sd=TRUE)
