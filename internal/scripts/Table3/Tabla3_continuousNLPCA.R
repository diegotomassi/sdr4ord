#############################################################################
require(homals)
source("cv.PCAhomals.R")
print("Running example for SES Index construction using NLPCA as implemented in package HOMALS")
print("This can take some time...")

set.seed(160480)
oldw <- getOption("warn")
options(warn = -1)

# BUENOS AIRES
traindata_GBA <- read.csv("../data/EPH/traindata_GBA.txt", header=FALSE)
XX = as.matrix(traindata_GBA[,-1])
Y = as.matrix(traindata_GBA[,1])/1000
outGBA = cv.PCAhomals(XX,Y,model='linear',nfold=10)

print('GBA done')

# HUMID PAMPAS
traindata_Pampeana <- read.csv("../data/EPH/traindata_Pampeana.txt", header=FALSE)
X = as.matrix(traindata_Pampeana[,-1])
Y = as.matrix(traindata_Pampeana[,1])/1000
outPampeana = cv.PCAhomals(X,Y,model='linear',nfold=10)

print('PAMPA done')

# NORTHwEST
traindata_NOA <- read.csv("../data/EPH/traindata_NOA.txt", header=FALSE)
X = as.matrix(traindata_NOA[,-1])
Y = as.matrix(traindata_NOA[,1])/1000
outNOA = cv.PCAhomals(X,Y,model='linear',nfold=10)

print('NOA done')

# NORTHEAST
traindata_NEA <- read.csv("../data/EPH/traindata_NEA.txt", header=FALSE)
X = as.matrix(traindata_NEA[,-1])
Y = as.matrix(traindata_NEA[,1])/1000
outNEA = cv.PCAhomals(X,Y,model='linear',nfold=10)

print('NEA done')


# PATAGONIA
traindata_Patagonia <- read.csv("../data/EPH/traindata_Patagonia.txt", header=FALSE)
X = as.matrix(traindata_Patagonia[,-1])
Y = as.matrix(traindata_Patagonia[,1])/1000
outPatagonia = cv.PCAhomals(X,Y,model='linear',nfold=10)

print('PATAGONIA done')

output = c(outGBA$mse,outPampeana$mse,outNOA$mse,outNEA$mse,outPatagonia$mse)
print(output)

output = c(outGBA$std,outPampeana$std,outNOA$std,outNEA$std,outPatagonia$std)
print(output)


options(warn = oldw)
