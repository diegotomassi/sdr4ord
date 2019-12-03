#############################################################################
require(homals)
source("cv.PCAhomals.R")
set.seed(160480)
oldw <- getOption("warn")
options(warn = -1)

# BUENOS AIRES
traindata_binGBA <- read.csv("../data/EPH/traindata_binGBA.txt", header=FALSE)
XX = as.matrix(traindata_binGBA[,-1])
Y = as.matrix(traindata_binGBA[,1])
outGBA = cv.PCAhomals(XX,Y,model='logit',nfold=10)

print('GBA done')

# HUMID PAMPAS
traindata_binPampeana <- read.csv("../data/EPH/traindata_binPampeana.txt", header=FALSE)
X = as.matrix(traindata_binPampeana[,-1])
Y = as.matrix(traindata_binPampeana[,1])
outPampeana = cv.PCAhomals(X,Y,model='logit',nfold=10)

print('PAMPA done')

# NORTHwEST
traindata_binNOA <- read.csv("../data/EPH/traindata_binNOA.txt", header=FALSE)
X = as.matrix(traindata_binNOA[,-1])
Y = as.matrix(traindata_binNOA[,1])
outNOA = cv.PCAhomals(X,Y,model='logit',nfold=10)

print('NOA done')

# NORTHEAST
traindata_binNEA <- read.csv("../data/EPH/traindata_binNEA.txt", header=FALSE)
X = as.matrix(traindata_binNEA[,-1])
Y = as.matrix(traindata_binNEA[,1])
outNEA = cv.PCAhomals(X,Y,model='logit',nfold=10)

print('NEA done')


# PATAGONIA
traindata_binPatagonia <- read.csv("../data/EPH/traindata_binPatagonia.txt", header=FALSE)
X = as.matrix(traindata_binPatagonia[,-1])
Y = as.matrix(traindata_binPatagonia[,1])
outPatagonia = cv.PCAhomals(X,Y,model='logit',nfold=10)

print('PATAGONIA done')

output = c(outGBA$mse,outPampeana$mse,outNOA$mse,outNEA$mse,outPatagonia$mse)
print(output)

output = c(outGBA$std,outPampeana$std,outNOA$std,outNEA$std,outPatagonia$std)
print(output)


options(warn = oldw)
