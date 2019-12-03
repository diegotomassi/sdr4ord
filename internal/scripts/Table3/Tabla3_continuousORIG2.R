# prediction error for INCOME PER CAPITA using full set of predictors and dummy variables.
# This strategy is named FULL-I in the mansucript.
# Compare obtained results with entries for LASSOord in Table 3.
##########################################################################
require(ordPens)
source("cv.ordSelect.R")

set.seed(160480)
# BUENOS AIRES
traindata_GBA <- read.csv("../data/EPH/traindata_GBA.txt", header=FALSE)
X = traindata_GBA[,-1]
Y = traindata_GBA[,1]/1000
lambda = 0
outGBA = cv.ordSelect(X,Y,lambda,model='linear',nfold=10)

# HUMID PAMPAS
traindata_Pampeana <- read.csv("../data/EPH/traindata_Pampeana.txt", header=FALSE)
X = traindata_Pampeana[,-1]
Y = traindata_Pampeana[,1]/1000
lambda = 0
outPampeana = cv.ordSelect(X,Y,lambda,model='linear',nfold=10)

# NORTHwEST
traindata_NOA <- read.csv("../data/EPH/traindata_NOA.txt", header=FALSE)
X = traindata_NOA[,-1]
Y = traindata_NOA[,1]/1000
lambda = 0
outNOA = cv.ordSelect(X,Y,lambda,model='linear',nfold=10)

# NORTHEAST
traindata_NEA <- read.csv("../data/EPH/traindata_NEA.txt", header=FALSE)
X = traindata_NEA[,-1]
Y = traindata_NEA[,1]/1000
lambda = 0
outNEA = cv.ordSelect(X,Y,lambda,model='linear',nfold=10)


# PATAGONIA
traindata_Patagonia <- read.csv("../data/EPH/traindata_Patagonia.txt", header=FALSE)
X = traindata_Patagonia[,-1]
Y = traindata_Patagonia[,1]/1000
lambda = 0
outPatagonia = cv.ordSelect(X,Y,lambda,model='linear',nfold=10)

output = c(outGBA$mse,outPampeana$mse,outNOA$mse,outNEA$mse,outPatagonia$mse)
print(output)

output = c(outGBA$std,outPampeana$std,outNOA$std,outNEA$std,outPatagonia$std)
print(output)