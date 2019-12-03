# prediction error for POVERTY using PCApoly
# Compare obtained results with entries for LASSOord in Table 3.
##########################################################################

library(polycor)
source("./cv.PCApoly.R")

set.seed(160480)
oldw <- getOption("warn")
options(warn = -1)

# BUENOS AIRES
traindata_binGBA <- read.csv("../data/EPH/traindata_binGBA.txt", header=FALSE)
X = as.matrix(traindata_binGBA[,-1])
Y = as.matrix(traindata_binGBA[,1])
outGBA = cv.PCApoly(X,Y,model='logit',nfold=10)

# HUMID PAMPAS
traindata_binPampeana <- read.csv("../data/EPH/traindata_binPampeana.txt", header=FALSE)
X = as.matrix(traindata_binPampeana[,-1])
Y = as.matrix(traindata_binPampeana[,1])
outPampeana = cv.PCApoly(X,Y,model='logit',nfold=10)

# NORTHwEST
traindata_binNOA <- read.csv("../data/EPH/traindata_binNOA.txt", header=FALSE)
X = as.matrix(traindata_binNOA[,-1])
Y = as.matrix(traindata_binNOA[,1])
outNOA = cv.PCApoly(X,Y,model='logit',nfold=10)

# NORTHEAST
traindata_binNEA <- read.csv("../data/EPH/traindata_binNEA.txt", header=FALSE)
X = as.matrix(traindata_binNEA[,-1])
Y = as.matrix(traindata_binNEA[,1])
outNEA = cv.PCApoly(X,Y,model='logit',nfold=10)


# PATAGONIA
traindata_binPatagonia <- read.csv("../data/EPH/traindata_binPatagonia.txt", header=FALSE)
X = as.matrix(traindata_binPatagonia[,-1])
Y = as.matrix(traindata_binPatagonia[,1])
outPatagonia = cv.PCApoly(X,Y,model='logit',nfold=10)

output = c(outGBA$mse,outPampeana$mse,outNOA$mse,outNEA$mse,outPatagonia$mse)
print(output)

output = c(outGBA$std,outPampeana$std,outNOA$std,outNEA$std,outPatagonia$std)
print(output)


options(warn = oldw)
