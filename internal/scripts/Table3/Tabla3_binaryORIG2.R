# prediction error for POVERTY using full set of predictors with dummy variables.
# This strategy is named FULL-I in the mansucript.
# Compare obtained results with entries for LASSOord in Table 4
##########################################################################
require(ordPens)
lambda = 0

# BUENOS AIRES
traindata_binGBA <- read.csv("../data/EPH/traindata_binGBA.txt", header=FALSE)
X = traindata_binGBA[,-1]
Y = traindata_binGBA[,1]
outbinGBA = cv.ordSelect(X,Y,lambda,model='logit',nfold=10)

# HUMID PAMPAS
traindata_binPampeana <- read.csv("../data/EPH/traindata_binPampeana.txt", header=FALSE)
X = traindata_binPampeana[,-1]
Y = traindata_binPampeana[,1]
lambda = 0
outbinPampeana = cv.ordSelect(X,Y,lambda,model='logit',nfold=10)

# NORTHwEST
traindata_binNOA <- read.csv("../data/EPH/traindata_binNOA.txt", header=FALSE)
X = traindata_binNOA[,-1]
Y = traindata_binNOA[,1]
lambda = 0
outbinNOA = cv.ordSelect(X,Y,lambda,model='logit',nfold=10)

# NORTHEAST
traindata_binNEA <- read.csv("../data/EPH/traindata_binNEA.txt", header=FALSE)
X = traindata_binNEA[,-1]
Y = traindata_binNEA[,1]
lambda = 0
outbinNEA = cv.ordSelect(X,Y,lambda,model='logit',nfold=10)


# PATAGONIA
traindata_binPatagonia <- read.csv("../data/EPH/traindata_binPatagonia.txt", header=FALSE)
X = traindata_binPatagonia[,-1]
Y = traindata_binPatagonia[,1]
lambda = 0
outbinPatagonia = cv.ordSelect(X,Y,lambda,model='logit',nfold=10)

output = c(outbinGBA$mse,outbinPampeana$mse,outbinNOA$mse,outbinNEA$mse,outbinPatagonia$mse)
print(output)

output = c(outbinGBA$std,outbinPampeana$std,outbinNOA$std,outbinNEA$std,outbinPatagonia$std)
print(output)