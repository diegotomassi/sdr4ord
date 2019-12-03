# This script computes SES INDEX for differents regions in Argentina using NLPCA
# Compare obtained coefficients with entries for NLPCA in Table 5 and Table 6
###########################################################################

require(homals)

traindata_GBA <- read.csv("../data/EPH/traindata_GBA.txt", header=FALSE)
m = homals(traindata_GBA[,-1],ndim=1,rank=1,level='ordinal',active=TRUE)
NLPCA_GBA= as.numeric(m$loadings)

traindata_Pampeana <- read.csv("../data/EPH/traindata_Pampeana.txt", header=FALSE)
m = homals(traindata_Pampeana[,-1],ndim=1,rank=1,level='ordinal',active=TRUE)
NLPCA_Pampeana= as.numeric(m$loadings)

traindata_NOA <- read.csv("../data/EPH/traindata_NOA.txt", header=FALSE)
m = homals(traindata_NOA[,-1],ndim=1,rank=1,level='ordinal',active=TRUE)
NLPCA_NOA= as.numeric(m$loadings)

traindata_NEA <- read.csv("../data/EPH/traindata_NEA.txt", header=FALSE)
m = homals(traindata_NEA[,-1],ndim=1,rank=1,level='ordinal',active=TRUE)
NLPCA_NEA= as.numeric(m$loadings)

traindata_Patagonia <- read.csv("../data/EPH/traindata_Patagonia.txt", header=FALSE)
m = homals(traindata_Patagonia[,-1],ndim=1,rank=1,level='ordinal',active=TRUE)
NLPCA_Patagonia= as.numeric(m$loadings)


loadings = cbind(NLPCA_GBA,NLPCA_Pampeana,NLPCA_NOA,NLPCA_NEA,NLPCA_Patagonia)
print(loadings)
