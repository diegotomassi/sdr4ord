# This script computes SES INDEX for differents regions in Argentina using 
# PCAover polychoric correlations (PCApoly).
# Compare obtained coefficients with entries for PCApoly in Table 5 and Table 6
###########################################################################
source("auxfuns.R")

oldw <- getOption("warn")
options(warn = -1)

traindata_GBA <- read.csv("../data/EPH/traindata_GBA.txt", header=FALSE)
m = polycorPCA(traindata_GBA[,-1],1)
PCApoly_GBA= m

traindata_Pampeana <- read.csv("../data/EPH/traindata_Pampeana.txt", header=FALSE)
m = polycorPCA(traindata_Pampeana[,-1],1)
PCApoly_Pampeana= m

traindata_NOA <- read.csv("../data/EPH/traindata_NOA.txt", header=FALSE)
m = polycorPCA(traindata_NOA[,-1],1)
PCApoly_NOA= m

traindata_NEA <- read.csv("../data/EPH/traindata_NEA.txt", header=FALSE)
m = polycorPCA(traindata_NEA[,-1],1)
PCApoly_NEA= m

traindata_Patagonia <- read.csv("../data/EPH/traindata_Patagonia.txt", header=FALSE)
m = polycorPCA(traindata_Patagonia[,-1],1)
PCApoly_Patagonia = m


loadings = cbind(PCApoly_GBA,PCApoly_Pampeana,PCApoly_NOA,PCApoly_NEA,PCApoly_Patagonia)
print(loadings)

options(warn = oldw)

