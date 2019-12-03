cv.PCApoly <- function(X,Y,model="linear",nfold=10){
  N = length(Y)
  
  # start CV procedure
  foldsize = floor(N/nfold)
  msecv = matrix(0,nfold,1)
  
  for (fold in 1:nfold){
    # data partition for k-fold
    if (fold < nfold){idx = seq((fold-1)*foldsize+1,fold*foldsize)}
    else {idx = seq((fold-1)*foldsize+1,fold*foldsize)}
 
    evec = polycorPCA(X[-idx,],1)

    pcadata = data.frame(Y=Y[-idx],X=X[-idx,]%*%evec)
    newdatapca = data.frame(X=X[idx,]%*%evec)
    
    if (model=='logit'){ # logistic regression for binary response
      m=glm(Y~X,data=pcadata,family='binomial')
      yhat = (predict(m,newdata=newdatapca)>.5)
      msecv[fold] = mean(yhat!=Y[idx])
#
    }else{# linear model for continuous response
      m.pca = lm(Y[-idx]~(X[-idx,]%*%evec))
      yhat = cbind(rep(1,length(idx)),X[idx,]%*%evec)%*%coef(m.pca)
      msecv[fold] = sum((yhat-Y[idx])^2)/length(idx);
    }
  }
  mseav = mean(msecv)
  return(list(mse = mseav,std=sd(msecv)))
}


## Auxiliar functions
require(polycor)


# computes polychoric correlations
polycorMtx <- function(x){
  p = ncol(x)
  out = matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:i){
      val = polychor(x[,i],x[,j])
      out[i,j] = val
      out[j,i] = val
    }
  }
  return(out)
}


polycorPCA <- function(x,dim){
  m = polycorMtx(x)
  valvecs = eigen(m)
  return(valvecs$vectors[,1:dim])
}

