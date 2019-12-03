cv.PCAhomals <- function(X,Y,model="linear",nfold=10){
  N = length(Y)
  
  # start CV procedure
  foldsize = floor(N/nfold)
  msecv = matrix(0,nfold,1)
  
  for (fold in 1:nfold){
    # data partition for k-fold
    if (fold < nfold){idx = seq((fold-1)*foldsize+1,fold*foldsize)}
    else {idx = seq((fold-1)*foldsize+1,fold*foldsize)}
 
    evec = as.matrix(homals(as.data.frame(X[-idx,]),ndim=1,level='ordinal')$loadings)

    pcadata = data.frame(Y=Y[-idx],X=X[-idx,]%*%as.numeric(evec))
    newdatapca = data.frame(X=X[idx,]%*%as.numeric(evec))
    
    if (model=='logit'){ # logistic regression for binary response
      m=glm(Y~X,data=pcadata,family='binomial')
      yhat = (predict(m,newdata=newdatapca)>.5)
      msecv[fold] = mean(yhat!=Y[idx])
#
    }else{# linear model for continuous response
      m.pca = lm(Y[-idx]~(X[-idx,]%*%as.numeric(evec)))
      yhat = as.matrix(cbind(rep(1,length(idx)),X[idx,]%*%as.numeric(evec)))%*%coef(m.pca)
      msecv[fold] = sum((yhat-Y[idx])^2)/length(idx);
    }
  }
  mseav = mean(msecv)
  return(list(mse = mseav,std=sd(msecv)))
}


