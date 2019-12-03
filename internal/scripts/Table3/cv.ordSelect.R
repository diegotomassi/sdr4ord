cv.ordSelect <- function(x,y,lambda,model="linear",nfold=10){
  nlambda = length(lambda)
  lambda = sort(lambda,decreasing=TRUE)
  N = length(y)
  x = cbind(rep(1,N),x)
  
  # start CV procedure
  foldsize = floor(N/nfold)
  msecv = matrix(0,nfold,nlambda)
  print(lambda)
  for (fold in 1:nfold){
    if (fold < nfold){idx = seq((fold-1)*foldsize+1,fold*foldsize)}
    else {idx = seq((fold-1)*foldsize+1,fold*foldsize)}
    xtr = x[-idx,]; ytr = y[-idx];xts=x[idx,]; yts = y[idx]
    m = ordSelect(xtr,ytr,lambda=lambda,model=model)
    if (model=='logit'){
      yhat = (predict(m,xts)>.5)
      for (k in 1:nlambda){
        msecv[fold,k] = mean(yhat[,k]!=yts)
      }
    }else{
      yhat = predict(m,xts)
      for (k in 1:nlambda){
        msecv[fold,k] = mean((yhat[,k]-yts)^2)
      }
    }
  }
  mseav = apply(msecv,2,mean)
  varav = apply(msecv,2,sd)
  kopt=which.min(mseav)
  print(mseav[c(kopt,nlambda)])
  print(apply(msecv,2,sd)[c(kopt,nlambda)])
  
  lambdaopt = lambda[kopt]
  mopt = ordSelect(x,y,lambda=lambdaopt,model=model)
  return(list(coefficients = mopt$coefficients,lambda.opt = lambdaopt,mse = mseav[kopt],std=varav[kopt]))
}