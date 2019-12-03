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


