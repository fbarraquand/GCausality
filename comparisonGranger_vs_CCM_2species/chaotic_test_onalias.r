waldtest.lmbis <- function(object, ..., test = c("F", "Chisq"))
{
  waldtest.default(object, ..., test = match.arg(test))
} 



tmax=800
  X=Y=rep(1,tmax)## Problem if I start with 1
  X[1]=runif(1,0.1,0.7)
  Y[1]=runif(1,0.1,0.7)
  for (t in 1:tmax)
  {
 
  X[t+1]=X[t]*(3.8-3.8*X[t]-0.0*Y[t])
  Y[t+1]=Y[t]*(3.5-3.5*Y[t]-0.0*X[t])
  }
  x=log(X[501:800])
  y=log(Y[501:800])
  # centering
  x=x-mean(x)
  y=y-mean(y)



order=4
lag_X_tmp=lapply(1:order, function(k) lag(x, -k))
lag_Y_tmp=lapply(1:order, function(k) lag(y, -k))
lag_X=matrix(unlist(lag_X_tmp),ncol=order,byrow=T)
lag_Y=matrix(unlist(lag_Y_tmp),ncol=order,byrow=T)
all <- cbind(x, y, lag_X, lag_Y)
  colnames(all) <- c("x", "y", paste("x", 1:order, sep = "_"), paste("y", 1:order, sep = "_"))
  ybis <- as.vector(all[,2])
  lagX <- as.matrix(all[,(1:order + 2)])
  lagY <- as.matrix(all[,(1:order + 2 + order)])

  fm <- lm(ybis ~ lagY + lagX)

waldtest(fm,2)
