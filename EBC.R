EBC <- function (X, Y, N = 1e3){
  
  n<-nrow(X)
  m<-nrow(Y)
  delta2nm<-n^(1/2)*(n+m)^(-1/2)
  delta2mn<-m^(1/2)*(n+m)^(-1/2)
  
  meanX<-n^(1/2)*as.matrix(colMeans(X))
  meanY<-m^(1/2)*as.matrix(colMeans(Y))
  
  TestStatistic<-norm(delta2mn*meanX-delta2nm*meanY, type = "M")
  
  cX <- matrix(rep(colMeans(X),n),n,byrow = T)
  cY <- matrix(rep(colMeans(Y),m),m,byrow = T)
  
  SSS<-numeric(N)
  for(j in 1:N){
    s1 <- sample(1:n, n, replace = T)
    s2 <- sample(1:m, m, replace = T)
    X.cent <- X[s1,] - cX
    Y.cent <- Y[s2,] - cY
    SeX<-n^(1/2)*as.matrix(colMeans(X.cent)) 
    SeY<-m^(1/2)*as.matrix(colMeans(Y.cent))
    SSS[j]<-norm(delta2mn*SeX-delta2nm*SeY, type = "M")
  }
  
  pval<-mean(TestStatistic < SSS)
  
  return(pval)
}
