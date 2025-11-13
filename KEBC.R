KEBC <- function (X, N = 1e3){
  # X (list of matrices of real numbers) is a list of K (number of classes) elements, each contains a n_k x p matrix of observations
  # alpha (real number in (0,1)) is the level of significance
  # N (natural number) the number of bootstrap samples
  
  K <- length(X)
  Tn <- numeric(K*(K-1)/2)
  ind <- 1
  for(k in 1:(K-1)){
    for (l in (k+1):K){
      Tn[ind] <- KEBC_statistic(as.matrix(X[[k]]), as.matrix(X[[(l)]]))
      ind <- ind + 1
    }
  }
  TestStatistic <- max(Tn)
  
  cX <- list()
  for(k in 1:K){
    nk <- nrow(X[[k]])
    cX[[k]] <- matrix(rep(colMeans(X[[k]]),nk),nk,byrow = T)
  }
  
  SSS<-numeric(N)
  for(j in 1:N){
    SSSk <- numeric(K*(K-1)/2)
    ind <- 1
    for(k in 1:(K-1)){
      for (l in (k+1):K) {
        SSSk[ind] <- KEBC_bootstrap(as.matrix(X[[k]]), as.matrix(X[[l]]), as.matrix(cX[[k]]), as.matrix(cX[[l]]))
        ind <- ind + 1
      }
    }
    SSS[j] <- max(SSSk)
  }
  
  pval<-mean(TestStatistic < SSS)
  
  return(pval)
}

KEBC_statistic <- function(X, Y){
  
  n<-nrow(X)
  m<-nrow(Y)
  delta2nm<-n^(1/2)*(n+m)^(-1/2)
  delta2mn<-m^(1/2)*(n+m)^(-1/2)
  
  meanX<-n^(1/2)*as.matrix(colMeans(X))
  meanY<-m^(1/2)*as.matrix(colMeans(Y))
  
  TestStatistic <- as.numeric(norm(delta2mn*meanX-delta2nm*meanY, type = "M"))
  
  return(TestStatistic)
}

KEBC_bootstrap <- function(X, Y, cX, cY){
  
  n<-nrow(X)
  m<-nrow(Y)
  delta2nm<-n^(1/2)*(n+m)^(-1/2)
  delta2mn<-m^(1/2)*(n+m)^(-1/2)
  
  s1 <- sample(1:n, n, replace = T)
  s2 <- sample(1:m, m, replace = T)
  X.cent <- X[s1,] - cX
  Y.cent <- Y[s2,] - cY
  SeX<-n^(1/2)*as.matrix(colMeans(X.cent)) 
  SeY<-m^(1/2)*as.matrix(colMeans(Y.cent))
  SSS<-as.numeric(norm(delta2mn*SeX-delta2nm*SeY, type = "M"))
  
  return(SSS)
}
