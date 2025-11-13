Y1 <- matrix(0, n1, p)
Y2 <- matrix(0, n2, p)

Z1 <- matrix(0, n1, p)
Z2 <- matrix(0, n2, p)

###### mean
md <- c(runif(floor(p*s), -d/2,d/2), rep(0, p - floor(p*s)))
samp <- sample(1:p,p)

m2 <- numeric(p)
m1 <- rep(0,p)
m2[samp] <- md

pd <- p/4
###### variance
if(COV == 1 | COV == 3){
  r1 <- 0.999
  r2 <- 0.9
  r3 <- 0.7
  r4 <- 0.5
  
  A1 <- (1-r1)*diag(pd) + matrix(r1, pd, pd)
  A2 <- (1-r2)*diag(pd) + matrix(r2, pd, pd)
  A3 <- (1-r3)*diag(pd) + matrix(r3, pd, pd)
  A4 <- (1-r4)*diag(pd) + matrix(r4, pd, pd)
  
  D2 <- matrix(0,pd,pd)
  D1 <- matrix(1/pd,pd,pd)
  D3 <- matrix(-1/pd,pd,pd)
  
  C1 <- cbind(A1,D1,D2,D3)
  C2 <- cbind(D1,A2,D1,D2)
  C3 <- cbind(D2,D1,A3,D1)
  C4 <- cbind(D3,D2,D1,A4)
  
  EC1 <- rbind(C1,C2,C3,C4)
  
  EC1evec <- eigen(EC1)$vector
  EC1eval <- diag(eigen(EC1)$values)
  EC1eval.half <- sqrt(EC1eval)
  L1 <- EC1evec %*% EC1eval.half
  
  if(COV == 1){
    L2 <- L1
  }else{
    C1 <- cbind(A4,D1,D2,D3)
    C2 <- cbind(D1,A3,D1,D2)
    C3 <- cbind(D2,D1,A2,D1)
    C4 <- cbind(D3,D2,D1,A1)
    
    EC2 <- rbind(C1,C2,C3,C4)
    
    EC2evec <- eigen(EC2)$vector
    EC2eval <- diag(eigen(EC2)$values)
    EC2eval.half <- sqrt(EC2eval)
    L2 <- EC2evec %*% EC2eval.half
  }
}

if(COV == 2){
  r1 <- 0.5
  r2 <- 0.3
  r3 <- 0.1
  
  A1 <- toeplitz(r1^(0:(pd-1)))
  A2 <- toeplitz(r2^(0:(pd-1)))
  A3 <- toeplitz(r3^(0:(pd-1)))
  A4 <- diag(pd)
  
  D2 <- matrix(0,pd,pd)
  D1 <- matrix(1/pd,pd,pd)
  D3 <- matrix(-1/pd,pd,pd)
  
  C1 <- cbind(A1,D1,D2,D3)
  C2 <- cbind(D1,A2,D1,D2)
  C3 <- cbind(D2,D1,A3,D1)
  C4 <- cbind(D3,D2,D1,A4)
  
  EC1 <- rbind(C1,C2,C3,C4)
  
  EC1evec <- eigen(EC1)$vector
  EC1eval <- diag(eigen(EC1)$values)
  EC1eval.half <- sqrt(EC1eval)
  L1 <- EC1evec %*% EC1eval.half
  
  
  C1 <- cbind(A4,D1,D2,D3)
  C2 <- cbind(D1,A1,D1,D2)
  C3 <- cbind(D2,D1,A2,D1)
  C4 <- cbind(D3,D2,D1,A3)
  
  EC2 <- rbind(C1,C2,C3,C4)
  
  EC2evec <- eigen(EC2)$vector
  EC2eval <- diag(eigen(EC2)$values)
  EC2eval.half <- sqrt(EC2eval)
  L2 <- EC2evec %*% EC2eval.half
}


G <- diag(p) - matrix(1, p, p)/p


##### Data

if(SETUP == 1){
  for(i in 1:n1){
    Z1[i,] <- (m1 + L1 %*% rnorm(p))
    Y1[i,] <- G %*% Z1[i,]
  }
  for(i in 1:n2){
    Z2[i,] <- (m2 + L2 %*% rnorm(p))
    Y2[i,] <- G %*% Z2[i,]
  }
}

if(SETUP == 2){
  for(i in 1:n1){
    Z1[i,] <- (m1 + L1 %*% (rt(p, 5)/sqrt(5/3)))
    Y1[i,] <- G %*% Z1[i,]
  }
  for(i in 1:n2){
    Z2[i,] <- (m2 + L2 %*% (rt(p, 5)/sqrt(5/3)))
    Y2[i,] <- G %*% Z2[i,]
  }
}
