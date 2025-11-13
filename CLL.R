CLL <- function(Y1, Y2, alpha = 0.05){
  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  n <- n1+n2
  p <- ncol(Y1)
  
  Y1.mean <- apply(Y1, 2, mean)
  Y2.mean <- apply(Y2, 2, mean)
  
  pooled <- ((n1-1)*apply(Y1, 2, var) + (n2-1)*apply(Y2, 2, var))/n
  
  Mn.vec <- (n1*n2/n) * (Y1.mean - Y2.mean)^2 / pooled
  
  Mn <- max(as.numeric(Mn.vec))
  
  cutoff <- -log(pi) - 2*log(-log(1-alpha)) + 2*log(p) - log(log(p))
  
  indicator <- ifelse(Mn >= cutoff, 1, 0)
  
  return(indicator)
}