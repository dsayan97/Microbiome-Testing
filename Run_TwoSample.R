library(highmean) #The apval_Cai2014 function shares the same statistic as MECAF
source("CLL.R") #Cao-Lin-Li Test
source("EBC.R") #Proposed Test

  ########### Sample Parameters  
  n1 <- 100 #1st sample size #Table values: 50 and 100
  n2 <- 120 #2nd sample size #Table values: 60 and 120
  p <- 500 #dimension #Table values: 200, 500 and 1000
  s <- 0.1 #sparsity parameter 
  d <- 0.6 # separation parameter #Table values: 0 (size), 0.2, 0.4, 0.6, 0.8, 1
  
  ########### Simulation Parameters
  alpha <- 0.05 #level of significance
  ITER <- 1e3 #number of iterations
  
  ########### Simulated Data
  COV <- 1 #Table values: 1,2,3
  SETUP <- 1 ##Table values: 1(Gaussian), 2(t4)
  
  ########### Simulation
  r_1 <- r_2 <- r_3 <- 0
  for(iter in 1:ITER){
    source("Parameters_Two.R")
    cat("Iteration:", iter, "\n")
    
      res1 <- ifelse(EBC(Y1, Y2, N = 1e3)<alpha, 1, 0)
      r_1 <- c(r_1, res1)
    
      res2 <- CLL(Y1, Y2, alpha = alpha)
      r_2 <- c(r_2, res2)
      
      res3 <- ifelse(as.numeric(apval_Cai2014(Y1, Y2)$pval) < alpha, 1, 0)
      r_3 <- c(r_3, res3)
  }
  
  tab <- data.frame("EBC" = round(mean(r_1),3),
                   "CLL" = round(mean(r_2),3),
                   "MECAF" = round(mean(r_3),3))
  
  print(tab)