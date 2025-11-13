source("KEBC.R") #Proposed Test

########### Sample Parameters  
n1 <- 100 #1st sample size #Table values: 50 and 100
n2 <- 120 #2nd sample size #Table values: 60 and 120
n3 <- 140 #3rd sample size #Table values: 70 and 140
n4 <- 160 #4th sample size #Table values: 80 and 160
n5 <- 180 #5th sample size #Table values: 90 and 180
p <- 500 #dimension #Table values: 200, 500 and 1000
s <- 0.1 #sparsity parameter 
d <- 0.6 # separation parameter #Table values: 0 (size), 0.2, 0.4, 0.6, 0.8, 1

########### Simulation Parameters
alpha <- 0.05 #level of significance
ITER <- 1e3 #number of iterations

########### Simulated Data
SETUP <- 1 #Table values: 1(Gaussian), 2(t4)

########### Simulation
r_1 <- r_2 <- r_3 <- 0
for(iter in 1:ITER){
  source("Parameters_MANOVA.R")
  cat("Iteration:", iter, "\n")
  
  res1 <- ifelse(KEBC(list(Y1, Y2, Y3), N = 1e3)<alpha, 1, 0)
  r_1 <- c(r_1, res1)
  
  res2 <- ifelse(KEBC(list(Y1, Y2, Y3, Y4), N = 1e3)<alpha, 1, 0)
  r_2 <- c(r_2, res2)
  
  res3 <- ifelse(KEBC(list(Y1, Y2, Y3, Y4, Y5), N = 1e3)<alpha, 1, 0)
  r_3 <- c(r_3, res3)
}

tab <- data.frame("K_3" = round(mean(r_1),3),
                 "K_4" = round(mean(r_2),3),
                 "K_5" = round(mean(r_3),3))

print(tab)