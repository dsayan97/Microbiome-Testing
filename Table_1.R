library(parallel)
library(foreach)
library(doParallel)
set.seed(1)

library(highmean) 
source("CLL.R") 
source("EBC.R")

ncores <- detectCores()/2
registerDoParallel(ncores)

cov.sim <- 1
setup.sim <- c(1,2)
s.sim <- 0.1
p.sim <- c(200,500,1000)
d.sim <- seq(0,1,0.2)

sim.grid1 <- expand.grid(50, 60, cov.sim, s.sim, d.sim,  p.sim, setup.sim)
sim.grid2 <- expand.grid(100, 120, cov.sim, s.sim, d.sim,  p.sim, setup.sim)
sim.grid <- rbind(sim.grid1,sim.grid2)
colnames(sim.grid) <- c("samp1sim", "samp2sim", "covsim", "ssim", "dsim", "psim", "setupsim")

outerfun = function(v){
  n1 <<- as.numeric(v[1])
  n2 <<- as.numeric(v[2])
  COV <<- as.numeric(v[3])
  s <<- as.numeric(v[4])
  d <<- as.numeric(v[5])
  p <<- as.numeric(v[6])
  SETUP <<- as.numeric(v[7])
  
  ITER <- 1e3
  alpha <- 0.05

  r_1 <- r_2 <- r_3 <- r_4 <- 0
  ########### Simulation
  
  for(iter in 1:ITER){
    source("Parameters_Two.R")
    
    res1 <- ifelse(EBC(Y1, Y2, N = 1e3)<alpha, 1, 0)
    r_1 <- c(r_1, res1)
    
    res2 <- CLL(Y1, Y2, alpha = alpha)
    r_2 <- c(r_2, res2)
    
    res3 <- ifelse(as.numeric(apval_Cai2014(Y1, Y2)$pval) < alpha, 1, 0)
    r_3 <- c(r_3, res3)
  }
  
  tab = data.frame("Sample1" = n1,
                   "Sample2" = n2,
                   "Dimension" = p,
                   "Covariane" = COV,
                   "Distribution" = SETUP,
                   "sparcity" = s,
                   "delta" = d,
                   "EBC" = round(sum(r_1)/ITER,3),
                   "CLL" = round(sum(r_2)/ITER,3),
                   "MECAF" = round(sum(r_3)/ITER,3))
  return(tab)
}

res <- foreach(i = 1:nrow(sim.grid))%dopar%{
  outerfun(sim.grid[i,])
}

res.df <- data.frame(matrix(unlist(res), byrow = T, ncol = length(res[[1]])))
colnames(res.df) <- colnames(res[[1]])

print(res.df)
