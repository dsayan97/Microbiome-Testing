library(parallel)
library(foreach)
library(doParallel)
set.seed(1)

source("KEBC.R")

ncores <- detectCores()/2
registerDoParallel(ncores)

cov.sim <- 1
setup.sim <- c(1,2)
s.sim <- 0.1
p.sim <- c(200, 500, 1000)
d.sim <- seq(0,1,0.2)

sim.grid1 <- expand.grid(50, 60, 70, 80, 90, cov.sim, s.sim, d.sim,  p.sim, setup.sim)
sim.grid2 <- expand.grid(100, 120, 140, 160, 180, cov.sim, s.sim, d.sim,  p.sim, setup.sim)
sim.grid <- rbind(sim.grid1,sim.grid2)
colnames(sim.grid) <- c("samp1sim", "samp2sim","samp3sim", "samp4sim", "samp5sim", "covsim", "ssim", "dsim", "psim", "setupsim")

outerfun <- function(v){
  n1 <<- as.numeric(v[1])
  n2 <<- as.numeric(v[2])
  n3 <<- as.numeric(v[3])
  n4 <<- as.numeric(v[4])
  n5 <<- as.numeric(v[5])
  COV <<- as.numeric(v[6])
  s <<- as.numeric(v[7])
  d <<- as.numeric(v[8])
  p <<- as.numeric(v[9])
  SETUP <<- as.numeric(v[10])
  
  ITER <- 2
  alpha <- 0.05
  
  
  r_1 <- r_2 <- r_3 <- 0
  ########### Simulation
  
  for(iter in 1:ITER){
    source("Parameters_MANOVA.R")
    
    res1 <- ifelse(KEBC(list(Y1, Y2, Y3), N = 1e3)<alpha, 1, 0)
    r_1 <- c(r_1, res1)
    
    res2 <- ifelse(KEBC(list(Y1, Y2, Y3, Y4), N = 1e3)<alpha, 1, 0)
    r_2 <- c(r_2, res2)
    
    res3 <- ifelse(KEBC(list(Y1, Y2, Y3, Y4, Y5), N = 1e3)<alpha, 1, 0)
    r_3 <- c(r_3, res3)
  }
  
  tab <- data.frame("Sample1" = n1,
                   "Sample2" = n2,
                   "Sample3" = n3,
                   "Sample4" = n4,
                   "sample5" = n5,
                   "Dimension" = p,
                   "Covariane" = COV,
                   "Distribution" = SETUP,
                   "sparcity" = s,
                   "delta" = d,
                   "K_3" = round(sum(r_1)/ITER,3),
                   "K_4" = round(sum(r_2)/ITER,3),
                   "K_5" = round(sum(r_3)/ITER,3))
  return(tab)
}

res <- foreach(i = 1:nrow(sim.grid))%dopar%{
  outerfun(sim.grid[i,])
}

res.df <- data.frame(matrix(unlist(res), byrow = T, ncol = length(res[[1]])))
colnames(res.df) <- colnames(res[[1]])

print(res.df)