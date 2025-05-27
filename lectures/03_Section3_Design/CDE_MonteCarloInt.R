pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom
)

# goal is to estimate true CDE in complex DAG setting


expit <- function(x) {
  1/(1 + exp(-x))
  }

set.seed(123)

cde_function <- function(sample_size, exposure, mediator, exp_med_inter){
  ## starting here
  n = sample_size
  
  U <- rnorm(n)
  
  C <- rnorm(n)
  
  X <- exposure #rbinom(n, size = 1, p = expit(-1 + log(2) * C))
  
  L <- rnorm(n, mean = 0 + 1.5 * X + 1.5 * U)
  
  Z <- mediator #rbinom(n, size = 1, expit(-1.5 + log(2) * X + log(2) * L))
  
  Y <- rnorm(n, mean = 10 + 5*X + 5*C + 2.5*X*C + 5*Z + exp_med_inter*X*Z + 4*L + 4*U, sd = 5)
  
  mY <- mean(Y)
  
  res_out <- c(mY, sample_size, exp_med_inter)
  
  return(res_out)
  
  # CDE: E(Y^x=1,z=0) - E(Y^x=0,z=0)  
}

muY10 <- cde_function(sample_size = 1e7, exposure = 1, mediator = 0, exp_med_inter = 2.5)

muY10

muY00 <- cde_function(sample_size = 1e7, exposure = 0, mediator = 0, exp_med_inter = 2.5)

muY00

# CDE for scenario 1: 10.99163

muY10[1] - muY00[1]

muY11 <- cde_function(sample_size = 1e7, exposure = 1, mediator = 1, exp_med_inter = 2.5)

muY11

muY01 <- cde_function(sample_size = 1e7, exposure = 0, mediator = 1, exp_med_inter = 2.5)

muY01

# CDE for scenario 2: 13.50302

muY11[1] - muY01[1]