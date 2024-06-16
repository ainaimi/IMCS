pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom
)


# define the inverse logistic function
expit <- function(x){
  1/(1 + exp(-x))
} 

set.seed(123)

n = 500000

z <- rnorm(n, mean = 0, sd = 1)

x <- rbinom(n, size = 1, p = expit(-1 + log(2)*z))

param_list <- c(.25, .5, .8, 1, 1.5, 2, 2.5, 3)

res <- NULL
for(i in param_list){
  
  y25 <- rbinom(n, size = 1, 
              p = expit(-log(1/.25 - 1) - log(i)*mean(x) - log(2)*0 
                        + log(i)*x + log(2)*z))
  
  y5 <- rbinom(n, size = 1, 
              p = expit(-log(1/.5 - 1) - log(i)*mean(x) - log(2)*0 
                        + log(i)*x + log(2)*z))
  
  y75 <- rbinom(n, size = 1, 
              p = expit(-log(1/.75 - 1) - log(i)*mean(x) - log(2)*0 
                        + log(i)*x + log(2)*z))
  
  y0 <- rbinom(n, size = 1, 
               p = expit(-1 + log(i)*x + log(2)*z))
  
  res <- rbind(res, 
               cbind(mean(y25), 
                     mean(y5), 
                     mean(y75), 
                     mean(y0)))
}

res