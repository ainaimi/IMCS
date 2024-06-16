pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

expit<-function(a){1/(1+exp(-a))}

set.seed(123)

collapsibility_function <- function(sample_size, intercept){
  n = sample_size
  
  C <- rbinom(n,size = 1,p = .5)
  
  theta <- c(0,log(2))
  pi_a <- expit(theta[1]+theta[2]*C) ## note the correction!
  
  A_obs <- rbinom(n, size = 1, p = pi_a)
  
  beta <- c(intercept, log(1.5), log(.5), log(3))
  mu <- expit(beta[1] + beta[2]*A_obs + beta[3]*C + beta[4]*A_obs*C)
  
  mu1_A1 <- mean(mu[A_obs==1])
  mu0_A1 <- mean(expit(beta[1] + beta[2]*0 + beta[3]*C[A_obs==1] + beta[4]*0*C[A_obs==1]))
  
  mu1_A0 <- mean(expit(beta[1] + beta[2]*1 + beta[3]*C[A_obs==0] + beta[4]*1*C[A_obs==0]))
  mu0_A0 <- mean(mu[A_obs==0])
  
  mu1 <- mean(expit(beta[1] + beta[2]*1 + beta[3]*C + beta[4]*1*C))
  mu0 <- mean(expit(beta[1] + beta[2]*0 + beta[3]*C + beta[4]*0*C))
  
  ## compute the odds ratios
  OR_A1 <- (mu1_A1/(1-mu1_A1))/(mu0_A1/(1-mu0_A1))
  OR_A0 <- (mu1_A0/(1-mu1_A0))/(mu0_A0/(1-mu0_A0))
  OR <- (mu1/(1-mu1))/(mu0/(1-mu0))
  
  ## output the marginally adjusted odds
  return(c(mean(mu),OR_A1, OR_A0, OR))
}


collapsibility_function(sample_size = 1e6, intercept = -2.5)
collapsibility_function(sample_size = 5e6, intercept = -2.5)
collapsibility_function(sample_size = 1e7, intercept = -2.5)
collapsibility_function(sample_size = 1e8, intercept = -2.5)