## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
library(formatR)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR",
               "gridExtra","skimr","here","RColorBrewer","survival")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

remotes::install_github("rstudio/fontawesome")

library(fontawesome)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)


## ----tidy = F, warning = F, message = F, eval = F-----------------------------
#  
#  expit <- function(a){1/(1+exp(-a))}
#  
#  set.seed(123)
#  collapsibility_function <- function(index, intercept, exposure){
#      n = 500000000
#  
#      C <- rnorm(n,0,1)
#  
#      theta <- c(0,log(2))
#      pi <- expit(theta[1]+theta[1]*C)
#  
#      A <- exposure
#  
#      beta <- c(intercept,log(2),log(2))
#      mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
#      Y <- rbinom(n,1,mu)
#  
#      res <- mean(Y)/(1 - mean(Y))
#  
#      return(res)
#  }
#  
#  odds1 <- collapsibility_function(index=1, intercept = 0, exposure = 1)
#  odds0 <- collapsibility_function(index=1, intercept = 0, exposure = 0)
#  
#  or_marg <- odds1/odds0
#  
#  # true marginal OR for intercept = 0
#  # 1.871259
#  

## -----------------------------------------------------------------------------

library(parallel)
library(boot)

expit <- function(a){1/(1+exp(-a))}

set.seed(123)

collapsibility_function <- function(index, intercept, true_m, true_c){
    n = 500
    
    C <- rnorm(n,0,1)
  
    theta <- c(0,log(2))
    pi <- expit(theta[1]+theta[1]*C)
    
    A <- rbinom(n,1,pi)
  
    beta <- c(intercept,log(2),log(2))
    mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
    Y <- rbinom(n,1,mu)
  
    glm.res0 <- mean(Y)
  
    m1 <- glm(Y~A+C,family=binomial(link="logit"))
    glm.res1 <- summary(m1)$coefficients[2,1:2]
    
    muhat1 <- mean(predict(m1,newdata=data.frame(A=1,C),type="response"))
    muhat0 <- mean(predict(m1,newdata=data.frame(A=0,C),type="response"))
    
    ## compute the odds from these average probabilities
    odds1 <- muhat1/(1-muhat1)
    odds0 <- muhat0/(1-muhat0)
    
    ## glm.res2 is the marginal log-odds ratio
    glm.res2 <- log(odds1/odds0)
    
    # bootstrap SEs
    c_data <- data.frame(Y, A, C)
    boot_func <- function(data, index){ ## order matters!
      
      boot_dat <- data[index, ]
      
      m1_ <- glm(Y ~ A + C, data = boot_dat, family=binomial(link="logit"))
      glm.res1 <- summary(m1_)$coefficients[2,1:2]
    
      muhat1_ <- mean(predict(m1_, newdata = transform(boot_dat, A=1), type="response"))
      muhat0_ <- mean(predict(m1_, newdata = transform(boot_dat, A=0), type="response"))
    
      ## compute the odds from these average probabilities
      odds1_ <- muhat1_/(1-muhat1_)
      odds0_ <- muhat0_/(1-muhat0_)
      
      return(odds1_/odds0_)
    }
    
    boot_obj <- boot(data = c_data,
                     statistic = boot_func, 
                     R = 200, 
                     parallel = "no") # can only parallelize one level
    
    glm.res2 <- c(glm.res2, sd(boot_obj$t))
    
    res <- data.frame(intercept = intercept,
                      t(glm.res1), 
                      t(glm.res2), 
                      true_m = true_m,
                      true_c = true_c)
    
    return(res)
}

# choose the number of cores to use
# num_cores <- detectCores() - 2

# how many cores?
# num_cores

# special mclapply seed ... 
# RNGkind("L'Ecuyer-CMRG")

# set the seed
set.seed(123)

# run the function
sim_res <- lapply(1:500, 
                  function(x) collapsibility_function(index = x, 
                                                      intercept = 0,
                                                      true_m = 1.871259,
                                                      true_c = 2))#, 
                    #mc.cores = num_cores)

sim_res <- do.call(rbind, sim_res)


## -----------------------------------------------------------------------------

names(sim_res) <- c("intercept",
                    "cEstimate", "cSE",
                    "mEstimate", "mSE",
                    "true_m", "true_c")

head(sim_res, 3)

tail(sim_res, 3)


## ----tidy = FALSE-------------------------------------------------------------

mc_bias <- sim_res %>%
  summarize(bias_c = mean(cEstimate - log(true_c)),
            bias_m = mean(mEstimate - log(true_m)))

mc_bias


## -----------------------------------------------------------------------------

mc_se_bias <- function(x, n){
  sqrt(sum((x - mean(x))^2)/(n*(n-1)))
}

mc_bias_se <- c(mc_se_bias(sim_res$cEstimate, n = 500),
                mc_se_bias(sim_res$mEstimate, n = 500))

mc_bias_se

## ----tidy = FALSE-------------------------------------------------------------

mse_c = mean((sim_res$cEstimate - log(sim_res$true_c))^2)
mse_m = mean((sim_res$mEstimate - log(sim_res$true_m))^2)


## -----------------------------------------------------------------------------

mc_se_mse <- function(x, n){
  sqrt(sum(((x - mean(x))^2 - mean((x - mean(x))^2))^2)/(n*(n-1)))
}

mc_se_mse(sim_res$cEstimate, n = 500)

mc_se_mse(sim_res$mEstimate, n = 500)


## -----------------------------------------------------------------------------
head(sim_res, 3)

## -----------------------------------------------------------------------------

sd(sim_res$cEstimate) - mean(sim_res$cSE)

sd(sim_res$mEstimate) - mean(sim_res$mSE)

sd(sim_res$cEstimate)/mean(sim_res$cSE)

sd(sim_res$mEstimate)/mean(sim_res$mSE)


## ----tidy=FALSE---------------------------------------------------------------

sim_res <- sim_res %>% mutate(
  cLCL = cEstimate - 1.96*cSE,
  cUCL = cEstimate + 1.96*cSE,
  mLCL = mEstimate - 1.96*mSE,
  mUCL = mEstimate + 1.96*mSE
)


## -----------------------------------------------------------------------------

head(sim_res, 3)


## ----tidy=FALSE---------------------------------------------------------------

sim_res <- sim_res %>% mutate(
  cCoverage = cLCL < log(true_c) & log(true_c) < cUCL,
  mCoverage = mLCL < log(true_m) & log(true_m) < mUCL
)

mean(sim_res$cCoverage)

mean(sim_res$mCoverage)


## ----tidy=FALSE---------------------------------------------------------------

sim_res <- sim_res %>% mutate(
  cCI_length = cUCL - cLCL,
  mCI_length = mUCL - mLCL
)

mean(sim_res$cCI_length)

mean(sim_res$mCI_length)


## ----tidy = FALSE-------------------------------------------------------------

sim_res <- sim_res %>% mutate(
  cPower = cLCL > 0 | cUCL < 0,
  mPower = mLCL > 0 | mUCL < 0
)

mean(sim_res$cPower)

mean(sim_res$mPower)


