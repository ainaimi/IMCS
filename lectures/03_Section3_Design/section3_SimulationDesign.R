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


## ----tidy=F-------------------------------------------------------------------

rct_data <- matrix(
  c(53,193,139,350),
  ncol=2,
  byrow=T)

colnames(rct_data) <- c("event","nonevent")
rownames(rct_data) <- c("exposed","unexposed")
rct_data <- as.table(rct_data)
rct_data


## -----------------------------------------------------------------------------

risk_ratio <- (rct_data[1,1]/sum(rct_data[1,]))/(rct_data[2,1]/sum(rct_data[2,]))
round(risk_ratio,2)


## -----------------------------------------------------------------------------

SE_lnRR <- sqrt((1/rct_data[1,1] - 1/sum(rct_data[1,])) + (1/rct_data[2,1] - 1/sum(rct_data[2,])))


## -----------------------------------------------------------------------------

z <- (log(risk_ratio) - 0)/SE_lnRR
round(2*pnorm(-abs(z)),4)


## -----------------------------------------------------------------------------
rct_data
rct_data_long <- rbind(matrix(rep(c(1,1),rct_data[1,1]),ncol=2),
                       matrix(rep(c(0,1),rct_data[2,1]),ncol=2),
                       matrix(rep(c(1,0),rct_data[1,2]),ncol=2),
                       matrix(rep(c(0,0),rct_data[2,2]),ncol=2))
nrow(rct_data_long)

# re-shuffle rows
rct_data_long <- data.frame(rct_data_long[sample(nrow(rct_data_long)),])
names(rct_data_long) <- c("X","Y")

## -----------------------------------------------------------------------------
head(rct_data_long)

## ----tidy=F-------------------------------------------------------------------

set.seed(123)

rr_permuted <- NULL

permutations <- 20000

for(i in 1:permutations){
  
  permuted <- rct_data_long
  
  permuted$X <- permuted$X[sample(length(permuted$X))] # shuffle the exposure, N = 735
  
  res <- log(mean(subset(permuted,X==1 )$Y)/mean(subset(permuted,X==0)$Y)) # recalculate RR
  
  rr_permuted <- rbind(rr_permuted,res)
  
}

rr_permuted <- data.frame(rr_permuted)
names(rr_permuted) <- "estimates"

## ----warning=F, message=F, out.width = "5cm",fig.cap="Distribution of log risk ratios after 2,000 random permutations of the exposure variable in the 2x2 table data above. The solid blue density curve represents a nonparametric kernel density estimate of the distribution. The solid red density curve represents a normal density estimate of the distribution. The dashed red vertical line indicates the value of the log risk ratio estimated in the original unpermuted data.",echo=F----
ggplot(rr_permuted) +  
  geom_histogram(aes(estimates,y=..density..),color="gray",fill="white") + 
  geom_density(aes(estimates),color="blue") + 
  stat_function(
    fun = dnorm, 
    args = with(rr_permuted, c(mean = mean(estimates), sd = sd(estimates))),
    color="red"
  ) + 
  geom_vline(xintercept = log(risk_ratio),color="red",linetype=2)

## -----------------------------------------------------------------------------
sum(rr_permuted$estimate <= log(risk_ratio))

## -----------------------------------------------------------------------------
sum(rr_permuted$estimate <= log(risk_ratio))/permutations

## -----------------------------------------------------------------------------
sum(abs(rr_permuted$estimate) >= abs(log(risk_ratio)))/permutations

## ----out.width = "10cm",fig.cap="Simple confounding triangle, with exposure $A$, confounder $C$, and outcome $Y$.",echo=F----
knitr::include_graphics(here("_images","triangle_dag.pdf"))

## ----echo=T,fig.star=T,tidy=F,highlight=T-------------------------------------
  
  expit <- function(a){1/(1+exp(-a))}

  set.seed(123)

  n = 500
  # simulate confounder from a normal distribution
  ## first argument is sample size
  ## second argument is mean
  ## third argument is SD
  C <- rnorm(n,0,1)

  # propensity model
  ## theta is a 2-dimensional vector (list) of parameters for the exposure model
  ## theta[1] is the intercept and theta[2] is the log-OR for the confounder-exposure relation
  ## pi is the probability that the exposure is 1
  theta <- c(0,log(2))
  pi <- expit(theta[1]+theta[1]*C)
  
  # simulate exposure from binomial distribution
  ## second argument is number of trials
  ## third argument is probability that A = 1
  A <- rbinom(n,1,pi)

  # outcome model
  ## beta is a 3-dimensional vector (list) of parameters for the outcome model
  ## beta[1] is the intercept, beta[2] is the exp(OR) for the exposure-outcome relation
  ## beta[3] is the log-OR for the confounder-outcome relation
  beta <- c(-2.75,log(2),log(2))
  mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
  Y <- rbinom(n,1,mu)


## ----echo=T,fig.star=T,tidy=F,highlight=T-------------------------------------
  # glm.res0 is the outcome's prevalence
  ## we store this to ouput it from the function
  glm.res0 <- mean(Y)

  print(glm.res0)

## ----echo=T,fig.star=T,tidy=F,highlight=T-------------------------------------

head(data.frame(Y,A,C=round(C,2)),3)

tail(data.frame(Y,A,C=round(C,2)),3)


## ----echo=T,fig.star=T,tidy=F,highlight=T-------------------------------------

  # estimate the true exposure odds ratio using a conditonally adjusted logit model
  ## NB: conditionally adjusted logit model is not the same as "conditional logistic regression"
  ## glm.res1 is the conditional log-odds ratio
  
  m1 <- glm(Y~A+C,family=binomial(link="logit"))
  glm.res1 <- m1$coefficients[2]
  

## ----echo=T,fig.star=T,tidy=F,highlight=T-------------------------------------

  # estimate the true exposure odds ratio using a marginally adjusted logit model
  ## compute the average predicted probabilities under A = 1 and then A = 0
  muhat1 <- mean(predict(m1,newdata=data.frame(A=1,C),type="response"))
  muhat0 <- mean(predict(m1,newdata=data.frame(A=0,C),type="response"))
  
  ## compute the odds from these average probabilities
  odds1 <- muhat1/(1-muhat1)
  odds0 <- muhat0/(1-muhat0)
  
  ## glm.res2 is the marginal log-odds ratio
  glm.res2 <- log(odds1/odds0)


## ----tidy = F, warning = F, message = F---------------------------------------

  expit<-function(a){1/(1+exp(-a))}

  set.seed(123)
  
  collapsibility_function <- function(index, intercept){
      n=500
      
      C <- rnorm(n,0,1)
    
      theta <- c(0,log(2))
      pi <- expit(theta[1]+theta[1]*C)
      
      A <- rbinom(n,1,pi)
    
      beta <- c(intercept,log(2),log(2))
      mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
      Y <- rbinom(n,1,mu)
    
      glm.res0 <- mean(Y)
    
      m1 <- glm(Y~A+C,family=binomial(link="logit"))
      glm.res1 <- m1$coefficients[2]
      
      muhat1 <- mean(predict(m1,newdata=data.frame(A=1,C),type="response"))
      muhat0 <- mean(predict(m1,newdata=data.frame(A=0,C),type="response"))
      
      ## compute the odds from these average probabilities
      odds1 <- muhat1/(1-muhat1)
      odds0 <- muhat0/(1-muhat0)
      
      ## glm.res2 is the marginal log-odds ratio
      glm.res2 <- log(odds1/odds0)
      
      res <- data.frame(prevalenceY = glm.res0, 
                        conditionalOR = glm.res1, 
                        marginalOR = glm.res2)
      
      return(res)
  }
  
  sim_res <- lapply(1:2000, function(x) collapsibility_function(index = x, intercept = -2.75))
  
  sim_res <- do.call(rbind, sim_res)
  
  head(sim_res)
  
  mean(sim_res$conditionalOR)
  
  mean(sim_res$marginalOR)


## ----out.width = "10cm",fig.cap="Simple confounding triangle, with exposure $A$, confounder $C$, and outcome $Y$.",echo=F----
knitr::include_graphics(here("_images","triangle_dag.pdf"))

## ----eval=F-------------------------------------------------------------------
#  
#  collapsibility_function <- function(index, intercept){
#  
#        n=500
#  
#        C <- rnorm(n,0,1) ## FIRST!!!
#  
#        theta<- c(0,log(2))
#        pi <- expit(theta[1]+theta[1]*C)
#  
#        A <- rbinom(n,1,pi) ## SECOND!!
#  
#        beta <- c(intercept,log(2),log(2))
#        mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
#        Y <- rbinom(n,1,mu) ## THIRD!
#  
#        ...
#  
#  }
#  

## ----out.width = "10cm",fig.cap="Complex mediation diagram with unmeasured confounder $U$, baseline confounders $C$, mediator-outcome confounder affected by the exposure $L$, mediator $Z$, exposure $X$, and outcome $Y$.",echo=F----
knitr::include_graphics(here("_images","mediation_dag.pdf"))

## -----------------------------------------------------------------------------

expit <- function(x){1/(1 + exp(-x))}

set.seed(123)

n = 500

U <- rnorm(n)

C <- rnorm(n)

X <- rbinom(n, size = 1, p = expit(-1 + log(2)*C ))

L <- rnorm(n, mean = 0 + 1.5*X + 1.5*U)

Z <- rbinom(n, size = 1, expit(-1.5 + log(2)*X + log(2)*L ))

Y <- rnorm(n, mean = 10 + 5*X + 5*C + 5*Z + 4*L + 4*U, sd = 5)

med_data <- data.frame(Y, Z, L, X, C)

head(med_data)


## ----eval = F-----------------------------------------------------------------
#  
#    # outcome model
#    ## beta is a 3-dimensional vector (list) of parameters for the outcome model
#    ## beta[1] is the intercept, beta[2] is the exp(OR) for the exposure-outcome relation
#    ## beta[3] is the log-OR for the confounder-outcome relation
#    beta <- c(-2.75,log(2),log(2))
#    mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
#    Y <- rbinom(n,1,mu)
#  

## -----------------------------------------------------------------------------

expit<-function(a){1/(1+exp(-a))}

set.seed(123)
  
collapsibility_function <- function(intercept, exposure){
      n = 5e6 # five million observations
      
      C <- rnorm(n,0,1)
    
      theta <- c(0,log(2))
      pi <- expit(theta[1]+theta[1]*C)
      
      A <- exposure # set the exposure to a specific value
    
      beta <- c(intercept,log(2),log(2))
      mu <- expit(beta[1] + beta[2]*A + beta[3]*C)
      Y <- rbinom(n,1,mu)
    
      mu_ <- mean(Y) # output the mean of Y under the specific exposure value
    
      ## compute the odds from these average probabilities
      odds <- mu_/(1-mu_)

      ## output the marginally adjusted odds
      return(odds)
  }


## -----------------------------------------------------------------------------

odds1 <- collapsibility_function(intercept = -2.75, exposure = 1)

odds1


## -----------------------------------------------------------------------------

odds0 <- collapsibility_function(intercept = -2.75, exposure = 0)

odds0


## -----------------------------------------------------------------------------

true_ORm <- odds1/odds0

true_ORm


