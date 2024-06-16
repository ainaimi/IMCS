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


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

n = 5000

x <- rnorm(n, mean = 0, sd = 1)

y <- 5 + 2*x + rnorm(n, mean = 0, sd = 1)

a <- data.frame(x, y)

head(a)


## ----tidy = F, warning = F, message = F---------------------------------------

GGally::ggpairs(a)

summary(a$y)

summary(a$x)


## ----tidy = F, warning = F, message = F---------------------------------------

mod1 <- lm(y ~ x, data = a)

# mod1 <- glm(y ~ x, data = a, family = gaussian(link = "identity"))

summary(mod1)$coefficients


## ----tidy = F, warning = F, message = F---------------------------------------

# define the inverse logistic function
expit <- function(x){
  1/(1 + exp(-x))
} 

set.seed(123)

n = 5000

z <- rnorm(n, mean = 0, sd = 1)

x <- rbinom(n, size = 1, p = expit(-1 + log(2)*z))

y <- 100 + 10*x + 3*z + rnorm(n, mean = 0, sd = 10)

# use these variables to construct a dataset:

a <- data.frame(z, x, y)


## ----tidy = F, warning = F, message = F---------------------------------------

head(a)


## ----tidy = F, warning = F, message = F---------------------------------------

ggplot(a) + geom_histogram(aes(x = y))

ggplot(a) + geom_histogram(aes(x = z))

table(a$x)

mean(a$x)


## ----tidy = F, warning = F, message = F---------------------------------------

# mod1 <- lm(y ~ x, data = a)

mod1 <- glm(y ~ x + z, data = a, family = gaussian(link = "identity"))

summary(mod1)$coefficients


## ----warning = F, message = F, eval = F---------------------------------------
#  
#  mod1 <- glm(y ~ x + z, data = a, family = gaussian(link = "identity"))
#  
#  summary(mod1)
#  

## ----warning = F, message = F, eval = F---------------------------------------
#  
#  y <- 100 + 10*x + 3*z + rnorm(n, mean = 0, sd = 10)
#  

## ----tidy = F, warning = F, message = F---------------------------------------
a$propensity_score <- glm(x ~ z, data = a, family = binomial("logit"))$fitted.values

## ----tidy = F, warning = F, message = F---------------------------------------
x <- rbinom(n, size = 1, p = expit(-1 + log(2)*z))

## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

n = 5000

z <- rnorm(n, mean = 0, sd = 1)

x <- rbinom(n, size = 1, p = expit(-1 + log(2)*z))

y <- rpois(n, lambda = exp(2 + log(2)*x + log(1.5)*z))

# use these variables to construct a dataset:

b <- data.frame(z, x, y)

head(b)


## ----tidy = F, warning = F, message = F---------------------------------------

GGally::ggpairs(b)


## ----tidy = F, warning = F, message = F---------------------------------------

ggplot(b) + 
  geom_histogram(aes(x = y)) + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

mod1 <- glm(y ~ x + z, data = b, family = poisson(link = "log"))

summary(mod1)$coefficients


## ----tidy = F, warning = F, message = F---------------------------------------

library(boot)

#' Regress the outcome against the exposure and covariate
ms_model <- glm(y ~ x + z, data = a, family = gaussian(link = "identity"))

##' Generate predictions for everyone in the sample to obtain 
##' unexposed (mu0 predictions) and exposed (mu1 predictions) risks.
mu1 <- predict(ms_model, newdata = transform(a, x=1), type="response")
mu0 <- predict(ms_model, newdata = transform(a, x=0), type="response")

#' Mean difference in predicted outcomes
marg_stand_MD <- mean(mu1) - mean(mu0)

#' Using the bootstrap to obtain confidence intervals for the marginally adjusted 
#' mean difference.
bootfunc <- function(data,index){
  boot_dat <- data[index,]
  ms_model <- glm(y ~ x + z, data=boot_dat, family = gaussian(link = "identity"))
  mu1 <- predict(ms_model, newdata = transform(boot_dat,x=1), type="response")
  mu0 <- predict(ms_model, newdata = transform(boot_dat,x=0), type="response")
  
  marg_stand_MD_ <- mean(mu1) - mean(mu0)
  
  return(marg_stand_MD_)
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(a, bootfunc, R=2000)

boot_MD <- boot.ci(boot_res, type = "norm")

marg_stand_MD

boot_MD


## ----tidy = F, warning = F, message = F---------------------------------------

# create the propensity score in the dataset
a$propensity_score <- glm(x ~ z, data = a, family = binomial("logit"))$fitted.values

# stabilized inverse probability weights
a$sw <- (mean(a$x)/a$propensity_score)*a$x + 
  ((1-mean(a$x))/(1-a$propensity_score))*(1-a$x)

summary(a$sw)

head(a)


## ----tidy = F, warning = F, message = F---------------------------------------

mod_MD_weighted <- glm(y ~ x, data = a, weights=sw, family = gaussian("identity"))

summary(mod_MD_weighted)$coefficients


## ----tidy = F, warning = F, message = F---------------------------------------

library(lmtest)
library(sandwich)

coeftest(mod_MD_weighted, 
         vcov. = vcovHC(mod_MD_weighted, type = "HC3"))[2,]

coefci(mod_MD_weighted, 
       level = 0.95, 
       vcov. = vcovHC(mod_MD_weighted, type = "HC3"))[2,]


## ----tidy = F, warning = F, message = F---------------------------------------
#' Using the bootstrap to obtain confidence intervals for the IP weighted 
#' mean difference.
bootfunc <- function(data,index){
  
  boot_dat <- data[index,]

  boot_dat$propensity_score <- glm(x ~ z, data = boot_dat, family = binomial("logit"))$fitted.values

  # stabilized inverse probability weights
  boot_dat$sw <- (mean(boot_dat$x)/boot_dat$propensity_score)*boot_dat$x + 
    ((1-mean(boot_dat$x))/(1-boot_dat$propensity_score))*(1-boot_dat$x)  
  
  mod_MD_weighted_ <- glm(y ~ x, data = boot_dat, weights=sw, family = gaussian("identity"))

  res <- summary(mod_MD_weighted_)$coefficients[2,1]
  
  return(res)
  
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(a, bootfunc, R = 2000)

boot_IP_weight <- boot.ci(boot_res, type = "norm")

summary(mod_MD_weighted)$coefficients[2,1]

boot_IP_weight

## ----tidy = F, warning = F, message = F, include = F--------------------------

## marginal standardization

marg_stand_res <- c(marg_stand_MD, boot_MD$normal[2:3])

ip_weighted_res <- c(summary(mod_MD_weighted)$coefficients[2,1], 
                     coefci(mod_MD_weighted, 
                            vcov. = vcovHC(mod_MD_weighted, type = "HC3"))[2,])

ip_weighted_boot <- c(summary(mod_MD_weighted)$coefficients[2,1], 
                      boot_IP_weight$normal[2:3])

## ----tidy = F, message = F, warning = F, echo = F-----------------------------
kable(
  rbind(marg_stand_res,
      ip_weighted_res,
      ip_weighted_boot),
      col.names = c("Estimate", "LCL", "UCL")
)

