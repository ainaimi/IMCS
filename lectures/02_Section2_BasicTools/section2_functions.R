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


## -----------------------------------------------------------------------------

expit <- function(x){
  exp(x)/(1+exp(x))
  }


logit <- function(x){ 
  log(x/(1-x)) 
  }


## -----------------------------------------------------------------------------

expit(-2)

expit(-1.5)

logit(.2)

logit(.18)


## ----warning = F, message = F-------------------------------------------------

set.seed(123)

simulation_function <- function(nsim, sample_size, parameter){

  # data generation
  c <- rnorm(sample_size, mean = 0, sd = 1)
  
  p_y <- expit(-2 + log(parameter)*c)
  
  y <- rbinom(sample_size, size = 1, p = p_y)
  
  a_data <- data.frame(c, y)
  
  # analysis
  
  mY <- mean(a_data$y)
  
  glm_fit <- glm(y ~ c, data = a_data, family = binomial("logit"))
  
  glm_res <- summary(glm_fit)$coefficients[2,1]
  
  sim_res <- c(nsim, sample_size, parameter, mY, glm_res)
  
  return(sim_res)

}


## ----warning = F, message = F-------------------------------------------------

simulation_function(nsim = 1, sample_size = 500, parameter = 2)


## ----tidy = FALSE, warning = F, message = F-----------------------------------

simulation_results <- NULL

for(i in 1:10){
  simulation_results <- rbind(
    simulation_results, 
    simulation_function(nsim = i, sample_size = 500, parameter = 2)
  )
}

simulation_results


## ----tidy = FALSE, warning = F, message = F-----------------------------------

simulation_results <- list()

for(i in 1:10){
  simulation_results[[i]] <- simulation_function(nsim = i, sample_size = 500, parameter = 2)
}

simulation_results


## -----------------------------------------------------------------------------

library(data.table)

sim_res <- do.call(rbind, simulation_results)

sim_res


## ----tidy = FALSE, warning = F, message = F-----------------------------------

simulation_results <- lapply(1:10, function(x)
  simulation_function(nsim = x, sample_size = 500, parameter = 2)
)

simulation_results

do.call(rbind, simulation_results)


## ----tidy = FALSE, warning = F, message = F-----------------------------------

simulation_results <- sapply(1:10, function(x)
  simulation_function(nsim = x, sample_size = 500, parameter = 2), 
  simplify = T
)

simulation_results

# note the transpose!
t(simulation_results)


## ----tidy=F-------------------------------------------------------------------

simulation_results <- mapply(simulation_function,
                             nsim = 1:10, # number of simulations
                             sample_size = rep(c(250,500), each=10), # sample size
                             parameter = 2, # other function arguments
                             SIMPLIFY = T) # need to simplify for proper formatting

t(simulation_results)

sim_res <- t(simulation_results)

colnames(sim_res) <- c("index", "sample_size", "parameter", "meanY", "log_OR")

head(sim_res, 3)

tail(sim_res, 3)


## ----tidy=F-------------------------------------------------------------------

simulation_results <- mapply(simulation_function,
                             nsim = 1:10, # number of simulations
                             sample_size = rep(c(250, 500), each = 10), # sample size
                             parameter = rep(c(.5, 2), each = 20), # other function arguments
                             SIMPLIFY = T) # need to simplify for proper formatting

sim_res <- t(simulation_results)


## ----tidy=F-------------------------------------------------------------------

parm_data <- expand.grid(
  index = 1:10,
  n = c(250, 500, 1000),
  parms = c(.5, 2)
)


## ----tidy=F-------------------------------------------------------------------

simulation_results <- lapply(1:nrow(parm_data), function(x)
  simulation_function(nsim = parm_data[x,]$index,
                      sample_size = parm_data[x,]$n,
                      parameter = parm_data[x,]$parms)
)

sim_res <- do.call(rbind, simulation_results)

colnames(sim_res) <- c("index", "sample_size", "parameter", "meanY", "log_OR")

head(sim_res, 3)

tail(sim_res, 3)


## -----------------------------------------------------------------------------
set.seed(NULL)

## -----------------------------------------------------------------------------

a <- runif(5)

b <- rpois(5, lambda = 2)


## -----------------------------------------------------------------------------

a

b


## -----------------------------------------------------------------------------

set.seed(123)

a <- runif(5)

b <- rpois(5, lambda = 2)

a

b


## ----tidy = F, warning = F, message = F---------------------------------------

# re-initialize the seed 
set.seed(NULL)

set.seed(123)
a <- runif(50000)

set.seed(123)
b <- rpois(50000, lambda = 2)

## ----tidy = F, warning = F, message = F---------------------------------------
summary(lm(b ~ a))$coefficients

## ----tidy = F, warning = F, message = F---------------------------------------
plot_dat <- tibble(uniform = a, poisson = b)

ggplot(plot_dat, aes(x = uniform, y = poisson)) + 
  geom_point() +
  geom_smooth(method='lm', se = F)


## ----tidy = F, warning = F, message = F---------------------------------------

# re-initialize the seed 
set.seed(NULL)

set.seed(123)
a <- runif(50000)

#set.seed(123)
b <- rpois(50000, lambda = 2)

summary(lm(b ~ a))$coefficients

plot_dat <- tibble(uniform = a, poisson = b)

ggplot(plot_dat, aes(x = uniform, y = poisson)) + 
  geom_point() +
  geom_smooth(method='lm', se = F)


