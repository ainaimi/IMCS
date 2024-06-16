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


## ----tidy=F, warning = F, message = F-----------------------------------------

expit <- function(x){
  exp(x)/(1+exp(x))
  }

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


## ----tidy=F-------------------------------------------------------------------
simulation_function(nsim = 1, sample_size = 500, parameter = 2)

## ----tidy=F-------------------------------------------------------------------

#install.packages("microbenchmark")

library(microbenchmark)

microbenchmark(
  sapply_results1 = sapply(1:100, 
                           function(x) simulation_function(nsim = x, 
                                                           sample_size = 500, 
                                                           parameter = 2), 
                           simplify = T),
  lapply_results2 = lapply(1:100, 
                           function(x) simulation_function(nsim = x, 
                                                           sample_size = 500, 
                                                           parameter = 2)),
  times=10L
  )


## ----eval = F-----------------------------------------------------------------
#  lapply(1:100, function(x) simulation_function(nsim = x,
#                                                sample_size = 500,
#                                                parameter = 2))

## ----out.width = "10cm",fig.cap="Simple illustration of serial versus parallel computing.",echo=F----
knitr::include_graphics(here("_images","par_ser.pdf"))

## -----------------------------------------------------------------------------

library(parallel)

detectCores(logical = TRUE)

detectCores(logical = FALSE)


## ----tidy=FALSE---------------------------------------------------------------

expit<-function(a){1/(1+exp(-a))}

set.seed(123)

collapsibility_function <- function(index, intercept){
    n=500
    
    C<-rnorm(n,0,1)
  
    theta<- c(0,log(2))
    pi <- expit(theta[1]+theta[1]*C)
    
    A<-rbinom(n,1,pi)
  
    beta<-c(intercept,log(2),log(2))
    mu<-expit(beta[1] + beta[2]*A + beta[3]*C)
    Y<-rbinom(n,1,mu)
  
    glm.res0<-mean(Y)
  
    m1<-glm(Y~A+C,family=binomial(link="logit"))
    glm.res1<-m1$coefficients[2]
    
    muhat1<-mean(predict(m1,newdata=data.frame(A=1,C),type="response"))
    muhat0<-mean(predict(m1,newdata=data.frame(A=0,C),type="response"))
    
    ## compute the odds from these average probabilities
    odds1<-muhat1/(1-muhat1)
    odds0<-muhat0/(1-muhat0)
    
    ## glm.res2 is the marginal log-odds ratio
    glm.res2<-log(odds1/odds0)
    
    res <- data.frame(intercept = intercept,
                      prevalenceY = glm.res0, 
                      conditionalOR = glm.res1, 
                      marginalOR = glm.res2)
    
    return(res)
}

## ----tidy = FALSE-------------------------------------------------------------

# choose the number of cores to use
num_cores <- detectCores() - 2

# how many cores?
num_cores

# special mclapply seed ... 
RNGkind("L'Ecuyer-CMRG")

# set the seed
set.seed(123)

# run the function
sim_res <- mclapply(1:50000, 
                    function(x) collapsibility_function(index = x, 
                                                        intercept = -2.75), 
                    mc.cores = num_cores)

sim_res0 <- do.call(rbind, sim_res)

sim_res <- mclapply(1:50000, 
                    function(x) collapsibility_function(index = x, 
                                                        intercept = 0), # change the intercept! 
                    mc.cores = num_cores)

sim_res1 <- do.call(rbind, sim_res)

sim_res <- rbind(sim_res0,
                 sim_res1)

head(sim_res, 3)

tail(sim_res, 3)


## ----tidy=FALSE---------------------------------------------------------------

sim_res %>% 
  group_by(intercept) %>% 
  summarize(mY = mean(prevalenceY),
            m_cOR = mean(conditionalOR),
            m_mOR = mean(marginalOR))


## ----tidy = FALSE-------------------------------------------------------------

# define the cluster type
cluster_type <- "FORK"

# number of cores
n.cores <-  detectCores() - 2

# make the cluster
cluster = makeCluster(n.cores,
                      type = cluster_type)

clusterSetRNGStream(cl = cluster, iseed = 123)

sim_res <- parLapply(cluster,
                     1:50000, 
                     function(x) collapsibility_function(index = x, 
                                                         intercept = -2.75)
)

sim_res0 <- do.call(rbind, sim_res)

sim_res <- parLapply(cluster,
                     1:50000, 
                     function(x) collapsibility_function(index = x, 
                                                         intercept = 0)
)

sim_res1 <- do.call(rbind, sim_res)

stopCluster(cl = cluster)

sim_res <- rbind(sim_res0,
                 sim_res1)

head(sim_res, 3)

tail(sim_res, 3)


## ----tidy=FALSE---------------------------------------------------------------

sim_res %>% 
  group_by(intercept) %>% 
  summarize(mY = mean(prevalenceY),
            m_cOR = mean(conditionalOR),
            m_mOR = mean(marginalOR))


## ----tidy = FALSE-------------------------------------------------------------

cluster_type <- "PSOCK"

n.cores <- detectCores()

cluster = makeCluster(n.cores,
                      type = cluster_type)

clusterSetRNGStream(cl = cluster, iseed = 123)

clusterExport(cluster, 
              c("collapsibility_function",
                "expit")
              )

sim_res <- parLapply(cluster,
                      1:50000, 
                      function(x) collapsibility_function(index = x, 
                                                         intercept = 0)
)

sim_res <- do.call(rbind, sim_res)

stopCluster(cl = cluster)

head(sim_res, 3)



## ----eval = F-----------------------------------------------------------------
#  
#  parallel::clusterEvalQ(cluster, {
#    .libPaths("~/R/R_LIBS_USER")
#    library(SuperLearner)
#    library(ranger)
#    library(xgboost)
#    library(glmnet)
#  })
#  

