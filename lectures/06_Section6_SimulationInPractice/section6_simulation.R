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


## ----out.width = "10cm",fig.cap="Simple confounding triangle, with exposure $A$, confounder $C$, and outcome $Y$.",echo=F----
knitr::include_graphics(here("_images","triangle_dag.pdf"))

## ----eval = F-----------------------------------------------------------------
#  
#  p_model <- glm(x ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10, data = the_data, family = binomial(link = "logit"))
#  

## ----eval = F-----------------------------------------------------------------
#  
#  p_model$fitted.values
#  

## ----eval = F-----------------------------------------------------------------
#  
#  1 - p_model$fitted.values
#  

## ----eval = F-----------------------------------------------------------------
#  
#  p_model <- glm(x ~ ., data = the_data, family = binomial(link = "logit"))
#  

## ----eval = F-----------------------------------------------------------------
#  
#  sw_trim <- ifelse(sw > quantile(sw, .975),
#                    quantile(sw, .975),
#                    analysis_data$sw)
#  

## ----eval = F-----------------------------------------------------------------
#  
#  mu_model <- glm(y ~ ., data = the_data, family = binomial(link = "logit"))
#  

## ----eval = F-----------------------------------------------------------------
#  
#  muhat1 <- mean(predict(mu_model, newdata = transform(the_data, x=1), type="response"))
#  

## ----tidy = F, eval = F-------------------------------------------------------
#  
#   expit<-function(a){1/(1+exp(-a))}
#  
#   set.seed(123)
#   simulation_function_true <- function(intercept = -2,
#                                        exposure = 1,
#                                        c_number = 10,
#                                        cov_mat = 0,
#                                        diag_cov = 1){
#       n = 5e6
#  
#       # how many confounders?
#       p <- c_number
#  
#       ## confounder matrix
#       sigma <- matrix(cov_mat, nrow=p, ncol=p)
#       diag(sigma) <- diag_cov
#       c <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sigma)
#  
#       # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
#       piMat <- model.matrix(
#         as.formula(
#           paste("~(",
#                 paste("c[,",1:ncol(c),"]", collapse="+"),
#                 ")"
#           )
#         )
#       )[,-1]
#  
#       # parameters for confounder exposure relation
#       parmsC_pi <- rep(log(1.5), c_number)
#  
#       # simulate the exposure
#       x <- exposure #rbinom(n, size = 1, expit(-.5 + piMat%*%parmsC_pi))
#  
#       # parameters for the confounder outcome relation
#       parmsC_mu <- rep(log(2), c_number)
#  
#       # simulate the outcome
#       pY <- mean(expit(intercept + log(2)*x + piMat%*%parmsC_mu))
#  
#       print(paste("the mean of y is: ", pY))
#  
#       ## what would we change here if interested in marginally adjusted risk difference or risk ratio?
#       res <- pY/(1 - pY)
#  
#       return(res)
#   }
#  
#   odds1 <- simulation_function_true(intercept = -2, exposure = 1, c_number = 10)
#   odds0 <- simulation_function_true(intercept = -2, exposure = 0, c_number = 10)
#  
#   or_marg1 <- log(odds1/odds0)
#  
#   or_marg1
#  
#   odds1 <- simulation_function_true(intercept = -2, exposure = 1, c_number = 25)
#   odds0 <- simulation_function_true(intercept = -2, exposure = 0, c_number = 25)
#  
#   or_marg2 <- log(odds1/odds0)
#  
#   or_marg2
#  
#   true_log_or1 <- or_marg1
#   true_log_or2 <- or_marg2
#  
#   # the true log OR marg1 (c_number = 10) is: 0.4132932
#   #
#   # the true log OR marg2 (c_number = 25) is: 0.289803

## ----tidy = F, warning = F, message = F, eval = F-----------------------------
#  
#   library(parallel)
#   library(lmtest)
#   library(sandwich)
#  
#   expit<-function(a){1/(1+exp(-a))}
#  
#   set.seed(123)
#  
#   simulation_function <- function(index,
#                                   intercept=-2,
#                                   sample_size = 1000,
#                                   c_number = 10,
#                                   cov_mat = 0,
#                                   diag_cov = 1){
#  
#       # printing to console won't work with parallel processing
#       print(index)
#  
#       # DATA GENERATION
#       n <- sample_size
#       print(n)
#  
#       # how many confounders to simulate?
#       p <- c_number
#       print(p)
#  
#       ## confounder matrix
#       sigma <- matrix(cov_mat, nrow=p, ncol=p)
#       diag(sigma) <- diag_cov
#       c <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sigma)
#  
#       # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
#       piMat <- model.matrix(
#         as.formula(
#           paste("~(",
#                 paste("c[,",1:ncol(c),"]", collapse="+"),
#                 ")"
#                 )
#           )
#         )[,-1]
#  
#       # parameters for confounder exposure relation
#       parmsC_pi <- rep(log(1.5), c_number)
#  
#       # simulate the exposure
#       x <- rbinom(n, size = 1, expit(-.5 + piMat%*%parmsC_pi))
#  
#       # parameters for the confounder outcome relation
#       parmsC_mu <- rep(log(2), c_number)
#  
#       # simulate the outcome
#       y <- rbinom(n, size = 1, expit(intercept + log(2)*x + piMat%*%parmsC_mu))
#  
#       # construct dataset
#       analysis_data <- data.frame(y, x, c)
#  
#       # ANALYSIS
#  
#       # marginal standardization
#       meanY <- mean(analysis_data$y)
#  
#       m1 <- glm(y ~ ., data = analysis_data, family=binomial(link="logit"))
#  
#       muhat1 <- mean(predict(m1, newdata = transform(analysis_data, x = 1), type="response"))
#       muhat0 <- mean(predict(m1, newdata = transform(analysis_data, x = 0), type="response"))
#  
#       ## compute the odds from these average probabilities
#       odds1<-muhat1/(1-muhat1)
#       odds0<-muhat0/(1-muhat0)
#  
#       ## marginal log-odds ratio
#       marginal_standardization_estimate <- log(odds1/odds0)
#  
#       # IP WEIGHTING
#       # create the propensity score in the dataset
#       analysis_data$propensity_score <- glm(x ~ .,
#                                             data = analysis_data[,-y],
#                                             family = binomial("logit"))$fitted.values
#  
#       # stabilized inverse probability weights
#       analysis_data$sw <- (mean(analysis_data$x)/analysis_data$propensity_score)*analysis_data$x +
#         ((1-mean(analysis_data$x))/(1-analysis_data$propensity_score))*(1-analysis_data$x)
#  
#       summary(analysis_data$sw)
#       quantile(analysis_data$sw, .995)
#  
#       # trim the weights
#       analysis_data$sw <- ifelse(analysis_data$sw>quantile(analysis_data$sw, .995),
#                                  quantile(analysis_data$sw, .995),
#                                  analysis_data$sw)
#  
#       # outcome model
#       m2 <- glm(y ~ x, data = analysis_data, weights = sw, family = quasibinomial(link="logit"))
#  
#       ip_weighting_estimate <- coeftest(m2,
#                                         vcov. = vcovHC(m2, type = "HC3"))[2,1:2]
#  
#       # bootstrap SEs
#       boot_func <- function(nboot, data){
#  
#         a_dat <- data
#  
#         index <- sample(1:nrow(a_dat), nrow(a_dat), replace = T)
#  
#         boot_dat <- a_dat[index, ]
#  
#         # MARGINAL STANDARDIZATION
#         m1_ <- glm(y ~ ., data = boot_dat, family=binomial(link="logit"))
#  
#         muhat1_ <- mean(predict(m1_, newdata = transform(boot_dat, x=1), type="response"))
#         muhat0_ <- mean(predict(m1_, newdata = transform(boot_dat, x=0), type="response"))
#  
#         ## compute the odds from these average probabilities
#         odds1_ <- muhat1_/(1-muhat1_)
#         odds0_ <- muhat0_/(1-muhat0_)
#  
#         marginal_stand_ <- log(odds1_/odds0_)
#  
#         # IP WEIGHTING
#         # create the propensity score in the dataset
#         boot_dat$propensity_score <- glm(x ~ .,
#                                          data = boot_dat[,-y],
#                                          family = binomial("logit"))$fitted.values
#  
#         # stabilized inverse probability weights
#         boot_dat$sw <- (mean(boot_dat$x)/boot_dat$propensity_score)*boot_dat$x +
#           ((1-mean(boot_dat$x))/(1-boot_dat$propensity_score))*(1-boot_dat$x)
#  
#         # trim the weights
#         boot_dat$sw <- ifelse(boot_dat$sw>quantile(boot_dat$sw, .995),
#                               quantile(boot_dat$sw, .995),
#                               boot_dat$sw)
#  
#         # outcome model
#         m2_ <- glm(y ~ x,
#                    data = boot_dat,
#                    weights = sw,
#                    family = quasibinomial(link="logit"))
#  
#         ip_weighting_ <- summary(m2_)$coefficients[2,1]
#  
#         return(c(marginal_stand_, ip_weighting_))
#       }
#  
#       analysis_data <- analysis_data %>% select(-propensity_score, -sw)
#  
#       boot_res <- lapply(1:500, function(x) boot_func(x, data = analysis_data))
#  
#       boot_SE <- apply(do.call(rbind, boot_res), 2, sd)
#  
#       # SIMULATION FUNCTION OUTPUT
#  
#       res <- data.frame(intercept = intercept,
#                         sample_size = sample_size,
#                         c_number = c_number,
#                         cov_mat  = cov_mat,
#                         diag_cov = diag_cov,
#                         marginal_standardization_estimate,
#                         marginal_standardization_SE = boot_SE[1],
#                         ip_weighting_estimate = ip_weighting_estimate[1],
#                         ip_weighting_robust_SE = ip_weighting_estimate[2],
#                         ip_weighting_boot_SE = boot_SE[2])
#  
#       return(res)
#   }
#  
#   # simulation function parameters
#  
#   parm_data <- expand.grid(
#     index = 1:200,
#     sample_size = 1000,
#     intercept = -2,
#     c_number = c(10, 25),
#     cov_mat = 0,
#     diag_cov = 1
#   )
#  
#   head(parm_data, 3)
#  
#   tail(parm_data, 3)
#  
#   simulation_results <- mclapply(1:nrow(parm_data),
#                                function(x) simulation_function(index = parm_data[x,]$index,
#                                                                intercept = parm_data[x,]$intercept,
#                                                                sample_size = parm_data[x,]$sample_size,
#                                                                c_number = parm_data[x,]$c_number,
#                                                                cov_mat = parm_data[x,]$cov_mat,
#                                                                diag_cov = parm_data[x,]$diag_cov),
#                                mc.cores = detectCores() - 2)
#  
#   sim_res <- do.call(rbind, simulation_results)
#  
#   # # save the data to file!
#  
#   head(sim_res)
#  
#   write_csv(sim_res, here("lectures/06_Section6_SimulationInPractice", "simulation_results_MSIPW.csv"))
#  

## ----tidy = F, warning = F, message = F---------------------------------------

a <- read_csv(here("lectures/06_Section6_SimulationInPractice", "simulation_results_MSIPW.csv"))

head(a)

# to make things easier, let's create objects with our true values:

true_log_or1 <- 0.4132932
true_log_or2 <- 0.289803


## ----tidy = F, warning = F, message = F---------------------------------------

# the true log OR marg1 (c_number = 10) is: 0.4132932
#
# the true log OR marg2 (c_number = 25) is: 0.289803

a %>% 
  group_by(c_number) %>% 
  summarize(meanMS = mean(marginal_standardization_estimate),
            meanIPW = mean(ip_weighting_estimate))


## ----tidy = F, warning = F, message = F---------------------------------------

# the true log OR marg1 (c_number = 10) is: 0.4132932
#
# the true log OR marg2 (c_number = 25) is: 0.289803

a %>% 
  group_by(c_number) %>% 
  summarize(biasMS = mean(marginal_standardization_estimate - true_log_or1),
            biasIPW = mean(ip_weighting_estimate - true_log_or2))


## ----tidy = F, warning = F, message = F---------------------------------------
mc_se_bias <- function(x, n){
  sqrt(sum((x - mean(x))^2)/(n*(n-1)))
}

bias_results <- a %>% 
  group_by(c_number) %>% 
  summarize(biasMS = mean(marginal_standardization_estimate - true_log_or1),
            biasIPW = mean(ip_weighting_estimate - true_log_or2),
            
            biasMS_se = mc_se_bias(marginal_standardization_estimate - true_log_or1, n = 200),
            biasIPW_se = mean(ip_weighting_estimate - true_log_or2, n = 200),
            
            biasMS_p.value = round(2*(1 - pnorm(abs(biasMS/biasMS_se))),4),
            biasIPW_p.value = round(2*(1 - pnorm(abs(biasIPW/biasIPW_se))),4))

bias_results


## ----tidy = F, warning = F, message = F---------------------------------------

ggplot(a) +
  geom_histogram(aes(marginal_standardization_estimate - true_log_or1)) +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(~c_number)


## ----tidy = F, warning = F, message = F---------------------------------------

ggplot(a) +
  geom_histogram(aes(ip_weighting_estimate - true_log_or2)) +
  geom_vline(xintercept = 0, color = "red") +
  facet_wrap(~c_number)


