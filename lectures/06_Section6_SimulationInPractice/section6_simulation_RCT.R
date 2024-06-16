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


## ----out.width = "10cm",fig.cap="Simple causal diagram representing a randomized trial, with randomized treatmetn assignment $A$, outcome cause $C$, and outcome $Y$.",echo=F----
knitr::include_graphics(here("_images","rct_dag.pdf"))

## ----tidy = F, warning = F, message = F, eval = F-----------------------------
#  library(parallel)
#  library(lmtest)
#  library(sandwich)
#  
#  expit<-function(a){1/(1+exp(-a))}
#  
#  set.seed(123)
#  
#  simulation_function <- function(index,
#                                  sample_size = 500,
#                                  true_value = 2,
#                                  c_number = 10,
#                                  cov_mat = 0,
#                                  diag_cov = 1){
#  
#    # printing to console won't work with parallel processing
#    print(index)
#  
#    # DATA GENERATION
#    n <- sample_size
#    print(n)
#  
#    # how many confounders to simulate?
#    p <- c_number
#    print(p)
#  
#    ## confounder matrix
#    sigma <- matrix(cov_mat, nrow=p, ncol=p)
#    diag(sigma) <- diag_cov
#    c <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sigma)
#  
#    # DESIGN MATRIX FOR THE PROPENSITY SCORE MODEL
#    piMat <- model.matrix(
#      as.formula(
#        paste("~(",
#              paste("c[,",1:ncol(c),"]", collapse="+"),
#              ")"
#        )
#      )
#    )[,-1]
#  
#    # simulate the treatment via 50:50 randomization
#    x <- rbinom(n, size = 1, p = .5)
#  
#    # parameters for the covariate outcome relation
#    parmsC_mu <- rep(1.25, c_number)
#  
#    # simulate the outcome
#    y <- rnorm(n, mean = 100 + true_value*x + piMat%*%parmsC_mu, sd = 2)
#  
#    # construct dataset
#    analysis_data <- data.frame(y, x, c)
#  
#    # ANALYSIS
#  
#    mod1 <- lm(y ~ x, data = analysis_data)
#    unadjusted_effect <- summary(mod1)$coefficients[2,1:2]
#  
#    mod2 <- lm(y ~ ., data = analysis_data)
#    adjusted_effect <- summary(mod2)$coefficients[2,1:2]
#  
#    # SIMULATION FUNCTION OUTPUT
#  
#    res <- data.frame(sample_size = sample_size,
#                      c_number = c_number,
#                      cov_mat  = cov_mat,
#                      diag_cov = diag_cov,
#                      true_value = true_value,
#                      unadjusted_effect,
#                      adjusted_effect)
#  
#    return(res)
#  }
#  
#  simulation_results <- mclapply(1:5000,
#                                 function(x) simulation_function(index = x,
#                                                                 sample_size = 500,
#                                                                 true_value = 2,
#                                                                 c_number = 10,
#                                                                 cov_mat = 0,
#                                                                 diag_cov = 1),
#                                 mc.cores = detectCores() - 2)
#  
#  sim_res <- do.call(rbind, simulation_results)
#  
#  # # save the data to file!
#  
#  write_csv(sim_res, here("lectures/06_Section6_SimulationInPractice", "simulation_results_rct.csv"))
#  

## ----tidy = F, warning = F, message = F---------------------------------------

a <- read_csv(here("lectures/06_Section6_SimulationInPractice", "simulation_results_rct.csv"))

head(a)


## ----tidy = F, warning = F, message = F---------------------------------------

a %>% 
  group_by(sample_size, c_number, cov_mat) %>% 
  summarize(meanAdj = mean(adjusted_estimate),
            meanUnAdj = mean(unadjusted_estimate))


## ----tidy = F, warning = F, message = F---------------------------------------

a %>% 
  group_by(sample_size, c_number, cov_mat) %>% 
  summarize(meanAdjSE = mean(adjusted_estimate_se),
            meanUnAdjSE = mean(unadjusted_estimate_se))


## ----tidy = F, warning = F, message = F---------------------------------------

plot_res <- a %>% 
  group_by(sample_size, c_number, cov_mat) %>% 
  summarize(meanAdjSE = mean(adjusted_estimate_se),
            meanUnAdjSE = mean(unadjusted_estimate_se))

pacman::p_load(looplot)

p = nested_loop_plot(resdf = plot_res, 
                     x = "sample_size", steps = "c_number", grid_rows = "cov_mat",
                     steps_y_base = -.25, steps_y_height = .25, steps_y_shift = .1,
                     x_name = "Sample Size", y_name = "Standard Error",
                     spu_x_shift = 200,
                     steps_values_annotate = TRUE, steps_annotation_size = 3, 
                     hline_intercept = 0, 
                     y_expand_add = c(1, NULL), 
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 8)
                       ))
)

ggsave(here("_images", "nested_loop_plot_rct.pdf"), width = 8, height = 6)


## ----nested_loop_rct, out.width="10cm", fig.align='center', echo=F------------
knitr::include_graphics(here("_images", "nested_loop_plot_rct.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------
a %>% 
  group_by(sample_size, c_number, cov_mat) %>% 
  summarize(meanAdjSE = mean(adjusted_estimate_se),
            sdAdj = sd(adjusted_estimate)) %>% 
  mutate(se_bias = sdAdj - meanAdjSE)

a %>% 
  group_by(sample_size, c_number, cov_mat) %>% 
  summarize(meanUnadjSE = mean(unadjusted_estimate_se),
            sdUnadj = sd(unadjusted_estimate)) %>% 
  mutate(se_bias = sdUnadj - meanUnadjSE)

