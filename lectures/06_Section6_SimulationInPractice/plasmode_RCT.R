pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom,
  AIPW
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

?AIPW

?eager_sim_obs

data("eager_sim_obs")

a <- eager_sim_obs

head(eager_sim_obs)

mod_eager <- glm(sim_Y ~ sim_A + eligibility + loss_num + age + I(age^2), 
                 data = eager_sim_obs, 
                 family = binomial("logit"))

parms_obs <- summary(mod_eager)$coefficients[,1]

dim(a)

head(a)

a <- a %>% select(-sim_A, -sim_Y)

head(a)

### simulation code

expit<-function(a){1/(1+exp(-a))}

set.seed(123)

simulation_function <- function(index, 
                                sample_size = 500,
                                true_value = 2,
                                c_number = 10,
                                cov_mat = 0,
                                diag_cov = 1,
                                outcome_coef = parms_obs,
                                confounder_coef = parms_obs[4]){
  
  # printing to console won't work with parallel processing
  print(index)
  
  # DATA GENERATION
  n <- sample_size 
  print(n)
  
  # how many covariates to simulate?
  p <- c_number
  print(p)
  
  a_new <- a %>% slice_sample(n = sample_size, 
                              replace = T) %>% 
    select(eligibility, loss_num, age) %>% 
    mutate(loss_num = as.numeric(loss_num==2))
  
  a_new

  # simulate the treatment via 50:50 randomization
  x <- rbinom(n, size = 1, p = .5)
  
  # parameters for the covariate outcome relation
  parmsC_mu <- rep(1.25, c_number)
  
  # simulate the outcome
  
  mu <- outcome_coef[1] + true_value*x - outcome_coef[3]*a_new$eligibility + confounder_coef*a_new$loss_num + 
    outcome_coef[5]*scale(a_new$age) + outcome_coef[6]*I(scale(a_new$age)^2)
  
  y <- rnorm(n, mean = mu, sd = 2)
  
  # construct dataset
  analysis_data <- data.frame(y, x, a_new)
  
  head(analysis_data)
  
  # ANALYSIS
  
  mod1 <- lm(y ~ x, data = analysis_data)
  unadjusted_effect <- summary(mod1)$coefficients["x",1:2]
  
  mod2 <- lm(y ~ ., data = analysis_data)
  adjusted_effect <- summary(mod2)$coefficients["x",1:2]
  
  # SIMULATION FUNCTION OUTPUT
  
  res <- data.frame(index = index, 
                    sample_size = sample_size,
                    c_number = c_number,
                    cov_mat  = cov_mat,
                    diag_cov = diag_cov,
                    true_value = true_value,
                    unadjusted_effect = unadjusted_effect[1],
                    unadjusted_effect_se = unadjusted_effect[2],
                    adjusted_effect = adjusted_effect[1],
                    adjusted_effect_se = adjusted_effect[2])
  
  return(res)
}

simulation_function(index = 1, 
                    sample_size = 500,
                    true_value = 2,
                    c_number = 10,
                    cov_mat = 0,
                    diag_cov = 1,
                    outcome_coef = parms_obs)

library(parallel)

detectCores()

simulation_results <- mclapply(1:500,
                               function(x) simulation_function(index = x, 
                                                               sample_size = 500,
                                                               true_value = 2,
                                                               c_number = 10,
                                                               cov_mat = 0,
                                                               diag_cov = 1),
                               mc.cores = detectCores() - 2)

sim_res <- do.call(rbind, simulation_results)

dim(sim_res)

write_csv(sim_res, here("data", "RCT_sim_results.csv"))

### 

