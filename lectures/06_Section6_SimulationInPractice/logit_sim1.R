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

## simulate data from a logistic model with four variables:
#### one outcome
#### one exposure
#### two confounders
#### 
#### suppose there's no interactions in the true DGM (between 2 confounders)
#### 
#### what happens when we fit a model with versus without interactions?


expit <- function(x){1/(1+exp(-x))} # or use plogis()

set.seed(123)

simulation_function <- function(index, sample_size = 50){
  
  ## DATA GENERATING MECHANSIM CODE
  #index = 1
  print(index)
  
  n <- sample_size
  
  print("confounders")
  c1 <- rbinom(n, size = 1, prob = .5)
  
  c2 <- rnorm(n, mean = 0, sd = 1)
  
  print("exposure")
  x <- rbinom(n, size = 1, prob = expit(-1 + log(2)*c1 + log(1.25)*c2))
  
  print("outcome")
  y <- rbinom(n, size = 1, prob = expit(-2 + log(1.5)*x + log(2)*c1 + log(1.25)*c2))
  
  a <- data.frame(y, x, c1, c2)  
  
  ## ANALYSIS CODE 
  ## 
  print("model1")
  mod1 <- glm(y ~ x + c1 + c2, data = a, family = binomial(link = "logit"))
  
  mod1_fitted.values <- mod1$fitted.values
  
  print("model2")
  mod2 <- glm(y ~ x + c1 + c2 + c1:c2, data = a, family = binomial(link = "logit"))
  
  mod2_fitted.values <- mod2$fitted.values
  
  summary(mod2)
  
  ## SIMULATION RESULTS
  ##
  print("extract results")
  res_mod1 <- summary(mod1)$coefficients["x",1]
  
  res_mod2 <- summary(mod2)$coefficients["x",1]
  
  res1 <- data.frame(index, sample_size, res_mod1, res_mod2)
  
  res2 <- data.frame(index, mod1_fitted.values, mod2_fitted.values)

  
  if(index == 1){
    write_csv(x = res1, file = here("data", "parameter_results.csv"))
  } else{
    write_csv(x = res1, file = here("data", "parameter_results.csv"), append = T)
  }
  
  if(index ==1){
    write_csv(x = res2, file = here("data", "fitted_results.csv"))
  } else{
    write_csv(x = res2, file = here("data", "fitted_results.csv"), append = T)
  }
  
  
  # print("create output data")
  # res <- list(
  #   data.frame(index, sample_size, res_mod1, res_mod2), 
  #   data.frame(mod1_fitted.values, mod2_fitted.values)
  # )
  # 
  # return(res)
  
}

#simulation_function(index = 1, sample_size = 5000)

parms_data <- expand.grid(
  index = 1:1000,
  sample_size = c(100, 1000, 5000)
)

parms_data %>% 
  group_by(sample_size) %>% 
  count()

sim_res <- lapply(1:nrow(parms_data), 
                  function(x) simulation_function(index = parms_data[x,]$index, 
                                                  sample_size = parms_data[x,]$sample_size))


newParms_data <- read_csv(here("data", "parameter_results.csv"))

head(newParms_data)

parameter_results <- NULL
for(i in 1:1000){
  parameter_results <- rbind(parameter_results, 
                             sim_res[[i]][[1]])
}

head(parameter_results)

fitted_results <- NULL
for(i in 1:1000){
  fitted_results <- rbind(fitted_results, 
                          cbind(i, sim_res[[i]][[2]]))
}

head(fitted_results)
tail(fitted_results)

sim_res <- do.call(rbind, sim_res)

dim(sim_res)

head(sim_res)

## simulation analysis 

ggplot(sim_res) +
  geom_histogram(aes(x = noInt)) +
  geom_vline(xintercept = log(1.5), color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~sample_size)

ggplot(sim_res) +
  geom_histogram(aes(x = Int)) +
  geom_vline(xintercept = log(1.5), color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~sample_size)

## bias 

# sim_res %>% 
#   group_by(sample_size) %>% 
#   summarize(... ... ... )

mean(sim_res$noInt - log(1.5))
mean(sim_res$Int - log(1.5))

## mean squared error 

mean((sim_res$noInt - log(1.5))^2)
mean((sim_res$Int - log(1.5))^2)

## are the two estimators different
t.test(sim_res$noInt, sim_res$Int, paired = T)

