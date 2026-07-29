## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
library(formatR)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR",
               "gridExtra","skimr","here","RColorBrewer")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='https://cloud.r-project.org')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

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

# estimating pi
pi_est <- function(simulation_n) {

   x <- runif(simulation_n, 0, 1)
   y <- runif(simulation_n, 0, 1)

   # if point lies within radius, set to one, otherwise zero
   rad <- as.numeric(x^2 + y^2 <= 1)
   
   # proportion of points within quarter circle (circle area) to 
   # proportion of points within unit square (square area)
   res <- c(simulation_n, (sum(rad)/simulation_n)*4)
   
   return(res)

}

# set the seed value
set.seed(41)

n_list <- c(500, 1000, 5000, 50000, 100000, 500000, 10000000)

pi_data <- NULL
for(i in n_list){
  pi_data <- rbind(pi_data, pi_est(simulation_n = i))
}


## ----tidy = F, message = F, warning = F---------------------------------------
knitr::kable(pi_data, "simple",
             col.names = c("Simulation N", 
                           "Estimated pi"))


## ----tidy = F, warning = F, message = F---------------------------------------

# write a function that computes the statistical 
# mode. modified from: https://stackoverflow.com/a/8189441
mode_func <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  if(max(tab)==1){
    warning("no duplicates in data, mode does not exist")
    res <- sample(ux, size = 1)
  } else{
    res <- ux[which.max(tab)]
  }
  return(res)
}

# set the seed value
set.seed(123)

# how many observations?
n = 200

# start the simulation loop
res <- NULL
for(i in 1:1e4){
  
  # simulate a variable from a normal 
  # with mean 26.5 and SD 5.25
  y <- round(rnorm(n, mean = 26.5, sd = 5.25))
  
  # estimate mean, median, and mode
  mode_estimator <- mode_func(y)
  mean_estimator <- mean(y)
  median_estimator <- median(y)
  
  # store results
  res <- rbind(res, 
               cbind(mode_estimator, 
                     mean_estimator, 
                     median_estimator)
  )
  
}

head(res)

# convert results to data frame
res <- data.frame(res)

# plot results: one histogram per estimator, overlaid,
# with a colorblind-friendly palette
res %>%
  pivot_longer(everything(),
               names_to = "estimator",
               values_to = "estimate") %>%
  mutate(estimator = str_remove(estimator, "_estimator")) %>%
  ggplot() +
  geom_histogram(aes(estimate, fill = estimator),
                 alpha = .5, position = "identity") +
  scale_fill_manual(values = c(mean = "#56B4E9",
                               median = "#E69F00",
                               mode = "#009E73")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

res %>% 
  summarize(across(everything(),
                   list(mean = mean, 
                        std_dev = sd)))


