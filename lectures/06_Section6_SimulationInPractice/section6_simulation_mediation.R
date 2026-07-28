## ----setup, include=FALSE-----------------------------------------------------
# Load pacman (install if necessary)
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

# Use pacman to load (and install if missing) all required packages
pacman::p_load(
  tidyverse, ggplot2, formatR, gridExtra, here,
  lmtest, sandwich, knitr
)

# Set ggplot2 theme
thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

# Set global chunk options
opts_chunk$set(
  tidy.opts = list(width.cutoff = 60),
  tidy = TRUE,
  echo = TRUE
)

# Set knit hook for purling
knit_hooks$set(purl = hook_purl)

## ----out.width = "10cm",fig.cap="Causal diagram representing exposure induced mediator outcome confounding.",echo=F----
knitr::include_graphics(here("_images","mediation_dag.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

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


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

true_CDE <- function(sample_size, exposure, mediator){

  n <- sample_size

  U <- rnorm(n)

  C <- rnorm(n)

  X <- exposure # set the exposure to a specific value

  L <- rnorm(n, mean = 0 + 1.5*X + 1.5*U)

  Z <- mediator # set the mediator to a specific value

  return(mean(10 + 5*X + 5*C + 5*Z + 4*L + 4*U))

}

# CDE with the mediator fixed at z = 0
CDE_0 <- true_CDE(sample_size = 5e6, exposure = 1, mediator = 0) -
  true_CDE(sample_size = 5e6, exposure = 0, mediator = 0)

CDE_0

# CDE with the mediator fixed at z = 1
CDE_1 <- true_CDE(sample_size = 5e6, exposure = 1, mediator = 1) -
  true_CDE(sample_size = 5e6, exposure = 0, mediator = 1)

CDE_1


## -----------------------------------------------------------------------------

true_cde <- 11


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

simulation_function <- function(index, sample_size = 500){

  # DATA GENERATION
  n <- sample_size

  U <- rnorm(n)
  C <- rnorm(n)
  X <- rbinom(n, size = 1, p = expit(-1 + log(2)*C ))
  L <- rnorm(n, mean = 0 + 1.5*X + 1.5*U)
  Z <- rbinom(n, size = 1, expit(-1.5 + log(2)*X + log(2)*L ))
  Y <- rnorm(n, mean = 10 + 5*X + 5*C + 5*Z + 4*L + 4*U, sd = 5)

  med_data <- data.frame(Y, Z, L, X, C)

  # ANALYSIS 1: standard regression, adjusting for L
  model1 <- glm(Y ~ X + Z + C + L, data = med_data,
                family = gaussian("identity"))
  adjL <- summary(model1)$coefficients["X", 1:2]

  # ANALYSIS 2: standard regression, ignoring L
  model2 <- glm(Y ~ X + Z + C, data = med_data,
                family = gaussian("identity"))
  noL <- summary(model2)$coefficients["X", 1:2]

  # ANALYSIS 3: IP weighted marginal structural model
  ## stabilized weights for the exposure
  psX <- glm(X ~ C, data = med_data,
             family = binomial("logit"))$fitted.values
  swX <- (mean(med_data$X)/psX)*med_data$X +
    ((1 - mean(med_data$X))/(1 - psX))*(1 - med_data$X)

  ## stabilized weights for the mediator
  psZ <- glm(Z ~ X + C + L, data = med_data,
             family = binomial("logit"))$fitted.values
  swZ <- (mean(med_data$Z)/psZ)*med_data$Z +
    ((1 - mean(med_data$Z))/(1 - psZ))*(1 - med_data$Z)

  med_data$sw <- swX*swZ

  ## weighted (saturated) marginal structural model
  model3 <- glm(Y ~ X + Z + X:Z, data = med_data, weights = sw,
                family = gaussian("identity"))
  ipw <- coeftest(model3,
                  vcov. = vcovHC(model3, type = "HC3"))["X", 1:2]

  # SIMULATION FUNCTION OUTPUT: one row per estimator
  res <- data.frame(
    index = index,
    sample_size = sample_size,
    estimator = c("Regression with L",
                  "Regression without L",
                  "IP Weighted MSM"),
    estimate = c(adjL[1], noL[1], ipw[1]),
    se = c(adjL[2], noL[2], ipw[2]),
    row.names = NULL
  )

  return(res)
}


## ----tidy = F, warning = F, message = F---------------------------------------

nsim <- 2000

sim_res <- lapply(1:nsim, function(x) simulation_function(index = x))

sim_res <- do.call(rbind, sim_res)

head(sim_res)


## ----tidy = F, warning = F, message = F---------------------------------------

mc_se_bias <- function(x){
  n <- length(x)
  sqrt(sum((x - mean(x))^2)/(n*(n-1)))
}

perf <- sim_res %>%
  mutate(lcl = estimate - 1.96*se,
         ucl = estimate + 1.96*se) %>%
  group_by(estimator) %>%
  summarize(mean_estimate = mean(estimate),
            bias = mean(estimate - true_cde),
            bias_mcse = mc_se_bias(estimate - true_cde),
            emp_sd = sd(estimate),
            mean_se = mean(se),
            coverage = mean(lcl < true_cde & true_cde < ucl)) %>%
  arrange(desc(abs(bias)))

knitr::kable(perf, digits = 2,
             col.names = c("Estimator", "Mean Est.", "Bias",
                           "Bias MCSE", "Emp. SD", "Mean SE",
                           "Coverage"))


## ----medhist, tidy = F, warning = F, message = F, out.width="10cm", fig.align='center', fig.cap="Distribution of controlled direct effect estimates across all Monte Carlo runs for each of the three estimators. The dashed red line marks the true CDE of 11."----

ggplot(sim_res) +
  geom_histogram(aes(estimate), bins = 60) +
  geom_vline(xintercept = true_cde, color = "red", linetype = 2) +
  facet_wrap(~estimator, ncol = 1) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Estimated CDE") + ylab("Count")


## ----include = F--------------------------------------------------------------
# convenience objects for the inline numbers quoted below
p_adjL <- perf %>% filter(estimator == "Regression with L")
p_noL  <- perf %>% filter(estimator == "Regression without L")
p_ipw  <- perf %>% filter(estimator == "IP Weighted MSM")

