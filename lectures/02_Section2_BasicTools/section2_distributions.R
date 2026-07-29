## ----setup, include=FALSE-----------------------------------------------------

# Load pacman (install if necessary)
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = 'https://cloud.r-project.org')
}

# Use pacman to load (and install if missing) all required packages
pacman::p_load(
  data.table, tidyverse, ggplot2, ggExtra, formatR,
  gridExtra, skimr, here, RColorBrewer,
  knitr, remotes
)

# Ensure fontawesome is installed from GitHub if not available
if (!requireNamespace("fontawesome", quietly = TRUE)) {
  remotes::install_github("rstudio/fontawesome")
}
library(fontawesome)

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


## ----gaussianplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for Univariate Normal Distribution with Mean = 0 and Standard Deviation = 1 for 5000 Simulated Observations."----

set.seed(123)
ggplot() + 
  geom_histogram(aes(x = rnorm(n = 5000))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----warning = F, message = F-------------------------------------------------

set.seed(123)

n <- 5

y <- rnorm(n, mean = 0, sd = 1)

y


## ----warning = F, message = F-------------------------------------------------

n <- 5

set.seed(123)
x <- rnorm(n)

set.seed(123)
y_version1 <- rnorm(n, mean = 1 + 2*x, sd = 1)

set.seed(123)
y_version2 <- 1 + 2*x + rnorm(n, mean = 0, sd = 1)

y_version1

y_version2


## ----mvnplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Contour Plot for Multivariate Normal Distribution with Mean = (0,0), Standard Deviation = (1,1), and covariance (.5, .5) for 5000 Simulated Observations."----

set.seed(123)
m <- c(0, 0)
sigma <- matrix(c(1,.5,.5,1), nrow=2)
data.grid <- expand.grid(x = seq(-3, 3, length.out=200), y = seq(-3, 3, length.out=200))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))
ggplot(q.samp, aes(x=x, y=y, z=prob)) + 
    geom_contour() +
    coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3), ratio = 1) 


## ----warning = F, message = F-------------------------------------------------

  n <- 5

  set.seed(123)

  # create variance - covariance matrix:
  sigma <- matrix(0,nrow=3,ncol=3) 
  diag(sigma) <- 1
  
  sigma
  
  # create mean vector:
  mu <- rep(1, 3)
  
  mu
   
  # simulate variables
  c <- MASS::mvrnorm(n, mu=mu, Sigma=sigma)
  
  c


## ----warning = F, message = F-------------------------------------------------

  n <- 5

  set.seed(123)

  # create variance - covariance matrix:
  sigma <- matrix(0,nrow=3,ncol=3) 
  diag(sigma) <- 1
  
  sigma
  
  # create mean vector:
  mu <- rep(1, 3)
  
  mu
   
  # simulate variables
  c <- mvtnorm::rmvnorm(n, mean=mu, sigma=sigma)
  
  c


## ----tidy = F, warning = F, message = F---------------------------------------
  
library(mvtnorm)

n <- 5
p <- c_number <- 3

# Covariate matrix
sigma <- diag(p)
c <- rmvnorm(n, mean = rep(0, p), sigma = sigma)

# Outcome model
parmsC <- rep(1.5, c_number)
y <- 10 + c %*% parmsC + rnorm(n)

# Final data
data <- data.frame(y, c)


## ----uniformplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Uniform Distribution with Upper and Lower Bounds of 0 and 1, respectively for 5000 Simulated Observations."----

set.seed(123)
ggplot() + 
  geom_histogram(aes(x = runif(n = 5000))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## -----------------------------------------------------------------------------

set.seed(123)

n <- 5

y <- runif(n, min = 0, max = 1)

y


## ----binomplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Barplot for the Binomial (Bernoulli) Distribution with p = 0.25 for 5000 Simulated Observations."----

set.seed(123)
ggplot() + 
  geom_bar(aes(x = rbinom(n = 5000, size = 1, p = .25))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## -----------------------------------------------------------------------------

set.seed(123)

n <- 5

y <- rbinom(n, size = 1, p = .5)

y


## -----------------------------------------------------------------------------

set.seed(123)

n <- 5

y <- rbinom(n, size = 8, p = .5)

y


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

n <- 5

c <- rnorm(n)

expit <- function(x){1/(1+exp(-x))}

y <- rbinom(n, size = 1, p = expit(-2 + log(2)*c))

data <- data.frame(y, c)

data


## ----multinomplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Barplot for the Multinomial Distribution with Three Levels and p = {0.2, 0.1, 0.7} for 5000 Simulated Observations."----

set.seed(123)

mn_vars <- t(rmultinom(n = 5000, size = 1, p = c(.2, .1, .7)))
mn_vars <- do.call(rbind,
  lapply(1:nrow(mn_vars), function(x) which(mn_vars[x,]==1))
)

ggplot() +
  geom_bar(aes(x = mn_vars)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## -----------------------------------------------------------------------------

set.seed(123)

n <- 5

y <- rmultinom(n, size = 1, p = rep(1/6, 6))

y


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

n <- 5

y <- rmultinom(n, size = 1, 
               p = c(0.1,0.05,0.15, 1 - sum(0.1,0.05,0.15)))

y

t(y)


## ----poisplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Poisson Distribution with lamba = 5 for 5000 Simulated Observations."----

set.seed(123)
ggplot() + 
  geom_bar(aes(x = rpois(n = 5000, lambda = 5))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## -----------------------------------------------------------------------------

set.seed(123)

n <- 5

y <- rpois(n, lambda = 3)

y


## ----nbinomplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Negative Binomial Distribution with p = 0.5 and size = 100 for 5000 Simulated Observations."----

set.seed(123)
ggplot() +
  geom_bar(aes(x = rnbinom(n = 5000, prob = .2, size = 10) )) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

y <- rnbinom(n = 5, prob = .2, size = 2) 

y


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

y <- rnbinom(n = 5, mu = 5, size = 2) 

y


## ----geomplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Geometric Distribution with p = 0.5 for 5000 Simulated Observations."----

set.seed(123)
ggplot() +
  geom_bar(aes(x = rgeom(n = 5000, prob = .5) )) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

y <- rgeom(n = 5, prob = .5)

y


## ----exponentialplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Exponential Distribution with rate = 1 for 5000 Simulated Observations."----

set.seed(123)
ggplot() + 
  geom_histogram(aes(x = rexp(n = 5000, rate = 1)), bins = 50) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

n <- 5

y <- rexp(n, rate = 0.5)

y


## ----tidy = F, warning = F, message = F---------------------------------------

sim_dat <- data.frame(y_obs = ifelse(y<=1, y, 1), 
                      delta = as.numeric(y <= 1))

sim_dat


## ----weibullplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Weibull Distribution with a shape parameter = 2 and a scale parameter = 1 for 5000 Simulated Observations."----

set.seed(123)
ggplot() + 
  geom_histogram(aes(x = rweibull(n = 5000, shape = 2, scale = 1)), bins = 50) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)
n <- 5
shape <- 2
scale <- 1
y <- rweibull(n, shape = shape, scale = scale)
y


## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

n <- 5000
rate <- 2

y1 <- rweibull(n, shape = 1, scale = 1 / rate)

y2 <- rexp(n, rate = rate)


## ----tidy = F, warning = F, message = F, echo = F-----------------------------
df <- data.frame(
  value = c(y2, y1),
  distribution = factor(rep(c("Exponential", "Weibull"), each = n))
)

# Create histogram bins for each distribution
bin_width <- 0.2
max_val <- max(df$value)
breaks <- seq(0, ceiling(max_val), by = bin_width)

hist_df <- df %>%
  mutate(bin = cut(value, breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
  group_by(distribution, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    count = ifelse(distribution == "Weibull", -count, count),
    bin_mid = as.numeric(sub("\\[(.+),.*", "\\1", bin)) + bin_width / 2
  )

# Plot: Vertical bars with central x-axis
p1 <- ggplot(hist_df, aes(x = bin_mid, y = count, fill = distribution)) +
  geom_col(width = bin_width, color = "black", alpha = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(labels = abs) +
  labs(
    x = "Random Variable Value",
    y = "Count"
  ) +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.y = element_text(color = "black"),
    legend.title=element_blank()
  ) +
  scale_fill_manual(values = c("Exponential" = "steelblue", "Weibull" = "darkorange"))

ggsave(here("_images", "weibull_exponential_comparison.pdf"), width = 8, height = 8, units = "cm")


## ----expweib, out.width="8cm", fig.align='center', fig.cap="Comparison of the Distribution of a Random Variable Generated from the rexp() function with a rate parameter = 2, and the rweibull() function with shape = 1 and scale = 1 / rate.", echo=F----
knitr::include_graphics(here("_images", "weibull_exponential_comparison.pdf"))

## ----gengammaplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Histogram for the Generalized Gamma Distribution with a location parameter = 0, a sigma parameter = 0.5, and a Q parameter = 0.5 for 5000 Simulated Observations."----

pacman::p_load(flexsurv)
shape <- 2
scale <- 1
k <- 0.5
set.seed(123)
ggplot() + 
  geom_histogram(aes(x = rgengamma(n, mu = log(scale), sigma = 1/shape, Q = k)), bins = 50) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


## ----tidy = F, warning = F, message = F---------------------------------------

pacman::p_load(flexsurv)

set.seed(123)
n <- 5
shape <- 2
scale <- 1
q <- 0.5

y <- rgengamma(n, mu = log(scale), sigma = 1/shape, Q = q)

y


## -----------------------------------------------------------------------------

set.seed(123)
n <- 5000
rate <- 2
scale <- 1 / rate

# Generate generalized gamma with exponential parameters
y <- rgengamma(n, mu = log(scale), sigma = 1, Q = 1)

# Compare with exponential
exp_y <- rexp(n, rate = rate)

# Kolmogorov-Smirnov test
ks.test(y, exp_y)



## ----tidy = F, warning = F, message = F---------------------------------------

set.seed(123)

pareto_scale <- 10000 # this sets the minimum to the distribution
pareto_shape <- log(5, base = 4) # this sets the distribution to be Pareto 80:20

# now we use the inverse transform sampling method 
## generate uniform random variable
u <- runif(n = 1000)

## the equation that let's us apply inverse transform sampling from a uniform
## distribution to get a Pareto random variable
x <- pareto_scale*(1-u)^(-1/pareto_shape)

## Draw RV from same distribution but with the EnvStats Package function
## to compare the two random variables
pacman::p_load(EnvStats)
x_prime <- rpareto(n = 1000, 
                   location = pareto_scale, 
                   shape = pareto_shape)


## ----paretocomparison, out.width="5cm", fig.align='center', fig.margin=TRUE, fig.cap="Log of the complement of the CDF against the log of X for the random variable generated using the inverse transform method (red lines) and the rpareto function in the EnvStats package (blue lines).", echo=F----
knitr::include_graphics(here("_images", "pareto_comparison.pdf"))

