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


## ----mvnplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Contour Plot for Multivariate Normal Distribution with Mean = [0,0], and Standard Deviation = [1,1], and covariance [.5, .5] for 5000 Simulated Observations."----

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


## ----tidy = FALSE, attr.source='.numberLines'---------------------------------
  
  n = 5
  
  p <- c_number <- 3
  
  ## confounder matrix
  sigma <- matrix(0,nrow=p,ncol=p)
  diag(sigma) <- 1
  c <- mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sigma)
  
  # DESIGN MATRIX FOR THE OUTCOME MODEL
  muMatT <- model.matrix(  
    as.formula(  
      paste("~(",  
            paste("c[,",1:ncol(c),"]", collapse="+"),  
            ")"  
            )  
      )  
    )[,-1]
  
  parmsC <- rep(1.5,c_number)
  
  y <- 10 + muMatT%*%parmsC + rnorm(n)
  
  data.frame(y,c)

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


## ----multinomplot, out.width="5cm", fig.align='center', fig.margin=TRUE, warning = F, message = F, echo=F, fig.cap="Barplot for the Multinomial Distribution with Three Levels and p = {0.2, 0.1, 0.7} for 5000 Simulated Observations."----

mn_vars <- t(rmultinom(n = 5000, size = 1, p = c(.2, .1, .7)))
mn_vars <- do.call(rbind, 
  lapply(1:nrow(mn_vars), function(x) which(mn_vars[x,]==1))
)

set.seed(123)
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

apply(t(y), 2, mean)


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


