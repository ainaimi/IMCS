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


## ----out.width = "10cm",fig.cap="Causal diagram representing exposure induced mediator outcome confounding.",echo=F----
knitr::include_graphics(here("_images","mediation_dag.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

expit <- function(x){1/(1 + exp(-x))}

set.seed(123)

true_CDM <- function(sample_size, exposure, mediator){
  
  n <- sample_size

  U <- rnorm(n)
  
  C <- rnorm(n)
  
  X <- exposure # rbinom(n, size = 1, p = expit(-1 + log(2)*C ))
  
  L <- rnorm(n, mean = 0 + 1.5*X + 1.5*U)
  
  Z <- mediator # rbinom(n, size = 1, expit(-1.5 + log(2)*X + log(2)*L ))
  
  return(mean(10 + 5*X + 5*C + 5*Z + 4*L + 4*U))

}



