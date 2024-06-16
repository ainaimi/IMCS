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


## ----tidy = F, warning = F, message = F---------------------------------------

pacman::p_load(
  tidyverse,     
  dplyr, 
  purr,
  magrittr
  )

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

set.seed(123)
params = list(
    samplesize = c(100, 200, 500),
    param1 = c(1, 2), 
    param2 = c(1, 2, 3), 
    param3 = c(1, 2, 3, 4)
)

design = expand.grid(params)

# add some "results"
design %<>% 
    mutate(method1 = rnorm(n = n(),
                           mean = param1 * (param2 * param3 + 1000 / samplesize), 
                           sd = 2), 
           method2 = rnorm(n = n(),
                           mean = param1 * (param2 + param3 + 2000 / samplesize), 
                           sd = 2), 
           method3 = rnorm(n = n(),
                           mean = param1 * (param2 + param3 + 3000 / samplesize), 
                           sd = 2))

knitr::kable(head(design, n = 10))

## ----tidy = F, warning = F, message = F, results='hide'-----------------------

remotes::install_github("matherealize/looplot")


## ----tidy = F, warning = F, message = F---------------------------------------

pacman::p_load(looplot)

p = nested_loop_plot(resdf = design, 
                     x = "samplesize", steps = c("param2", "param3"),
                     grid_rows = "param1", 
                     steps_y_base = -10, steps_y_height = 3, steps_y_shift = 10,
                     x_name = "Sample Size", y_name = "Error",
                     spu_x_shift = 200,
                     steps_values_annotate = TRUE, steps_annotation_size = 3, 
                     hline_intercept = 0, 
                     y_expand_add = c(10, NULL), 
                     post_processing = list(
                        add_custom_theme = list(
                            axis.text.x = element_text(angle = -90,
                                                       vjust = 0.5,
                                                       size = 8)
                        ))
                     )


ggsave(here("_images", "nested_loop_plot.pdf"), width = 8, height = 6)

## ----nestedloop1, out.width="12cm", fig.align='center', fig.cap="Example Nested Loop Plot of Hypothetical Simulation Results.", echo=F----
knitr::include_graphics(here("_images", "nested_loop_plot.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

pacman::p_load(looplot)

p = nested_loop_plot(resdf = design, 
                     x = "samplesize", steps = c("param1", "param2", "param3"),
                     #grid_rows = "param1", 
                     steps_y_base = -10, steps_y_height = 3, steps_y_shift = 10,
                     x_name = "Sample Size", y_name = "Error",
                     spu_x_shift = 200,
                     steps_values_annotate = TRUE, steps_annotation_size = 3, 
                     hline_intercept = 0, 
                     y_expand_add = c(10, NULL), 
                     post_processing = list(
                        add_custom_theme = list(
                            axis.text.x = element_text(angle = -90,
                                                       vjust = 0.5,
                                                       size = 8)
                        ))
                     )


ggsave(here("_images", "nested_loop_plot2.pdf"), width = 10, height = 6)

## ----nestedloop2, out.width="12cm", fig.align='center', echo=F----------------
knitr::include_graphics(here("_images", "nested_loop_plot2.pdf"))

## ----tidy = F, warning = F, message = F, echo = F-----------------------------
pacman::p_load(rsimsum)

data("relhaz", package = "rsimsum")

head(relhaz)

dim(relhaz)

table(relhaz$n)
table(relhaz$baseline)
table(relhaz$model)


## ----tidy = F, warning = F, message = F---------------------------------------

relhaz %>% 
  group_by(n, baseline, model) %>% 
  count()


## ----tidy = F, warning = F, message = F---------------------------------------

relhaz <- relhaz %>% 
  mutate(lcl = theta - 1.96*se,
         ucl = theta + 1.96*se)

head(relhaz)


## ----tidy = F, warning = F, message = F---------------------------------------
relhaz <- relhaz %>% 
  mutate(include_flag = if_else(lcl<-.5 & ucl>-.5, "Include", "Exclude"))

## ----tidy = F, warning = F, message = F---------------------------------------

p <- relhaz %>% 
  filter(n == 50, baseline == "Exponential") %>% 
  ggplot(.) + 
  geom_hline(yintercept = -.5, lty = 2) +
  geom_pointrange(aes(x = dataset, 
                      y = theta, 
                      ymin = lcl, 
                      ymax = ucl, color = include_flag), 
                  size = .2, 
                  alpha = .75) + 
  scale_color_manual(values=c("red","grey")) +
  ylab("log Hazard Ratio") + 
  xlab("Sample Number") +
  coord_flip() +
  theme(legend.position = "none", text=element_text(size=12)) +
  facet_wrap(~model)

ggsave(here("_images", "zip_plot_version1.pdf"), p)


## ----zipper1, out.width="10cm", fig.align='center', fig.cap="Zipper plot displaying the distribution of normal-interval (Wald) confidence intervals in the relhaz data.", echo=F----
knitr::include_graphics(here("_images", "zip_plot_version1.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

relhaz <- relhaz %>% 
  mutate(test_statistic = abs(theta/se))


## ----tidy = F, warning = F, message = F---------------------------------------

p <- relhaz %>% 
  filter(n == 50, baseline == "Exponential") %>% 
  ggplot(.) + 
  geom_hline(yintercept = -.5, lty = 2) +
  geom_pointrange(aes(x = test_statistic, 
                      y = theta, 
                      ymin = lcl, 
                      ymax = ucl, color = include_flag), 
                  size = .2, 
                  alpha = .75) + 
  scale_color_manual(values=c("red","grey")) +
  ylab("log Hazard Ratio") + 
  xlab("Wald Null Test Statistic") +
  coord_flip() +
  theme(legend.position = "none", text=element_text(size=12)) +
  facet_wrap(~model)

ggsave(here("_images", "zip_plot_version2.pdf"), p)


## ----zipper2, out.width="10cm", fig.align='center', fig.cap="Zipper plot displaying the distribution of normal-interval (Wald) confidence intervals in the relhaz data. Bounds are ranked according to the magnitude of the Wald test statistic for each point estimate.", echo=F----
knitr::include_graphics(here("_images", "zip_plot_version2.pdf"))

