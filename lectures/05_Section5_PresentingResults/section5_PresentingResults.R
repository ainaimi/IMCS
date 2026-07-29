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

pacman::p_load(
  tidyverse,     
  dplyr, 
  purrr,
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

if (!requireNamespace("looplot", quietly = TRUE)) {
  remotes::install_github("matherealize/looplot")
}


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


ggsave(here("_images", "nested_loop_plot.pdf"), plot = p, width = 8, height = 6)

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


ggsave(here("_images", "nested_loop_plot2.pdf"), plot = p, width = 10, height = 6)

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
  mutate(include_flag = if_else(lcl < -.5 & ucl > -.5, "Include", "Exclude"))

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

## ----tidy = F, warning = F, message = F---------------------------------------

perf <- relhaz %>%
  group_by(n, baseline, model) %>%
  summarize(bias = mean(theta - (-0.5)),
            coverage = mean(include_flag == "Include"),
            .groups = "drop")

perf


## ----tidy = F, warning = F, message = F---------------------------------------

p <- perf %>%
  ggplot(.) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_segment(aes(x = model, xend = model, y = 0, yend = bias)) +
  geom_point(aes(x = model, y = bias), size = 2) +
  ylab("Bias in the log Hazard Ratio") +
  xlab("Method") +
  coord_flip() +
  facet_grid(n ~ baseline)

ggsave(here("_images", "lollipop_plot_bias.pdf"), p, width = 8, height = 5)


## ----lollipop1, out.width="10cm", fig.align='center', fig.cap="Lollipop plot of the bias in the log hazard ratio for each method, by sample size (rows) and true baseline hazard (columns).", echo=F----
knitr::include_graphics(here("_images", "lollipop_plot_bias.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

p <- perf %>%
  ggplot(.) +
  geom_hline(yintercept = 0.95, lty = 2) +
  geom_segment(aes(x = model, xend = model, y = 0.95, yend = coverage)) +
  geom_point(aes(x = model, y = coverage), size = 2) +
  ylab("95% Confidence Interval Coverage") +
  xlab("Method") +
  coord_flip() +
  facet_grid(n ~ baseline)

ggsave(here("_images", "lollipop_plot_coverage.pdf"), p, width = 8, height = 5)


## ----lollipop2, out.width="10cm", fig.align='center', fig.cap="Lollipop plot of 95\\% confidence interval coverage for each method, by sample size (rows) and true baseline hazard (columns). The dashed line marks the nominal 0.95 level.", echo=F----
knitr::include_graphics(here("_images", "lollipop_plot_coverage.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

p <- relhaz %>%
  filter(n == 50, baseline == "Exponential") %>%
  ggplot(.) +
  geom_histogram(aes(x = theta), bins = 20,
                 color = "black", fill = "grey") +
  geom_vline(xintercept = -.5, lty = 2, color = "red") +
  xlab("Estimated log Hazard Ratio") +
  ylab("Count") +
  facet_wrap(~model)

ggsave(here("_images", "estimate_histogram.pdf"), p, width = 8, height = 4)


## ----histogram1, out.width="12cm", fig.align='center', fig.cap="Histograms of the estimated log hazard ratio across 100 simulated datasets with $n = 50$ and an exponential baseline hazard. The dashed red line marks the true value of $-0.5$.", echo=F----
knitr::include_graphics(here("_images", "estimate_histogram.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

p <- relhaz %>%
  filter(baseline == "Exponential") %>%
  ggplot(.) +
  geom_density(aes(x = theta, color = model)) +
  geom_vline(xintercept = -.5, lty = 2) +
  scale_color_brewer(palette = "Set1") +
  xlab("Estimated log Hazard Ratio") +
  ylab("Density") +
  facet_wrap(~n, labeller = label_both)

ggsave(here("_images", "estimate_density.pdf"), p, width = 8, height = 4)


## ----density1, out.width="12cm", fig.align='center', fig.cap="Kernel density estimates of the estimated log hazard ratio for each method under an exponential baseline hazard, with 50 (left) and 250 (right) observations per dataset. The dashed line marks the true value of $-0.5$.", echo=F----
knitr::include_graphics(here("_images", "estimate_density.pdf"))

## ----tidy = F, warning = F, message = F---------------------------------------

relhaz_wide <- relhaz %>%
  select(dataset, n, baseline, model, theta) %>%
  pivot_wider(names_from = model, values_from = theta)

head(relhaz_wide)


## ----tidy = F, warning = F, message = F---------------------------------------

p <- relhaz_wide %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_point(aes(x = Cox, y = `RP(2)`), alpha = .5) +
  xlab("Estimated log Hazard Ratio: Cox") +
  ylab("Estimated log Hazard Ratio: Royston-Parmar") +
  facet_grid(n ~ baseline)

ggsave(here("_images", "estimate_scatter.pdf"), p, width = 8, height = 6)


## ----scatter1, out.width="10cm", fig.align='center', fig.cap="Estimated log hazard ratios from the Cox model plotted against those from the Royston-Parmar flexible parametric model, one point per simulated dataset. The dashed line indicates perfect agreement.", echo=F----
knitr::include_graphics(here("_images", "estimate_scatter.pdf"))

