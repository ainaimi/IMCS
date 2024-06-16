# purpose: create nested loop plots for the HTE simulation project
# by Echo Jan 25, 2024
# plan: 
# step 1:methods: compare DR learner and causal forest, both with super learner
# step 2:create the table for simulation results with confounder number, sample size and true effects
#        for bias, MSE and coverage 
# step 3:use the code in the link https://matherealize.github.io/looplot_demo.html
#                                 https://matherealize.github.io/looplot_gallery.html


# install.packages("devtools")
# devtools::install_github("matherealize/looplot", force = TRUE)

pacman::p_load(tidyverse, ggplot2, looplot,reshape, data.table) #other package needed: here

nested_plot <- function(metric = "coverage_prob", y_name = "Coverage",estimate = "cate_m1"){
  # step 2: create the table
  # read in the output files using regular expression
  setwd("H:/RA/HTE/output/sim1000_n10000")
  results_big <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv)
  setwd("H:/RA/HTE/output/sim1000_n1000")
  results_small <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv)
  results <- rbind(results_big, results_small)
  
  # read in the truth file
  truth_data <- read_csv(
    here::here("data", "truth_data_simulation_all_update.csv")
  )
  
  # get the truth data for each group
  truth <- truth_data %>%
    filter(c_dim == 15) %>%
    mutate(group = row_number(), RD_m0 = round(risk_diff_m0, 3), RD_m1 = round(risk_diff_m1, 3)) %>%
    select(group, RD_m0, RD_m1)
  
  # merge the truth data with the result
  plot_result <- results %>%
    mutate(truth = round(truth, 3)) %>%
    left_join(truth, by = "group")
  
  # double check if the truth is the same as we expected for each group
  #table(plot_result[(grep("_cate_m0$", plot_result$estimator)),]$truth) # not exactly same
  #table(plot_result[(grep("_cate_m1$", plot_result$estimator)),]$truth) # exactly same
  
  # subset data to only the targeted estimate (ate vs. cate) and evaluation metrics 
  plot_result2 <- plot_result[(grep(paste0(estimate, "$"), plot_result$estimator)),] %>%
    select(estimator, group, RD_m0, RD_m1, sample_size, confounder_number, all_of(metric)) %>%
    mutate(sample_size = factor(sample_size))%>%
    filter(complete.cases(.)) # oracle doesn;t have value for coverage 
  # transform data from long to wide for plots
  plot_data <- data.table::dcast(setDT(plot_result2), RD_m0+RD_m1+sample_size+confounder_number~estimator, value.var = metric)
  
  
  # plot_all <- plot_result %>%
  #   mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
  #          ifelse(grepl("_cate_m0$", estimator), "cate_m0",
  #          "cate_m1"))) 
  
  
  # step 3:
  if(metric=="coverage_prob"){
    p = nested_loop_plot(resdf = plot_data,
                         x = "sample_size", steps = c("confounder_number","RD_m1","RD_m0"),
                         steps_y_base = -0.1, steps_y_height = 0.05, steps_y_shift = 0.05,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(0, 1, by = 0.25),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                         hline_intercept = 0,
                         y_expand_add = c(0.1, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3"), # https://r-graph-gallery.com/ggplot2-color.html
                         line_linetypes = c(1, 2, 3),
                         point_shapes = c(1,2,3),
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 5)
                           )
                         ))
  }else if(metric=="bias_mean"){
    p = nested_loop_plot(resdf = plot_data,
                         x = "sample_size", steps = c("confounder_number","RD_m1","RD_m0"),
                         steps_y_base = -0.045, steps_y_height = 0.003, steps_y_shift = 0.003,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(-0.02, 0.02, by = 0.01),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                         hline_intercept = -0.04,
                         y_expand_add = c(0.005, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3", "yellow"),
                         line_linetypes = c(1,2,3,4),
                         point_shapes = c(1,2,3,4),
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 5)
                           )
                         ))
  } else{
    p = nested_loop_plot(resdf = plot_data,
                         x = "sample_size", steps = c("confounder_number","RD_m1","RD_m0"),
                         steps_y_base = -0.0006, steps_y_height = 0.0002, steps_y_shift = 0.0005,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(0, 0.005, by = 0.001),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                         hline_intercept = 0,
                         y_expand_add = c(0.0005, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3", "yellow"),
                         line_linetypes = c(1,2,3,4),
                         point_shapes = c(1,2,3,4),
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 5)
                           )
                         ))
  }
  
  
  # figure in a panel
  # p2 = nested_loop_plot(resdf = plot_data, 
  #                      x = "sample_size", steps = c("confounder_number"),
  #                      grid_rows = "RD_m1", grid_cols = "RD_m0",
  #                      steps_y_base = -0.1, steps_y_height = 0.1, 
  #                      x_name = "Sample Size", y_name = y_name, ylim = c(-0.5,1),
  #                      spu_x_shift = 1000,
  #                      steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
  #                      hline_intercept = 0, 
  #                      y_expand_add = c(0.2, NULL), 
  #                      post_processing = list(
  #                        add_custom_theme = list(
  #                          axis.text.x = element_text(angle = -90, 
  #                                                     vjust = 0.5, 
  #                                                     size = 5) 
  #                        )
  #                      ))
  # p3 = nested_loop_plot(resdf = plot_data,
  #                       x = "sample_size", steps = c("RD_m1","RD_m0"),
  #                       grid_cols = "confounder_number",
  #                       steps_y_base = -0.1, steps_y_height = 0.1,
  #                       x_name = "Sample Size", y_name = y_name,
  #                       spu_x_shift = 1000,
  #                       steps_values_annotate = TRUE, steps_annotation_size = 2.5,
  #                       hline_intercept = 0,
  #                       y_expand_add = c(0.2, NULL),
  #                       post_processing = list(
  #                         add_custom_theme = list(
  #                           axis.text.x = element_text(angle = -90,
  #                                                      vjust = 0.5,
  #                                                      size = 5)
  #                         )
  #                       ))
  ggsave(paste0(y_name, "_", estimate, ".png"), p, path = here::here("output"), width = 14, height = 7)
}

nested_plot(metric = "coverage_prob", y_name = "Coverage",estimate = "ate")
nested_plot(metric = "coverage_prob", y_name = "Coverage",estimate = "cate_m0")
nested_plot(metric = "coverage_prob", y_name = "Coverage",estimate = "cate_m1")

nested_plot(metric = "bias_mean", y_name = "Bias",estimate = "ate")
nested_plot(metric = "bias_mean", y_name = "Bias",estimate = "cate_m0")
nested_plot(metric = "bias_mean", y_name = "Bias",estimate = "cate_m1")

nested_plot(metric = "MSE", y_name = "MSE",estimate = "ate")
nested_plot(metric = "MSE", y_name = "MSE",estimate = "cate_m0")
nested_plot(metric = "MSE", y_name = "MSE",estimate = "cate_m1")
 