# What is climate anomaly in terms of degrees C?
rm(list=ls())
source("analysis/format_data/format_functions.R")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)

# get climate data
clim        <- read.csv("data/prism_point_reyes_87_18.csv")

# calculate ABSOLUTE Celsius yearly temperatures
mean_t_df <- tmp_mat %>% 
               group_by( year) %>% 
               summarise( mean_t = mean(clim_value) ) %>% 
               ungroup

# calculate standardized (z) and absolute anomalies
t_anom_df <- data.frame( t_abs    = mean_t_df[,2,drop=T],
                       t_anom_z = mean_t_df[,2,drop=T] %>% scale %>% .[,1] ) %>% 
              mutate(  t_anom_a = t_abs - mean(mean_t_df[,2,drop=T]) )

# verdict: anomaly of 1 == 0.6 Celsius!