rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(ggplot2)
library(testthat)
library(bbmle)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")


# climate format ----------------------------------------------------------------
years     <- 1990:2018
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  # set names of climate variables
  clim_names <- paste0( var,c('_t0','_tm1') )
  
  mutate(x, 
         avgt0     = x %>% select(V1:V12) %>% rowSums,
         avgtm1    = x %>% select(V13:V24) %>% rowSums ) %>% 
    select(year, avgt0, avgtm1 ) %>% 
    setNames( c('year',clim_names) )
  
}

# format climate - need to select climate predictor first 
ppt_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs) %>% 
              year_anom('ppt')

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              year_anom('tmp')
  
enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_anom('oni')


tmp_mat %>% 
  select( -year ) %>% 
  cor()

tmp_mat$tmp_t0 %>% plot(type='l')

tmp_mat$tmp_t0 %>% acf 
ppt_mat$ppt_t0 %>% acf
enso_mat$oni_t0 %>% acf

# temperature is strongly auto-correlated with lag 1!!!


setwd('C:/cloud/Dropbox/lupine/results/climate')

tiff('autoregression_airt.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

tmp_mat$tmp_t0 %>% acf(main="Autocorrelation of temperature")

dev.off()
