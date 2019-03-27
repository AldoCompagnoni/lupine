# find out if it makes sense to model survival separating
# between adults and seedlings. Alternative is modeling everything together.
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"

all_indiv_sample <- c("BS (7)", 'DR (3)', 'NB (2)')

# data
lupine_df   <- read.csv( "data/lupine_all.csv") #%>% 
                  # subset( location %in% c("BS (7)") ) %>%
                  # subset( year > 2008 )
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")

# format climate data ----------------------------------------
years     <- c(2005:2018)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  # set names of climate variables
  clim_names <- paste0( var,c('_t0','_tm1','_t0_tm1','_t0_tm2') )
  
  mutate(x, 
         avgt0     = x %>% select(V1:V12) %>% rowSums,
         avgtm1    = x %>% select(V13:V24) %>% rowSums,
         avgt0_tm1 = x %>% select(V1:V24) %>% rowSums,
         avgt0_tm2 = x %>% select(V1:V36) %>% rowSums ) %>% 
    select(year, avgt0, avgtm1, avgt0_tm1, avgt0_tm2) %>% 
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

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat,tmp_mat,enso_mat) )


# vital rates format --------------------------------------------------------------

# survival of all individuals
surv_all    <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

# survival of adults
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

# survival of seedlings
seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2,
                          log_area_t03= log_area_t0^3 ) %>% 
                  left_join( clim_mat ) 

# models -----------------------------------------------------

# seedlings
mod_sl   <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                    (1 | year) + (0 + log_area_t0 | year),
                  data=seedl, family='binomial')

# only adults
mod_a    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) , 
                  data=surv, family='binomial')

# all individuals
mod_all  <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) , 
                  data=surv_all, family='binomial')


# calculate AICs --------------------------------------------

# extract log-likelihood
log_lik <- function(x) x %>% summary %>% .$logLik %>% as.numeric

# calculate AICs
aic_sep <- -2*(log_lik(mod_sl) + log_lik(mod_s)) + (2*(7+7))
aic_all <- AIC(mod_all)

# separate modeling is better, despite the 7 additional parameters
aic_sep < aic_all
