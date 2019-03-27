# IPM from data
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
surv_all    <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat )


surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2,
                          log_area_t03= log_area_t0^3 ) %>% 
                  left_join( clim_mat ) 

sl_grow     <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2 ) %>% 
                  # only seedling GROWTH
                  subset( surv_t1 == 1 ) %>% 
                  left_join( clim_mat ) 

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset( !(stage_t0 %in%  "SL") ) %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM", "SL")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2 ) %>% 
                  left_join( clim_mat ) 

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  left_join( clim_mat ) 

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) 


# models ---------------------------------------------------------
mod_s   <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + 
                (1 | year) + (1 | newid) + (1|location),
                data=surv_all, family='binomial')
mod_g    <- lmer( log_area_t1 ~ log_area_t0 + 
                  (1 | year) + (1 | location) + (1|newid),
                  data=grow)

mod_fl   <- glmer(flow_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | year) + (1 | location) + (1 | newid),
                  data=flow, family='binomial')
mod_fr   <- glmer(numrac_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | year) + (1 | location) + (1 | newid),
                  data=fert, family='poisson')

add_level <- function(x,name_re){
  ranef(x)$newid %>% 
    tibble::add_column(.before=1, newid=row.names(.)) %>% 
    setNames( c('newid', name_re) )
}

id_re <- list( add_level(mod_s,  'surv_re'), 
               add_level(mod_g,  'grow_re'), 
               add_level(mod_fl, 'flow_re'), 
               add_level(mod_fr, 'fert_re') ) %>% 
            Reduce(function(...) inner_join(...),.)

select(id_re,-newid) %>% cor

id_re$flow_re %>% hist
id_re$fert_re %>% hist
id_re$surv_re %>% hist
id_re$grow_re %>% hist
