# Stochastic LTRE: lam_s ~ temp. anomaly
rm(list=ls())
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(rstan)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)
library(lintr)
library(goodpractice)


# data
lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                  mutate( log_area_t0  = log(area_t0) ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3 )
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
pred_g      <- read_xlsx('data/post predation_lupinus tidestromii.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")
germ        <- read_xlsx('data/seedbaskets.xlsx') %>% 
                 select(g0:g2) %>% 
                 colMeans
germ_adj    <- read.csv('results/ml_mod_sel/germ/germ_adj.csv')

# format climate data ----------------------------------------
years     <- c(2005:2018)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var){
  
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

# first, format site/year combinations
site_df     <- select(lupine_df, year, location) %>% 
                 unique %>% 
                 # create 'Site' column to merge with consumption dat
                 mutate( Site = gsub(' \\([0-9]\\)','',location) ) %>% 
                 subset( year > 2004 ) %>% 
                 arrange( location, year ) %>% 
                 complete(location,year)

surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          surv_t1, tmp_t0 )

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          log_area_t1, tmp_t0 )

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          flow_t0, tmp_t0 ) %>% 
                  subset( year != 2018 )

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          numrac_t0, tmp_t0 ) %>% 
                  subset( year != 2018 )


# All vital rates at once
vr_all <- Reduce( function(...) full_join(...), list( surv, grow,
                                                      flow, fert ) )

vr_all %>% 
  subset( !is.na(surv_t1) & 
          !is.na(log_area_t1) & 
          !is.na(flow_t0) & 
          !is.na(numrac_t0) ) %>% dim

dat_stan <- list(
  nS = nrow(surv), 
  nG = nrow(grow), 
  nF = nrow(flow), 
  nR = nrow(fert),
  nYears = 13, 
  yearsS = surv$year - 2004, 
  yearsG = grow$year - 2004, 
  yearsF = flow$year - 2004, 
  yearsR = fert$year - 2004,
  n_loc  = 7,
  loc_s = surv$location %>% as.factor %>% as.numeric, 
  loc_g = grow$location %>% as.factor %>% as.numeric, 
  loc_f = flow$location %>% as.factor %>% as.numeric, 
  loc_r = fert$location %>% as.factor %>% as.numeric,
  yS = surv$surv_t1,
  yG = grow$log_area_t1,  
  yF = flow$flow_t0, 
  yR = fert$numrac_t0,
  xS = surv$log_area_t0, xS2 = surv$log_area_t02, xS3 = surv$log_area_t03, 
  xG = grow$log_area_t0, 
  xF = flow$log_area_t0, 
  xR = fert$log_area_t0,
  c_s = surv$tmp_t0,
  c_f = flow$tmp_t0,
  c_r = fert$tmp_t0
)

# test_na <- function(x) x %>% is.na %>% sum

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )


# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 1
)

# quote a series of bare names
quote_bare <- function( ... ){
    substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}


# NULL model (model of the mean)
fit_try <- stan(
  file = 'analysis/vital_rates/bayes/vr_cmpl_smpl.stan',
  data = dat_stan,
  pars = quote_bare(yr_a_s,  yr_a_g,  yr_a_f,  yr_a_r,
                    yr_b_s,  yr_b_g,  yr_b_f,  yr_b_r,
                    loc_a_s, loc_a_g, loc_a_f, loc_a_r,
                    loc_b_s, loc_b_g, loc_b_f, loc_b_r,
                    b_c_s  ,          b_c_f  , b_c_r,
                    b_s2,    b_s3,    sigma_y),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

surv$location %>% is.na %>% sum
grow$location %>% is.na %>% sum
fert$location %>% is.na %>% sum
flow$location %>% is.na %>% sum

# fuck... So much missing data
vr_all$surv_t1 %>% is.na %>% sum
vr_all$log_area_t1 %>% is.na %>% sum
vr_all$flow_t0 %>% is.na %>% sum
vr_all$numrac_t0 %>% is.na %>% sum


#
cat("
  int n123 = 1;
  int n12  = 1;
  int n23  = 1;
  int n31  = 1;
  int n1   = 1;
  int n2   = 1;
  int n3   = 1;

  for (n in 1:N) {
    if (y_observed[n, 1] && y_observed[n, 2] && 
        y_observed[n, 3] ) {
      ns123[n123] = n;
      n123 = n123 + 1;
    } else if (y_observed[n, 1] && y_observed[n, 2]) {
      ns12[n12] = n;
      n12 = n12 + 1;
    } else if (y_observed[n, 2] && y_observed[n, 3]) {
      ns23[n23] = n;
      n23 = n23 + 1;
    } else if (y_observed[n, 1] && y_observed[n, 3]) {
      ns13[n13] = n;
      n13 = n13 + 1;
    }

  }

  
    
    ")