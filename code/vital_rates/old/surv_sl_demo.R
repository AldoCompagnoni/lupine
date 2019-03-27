setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(brms)
library(lme4)
options(stringsAsFactors = F)
# source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# set rstan options to parallel cores
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# data
lupine_df   <- read.csv("data/lupine_all.csv")
enso        <- read.csv("data/enso_data.csv")
clim        <- read.csv("data/lupine_fc_vars.csv")

clim_var    <- 'prec'

# data format --------------------------------------------------------------

# seedling data
seedl       <- dplyr::select(lupine_df, location, year, stage_t0, surv_t1)

seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0) )



# climate format ----------------------------------------------------------------
years     <- unique(seedl$year)
m_obs     <- 7
m_back    <- 36
expp_beta <- 20 # this for the 

# format climate - need to select climate predictor first 
if(clim_var == 'prec'){
clim_mat <- subset(clim, clim_var == "prate") %>%
              mutate( value = replace(value, value < 0, 0) ) %>% 
              month_clim_form("precip", years, m_back, m_obs)
}
if(clim_var == 'airt'){
clim_mat <- subset(clim, clim_var == "airt") %>%
              month_clim_form("airt", years, m_back, m_obs)
}

# seedling data
sl_clim   <- left_join(seedl, clim_mat) %>%
                subset( !is.na(location) ) %>%
                # remove locations from first 3 years ( for now)
                subset( !(location %in% c('R1','R3','stake 25')) )
   
# indices for STAN models
sl_df      <- sl_clim %>% 
                  mutate( year_i = year %>% as.factor %>% as.numeric,
                          site_i = location %>% as.factor %>% as.numeric,
                          avgt0   = sl_clim %>% dplyr::select(V1:V12) %>% rowSums,
                          avgtm1  = sl_clim %>% dplyr::select(V13:V24) %>% rowSums,
                          avgtm2  = sl_clim %>% dplyr::select(V25:V36) %>% rowSums )

# climate data
clim_pred <- dplyr::select(sl_df, year) %>%
                inner_join( clim_mat ) %>% 
                unique %>%
                arrange( year ) %>% 
                dplyr::select(-year)


# fit models ----------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n      = nrow(sl_df),
  n_year = sl_df$year_i %>% unique %>% length,
  yr_bck = m_back / 12,
  n_site = sl_df$site_i %>% unique %>% length,
  n_lag  = ncol(clim_pred),
  y      = sl_df$surv_t1,
  x_size = sl_df$log_area_t0,
  clim   = clim_pred,
  clim_means = rowMeans(clim_pred),
  year_i = sl_df$year_i,
  site_i = sl_df$site_i,
  expp_beta = expp_beta,
  
  # climate variables
  clim1  = t(clim_pred)[1:12 ,],
  clim2  = t(clim_pred)[13:24,],
  clim3  = t(clim_pred)[25:36,],
  clim1_means = rowMeans( clim_pred[,1:12] ),
  clim2_means = rowMeans( clim_pred[,13:24] ),
  clim3_means = rowMeans( clim_pred[,25:36] ),
  K      = ncol(clim_pred) / 12,
  M      = 12
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter   = 4000, 
  thin   = 2, 
  chains = 4
)

# Average of previous 3 years
fit_avg <- stan(
  file = paste0("analysis/stan/sl/bernoulli_null.stan"),
  data = dat_stan,
  pars = c('b0', 'b_size'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# update data list
dat_stan$clim1        <- t(clim_pred)[1:12 ,]
dat_stan$clim2        <- t(clim_pred)[13:24,]
dat_stan$clim3        <- t(clim_pred)[25:36,]
dat_stan$K            <- ncol(dat_stan$clim) / 12
dat_stan$M            <- 12

# power exponential moving window
fit_36_nest <- stan(
  file = paste0("analysis/stan/surv/bernoulli_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'b0', 'b_size', 'b_c'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

out <- rstan::extract(fit_null) %>% bind_cols
write.csv(out, 'surv_sl_36months.csv', row.names=F)
