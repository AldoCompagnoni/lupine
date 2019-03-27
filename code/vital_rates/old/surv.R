setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(brms)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# set rstan options to parallel cores
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# data
lupine_df   <- read.csv("data/lupine_all.csv")
enso        <- read.csv("data/enso_data.csv")
# clim        <- read.csv("data/lupine_fc_vars.csv")
clim        <- read.csv("data/prism_point_reyes_87_18.csv")

# 
clim_var    <- 'ppt'

# data format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( stage_t0 != 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 )

# climate format ----------------------------------------------------------------
years     <- unique(surv$year)
m_obs     <- 5
m_back    <- 36
expp_beta <- 20 # this for the 

# format climate - need to select climate predictor first 
clim_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs)
# clim_mat <- subset(enso, clim_var == "oni") %>%
#               oni_form(years, m_back, m_obs)


# seedling data
surv_clim <- surv %>% 
                # add 1 to get year_t1
                left_join( clim_mat ) %>%
                subset( !is.na(location) )
                # indices for STAN models

surv_df   <- surv_clim %>% 
                mutate( year_i  = year %>% as.factor %>% as.numeric,
                        site_i  = location %>% as.factor %>% as.numeric,
                        avgt0   = surv_clim %>% select(V1:V12) %>% rowSums,
                        avgtm1  = surv_clim %>% select(V13:V24) %>% rowSums,
                        avgtm2  = surv_clim %>% select(V25:V36) %>% rowSums )

# climate data
clim_pred <- dplyr::select( surv_df, year ) %>%
                inner_join( clim_mat ) %>% 
                unique %>%
                arrange( year ) %>% 
                dplyr::select( -year )

# test correspondence between years in clim_pred and years in surv_df.
surv_clim %>% 
  mutate( year_i  = year %>% as.factor %>% as.numeric,
          year_r  = year %>% as.factor ) %>% 
  select( year_i, year_r ) %>% 
  unique

# ML model selection -------------------------------------------------
mod             <- glmer(surv_t1 ~ log_area_t0 + (1 | year_i), data=surv_df, family='binomial')
mod2            <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year_i), data=surv_df, family='binomial')
mod_size        <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year_i), data=surv_df, family='binomial')
mod_size2       <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t02 | year_i), data=surv_df, family='binomial')
mod_size_size2  <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 + log_area_t02 | year_i), data=surv_df, family='binomial')

# best model is most complex
AIC(mod,
    mod2,
    mod_size,
    mod_size2,
    mod_size_size2)

# fit stan models ----------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n       = nrow(surv_df),
  n_year  = surv_df$year_i %>% unique %>% length,
  yr_bck  = m_back / 12,
  n_site  = surv_df$site_i %>% unique %>% length,
  n_lag   = ncol(clim_pred),
  y       = surv_df$surv_t1,
  x_size  = surv_df$log_area_t0,
  x_size2 = surv_df$log_area_t0^2,
  clim    = clim_pred,
  clim_means = rowMeans(clim_pred),
  year_i  = surv_df$year_i,
  site_i  = surv_df$site_i,
  expp_beta = expp_beta,
  
  # climate variables
  clim1  = t(clim_pred)[1:12 ,],
  clim2  = t(clim_pred)[13:24,],
  clim3  = t(clim_pred)[25:36,],
  clim1_means = rowMeans( clim_pred[,1:12] ),
  clim2_means = rowMeans( clim_pred[,13:24] ),
  clim3_means = rowMeans( clim_pred[,25:36] ),
  K      = ncol( clim_pred ) / 12,
  M      = 12
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter   = 4000, 
  thin   = 2, 
  chains = 3
)

# update data list
dat_stan$clim1        <- t(clim_pred)[1:12 ,]
dat_stan$clim2        <- t(clim_pred)[13:24,]
dat_stan$clim3        <- t(clim_pred)[25:36,]
dat_stan$K            <- ncol(dat_stan$clim) / 12
dat_stan$M            <- 12

# power exponential moving window
fit_36_nest <- stan(
  # file = paste0("analysis/stan/surv/bernoulli_dirichlet_nest_12.stan"),
  file = paste0("analysis/stan/surv/bernoulli_dirichlet_nest.stan"),
  data = dat_stan,
  # pars = c('theta_m', 'b0', 'b_size', 'b_size2', 'b_c'),
  pars = c('theta_y', 'theta_m', 'b0', 'b_size', 'b_size2', 'b_c'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# power exponential moving window
fit_12 <- stan(
  # file = paste0("analysis/stan/surv/bernoulli_dirichlet_nest_12.stan"),
  # file = paste0("analysis/stan/surv/bernoulli_dirichlet_size2_re_12.stan"),
  file = paste0("analysis/stan/surv/bernoulli_dirichlet_size2_re_resize2_12.stan"),
  data = dat_stan,
  # pars = c('theta_m', 'b0', 'b_size', 'b_size2', 'b_c'),
  # pars = c('theta_m', 'b0', 'b_size', 'b_size2', 'b_c', 'b_yr'),
  pars = c('theta_m', 'b0','b_c', 'b_yr', 'b_size_yr', 'b_size2_yr'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


out <- rstan::extract(fit_null) %>% bind_cols
write.csv(out, 'surv_36months.csv', row.names=F)


# Multiple years, weighted with simplex
fit_3yr <- stan(
  file = paste0("analysis/stan/surv/bernoulli_mYears_simplex_site_re.stan"),
  data = dat_stan,
  pars = c('theta_k', 'b0', 's_yr', 'b_site', 'b_size',  'b_size2', 'b_c', 'b_yr'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

init_t <- Sys.time()
# Null model using size2 and site fixef
fit_avg <- stan(
  file = paste0("analysis/stan/surv/bernoulli_avg_re.stan"),
  data = dat_stan,
  pars = c( 'b_size', 'b_c', 'b_yr', 'b0', 's_yr'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)
Sys.time() - init_t

# store results
out_summ  <- summary(fit_3yr)$summary
out_post  <- extract(fit_3yr) %>% as.data.frame

write.csv(out_summ, paste0("results/surv_sam_3y_summ_",clim_var,".csv"))
write.csv(out_post, paste0("results/surv_sam_3y_post_",clim_var,".csv"),row.names=F)


# # power exponential moving window
# fit_null <- stan(
#   file = paste0("analysis/stan/surv/bernoulli_null.stan"),
#   data = dat_stan,
#   pars = c('b0', 'b_size'),
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains#,
#   #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
# )
# 
# out <- rstan::extract(fit_null) %>% bind_cols
# write.csv(out, 'surv_36months.csv', row.names=F)
