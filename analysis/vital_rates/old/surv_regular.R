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

clim_var    <- 'ppt'

# data format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( stage_t0 != 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0) ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3 )

        
# climate format ----------------------------------------------------------------
years     <- unique(surv$year) + 1 
m_obs     <- 1
m_back    <- 24
expp_beta <- 20 # this for the 

# format climate - need to select climate predictor first 
clim_mat <- subset(clim, clim_var == "tmean") %>% 
              prism_seas_form("tmean", years, m_back, m_obs)  

# if(clim_var == 'prec'){
#   clim_mat <- subset(clim, clim_var == "prate") %>%
#                 mutate( value = replace(value, value < 0, 0) ) %>% 
#                 month_clim_form("precip", years, m_back, m_obs)
# } else{
#   clim_var_input <- clim_var
#   
#   if( clim_var == 'airt') clim_input <- clim else clim_input <- enso
#   
#   clim_mat <- subset(clim_input, clim_var == clim_var_input) %>%
#                 month_clim_form(clim_var_input, years, m_back, m_obs)
# }

# seedling data
surv_clim <- surv %>% 
                # add 1 to get year_t1
                mutate( year = year + 1 ) %>% 
                left_join( clim_mat) %>%
                subset( !is.na(location) )
                # indices for STAN models
# surv_df   <- surv_clim %>% 
#                 mutate( year_i = year %>% as.factor %>% as.numeric,
#                         site_i = location %>% as.factor %>% as.numeric,
#                         avgt0   = surv_clim %>% select(V1:V12) %>% rowSums,
#                         avgtm1  = surv_clim %>% select(V13:V24) %>% rowSums,
#                         avgtm2  = surv_clim %>% select(V25:V36) %>% rowSums )
surv_df   <- surv_clim %>%
                mutate( year_i = year %>% as.factor %>% as.numeric,
                        site_i = location %>% as.factor %>% as.numeric ) 



# climate data
clim_pred <- dplyr::select(surv_df, year) %>%
                inner_join( clim_mat ) %>% 
                unique %>%
                arrange( year ) %>% 
                dplyr::select(-year)




# LASSO
norm_lasso <- brm(surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (log_area_t0 | location) +
                    V1  + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + 
                    V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + V23 + V24, 
                  data=surv_df,
                  family=bernoulli(link='logit'),
                  prior=c(prior(normal(0,2),
                                class="Intercept"),
                          prior(lasso(), 
                                class = "b")),
                  chains=2,iter=4000, cores=3 )

surv_v <- surv_clim$surv_t1
preds  <- select(surv_clim, 
                 log_area_t0, log_area_t02, log_area_t03, 
                 V1:V12) 

fit   <- glmnet(x= as.matrix(preds), y=surv_v, family = "binomial",alpha = 0)
cvfit <- cv.glmnet(x= as.matrix(preds), y=surv_v, family = "binomial", n=15)
plot(cvfit)

opt.lam <- c(cvfit$lambda.min, cvfit$lambda.1se) 
coef(cvfit, s = opt.lam)

coef(cvfit, s = opt.lam)[,'2'][-c(1:4)] %>% plot(type='b')
abline(h=0)

# fit stan models ----------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n      = nrow(surv_df),
  n_year = surv_df$year_i %>% unique %>% length,
  yr_bck = m_back / 4,
  # yr_bck = m_back / 12,
  n_site = surv_df$site_i %>% unique %>% length,
  n_lag  = ncol(clim_pred),
  y      = surv_df$surv_t1,
  x_size = surv_df$log_area_t0,
  x_size2= surv_df$log_area_t0^2,
  clim   = clim_pred,
  clim_means = rowMeans(clim_pred),
  year_i = surv_df$year_i,
  site_i = surv_df$site_i,
  expp_beta = expp_beta,
  
  # climate variables
  clim1  = t(clim_pred)[1:4 ,],
  clim2  = t(clim_pred)[5:8,],
  clim3  = t(clim_pred)[9:12,],
  K      = ncol(clim_pred) / 4, # four seasons
  M      = 4
  # K      = ncol(clim_pred) / 12,
  # M      = 12
)

# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter   = 4000, 
  thin   = 2, 
  chains = 4
)

# update data list
dat_stan$clim1        <- t(clim_pred)[1:4 ,]
dat_stan$clim2        <- t(clim_pred)[5:8,]
dat_stan$clim3        <- t(clim_pred)[9:12,]
dat_stan$M            <- 4
dat_stan$K            <- ncol(dat_stan$clim) / dat_stan$M


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


out <- rstan::extract(fit_null) %>% bind_cols
write.csv(out, 'surv_36months.csv', row.names=F)

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
