library(dismo)
library(dplyr)
library(tidyr)
library(rstan)
options(stringsAsFactors = F)

# arguments from command line
args  <- commandArgs()

# load arguments
code_dir  <- args[7]
data_dir  <- args[8]
out_dir   <- args[9]
clim_var  <- args[10]

print(args)

# set rstan options to parallel cores
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
 
# read format data functions
source('/home/compagna/lupine/format_functions.R')

# data
lupine_df   <- read.csv( paste0(data_dir, "lupine_all.csv") )
enso        <- read.csv( paste0(data_dir, "enso_data.csv") )
clim        <- read.csv( paste0(data_dir, "lupine_fc_vars.csv") )

# data format --------------------------------------------------------------
grow        <- lupine_df %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1 = log(area_t1),
                          log_area_t0 = log(area_t0),
                          year        = year + 1 ) 
                   


# climate format ----------------------------------------------------------------
years     <- unique(grow$year)
m_obs     <- 5
m_back    <- 36
expp_beta <- 20 # this for the 

# format climate - need to select climate predictor first 
clim_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs)
# clim_mat <- subset(enso, clim_var == "oni") %>%
#               oni_form(years, m_back, m_obs)


# growth data
grow_clim <- left_join(grow, clim_mat) %>%
                subset( !is.na(location) )
                # indices for STAN models
grow_df   <- grow_clim %>% 
                mutate( year_i = year %>% as.factor %>% as.numeric,
                        site_i = location %>% as.factor %>% as.numeric,
                        avgt0   = grow_clim %>% select(V1:V12) %>% rowSums,
                        avgtm1  = grow_clim %>% select(V13:V24) %>% rowSums,
                        avgtm2  = grow_clim %>% select(V25:V36) %>% rowSums )

# climate data
clim_pred <- dplyr::select(grow_df, year) %>%
                inner_join( clim_mat ) %>% 
                unique %>%
                arrange( year ) %>% 
                dplyr::select(-year)

# fit stan models ----------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n      = nrow(grow_df),
  n_year = grow_df$year_i %>% unique %>% length,
  yr_bck = m_back / 12,
  n_site = grow_df$site_i %>% unique %>% length,
  n_lag  = ncol(clim_pred),
  y      = grow_df$log_area_t1,
  x_size = grow_df$log_area_t0,
  x_size2= grow_df$log_area_t0^2,
  clim   = clim_pred,
  clim_means = rowMeans(clim_pred),
  year_i = grow_df$year_i,
  site_i = grow_df$site_i,
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

# update data list
dat_stan$clim1        <- t(clim_pred)[1:12 ,]
dat_stan$clim2        <- t(clim_pred)[13:24,]
dat_stan$clim3        <- t(clim_pred)[25:36,]
dat_stan$K            <- ncol(dat_stan$clim) / 12
dat_stan$M            <- 12

# power exponential moving window
fit_36_nest <- stan(
  file = paste0("analysis/stan/grow/gaussian_dirichlet_nest.stan"),
  data = dat_stan,
  pars = c('theta_y', 'theta_m', 'b0', 'b_size', 'b_size2', 'b_c'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

fit_12 <- stan(
  file = paste0("analysis/stan/grow/gaussian_dirichlet_12_re_resize.stan"),
  data = dat_stan,
  pars = c('b0', 'b0_size', 'b_yr', 'b_size_yr', 
           'theta_m', 'b_c', 's_yr', 's_size_yr', 's'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)


save.image( paste0(out_dir, '_grow_36_', clim_var, '.Rdata') )

out <- rstan::extract(fit_mod) %>% as.data.frame
# write.csv(out, paste0(out_dir, '_surv_36.csv'), row.names=F)
