rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(mgcv)
library(lme4)
library(bbmle)
library(ggplot2)
library(parallel)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')

# data
lupine_df   <- read.csv("data/lupine_all.csv")
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")

# data format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  # subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0  = log_area_t0_z ) %>% 
                  mutate( log_area_t02 = log_area_t0_z^2,
                          log_area_t03 = log_area_t0_z^3 )

# climate format ----------------------------------------------------------------
years     <- 1990:2018
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies from monthly anomalies
year_m_anom <- function(x, var ){

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

# calculate yearly anomalies
year_anom <- function(clim_x, clim_var = "ppt", 
                      years, m_back, m_obs ){
  
  # "spread" the 12 months
  clim_m <- select(clim_x, -clim_var )
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(years, clim_back, m_obs, clim_m)
  x_clim    <- year_by_month_mat(clim_x_l, years) %>% 
                  gather(month,t0,V1:V12) %>% 
                  select(year,month,t0) %>% 
                  mutate( month = gsub('V','',month) ) %>%
                  mutate( month = as.numeric(month) )
  
  if( clim_var == 'ppt'){
    raw_df <- x_clim %>% 
                group_by(year) %>% 
                summarise( ppt_t0 = sum(t0) ) %>% 
                ungroup %>% 
                arrange( year ) %>% 
                mutate( ppt_t0  = scale(ppt_t0)[,1] ) %>% 
                mutate( ppt_tm1 = lag(ppt_t0) )
  }
  if( clim_var == 'tmp'){
    raw_df <- x_clim %>% 
                group_by(year) %>% 
                summarise( tmp_t0 = mean(t0) ) %>% 
                ungroup %>% 
                arrange( year ) %>% 
                mutate( tmp_t0  = scale(tmp_t0)[,1] ) %>% 
                mutate( tmp_tm1 = lag(tmp_t0) )
  }
  
  raw_df
  
}


# format climate - need to select climate predictor first 
ppt_mat <- subset(clim, clim_var == "ppt") %>%
              year_anom("ppt", years, m_back, m_obs)
              
tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              year_anom("tmp", years, m_back, m_obs) 
  
enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_m_anom('oni')

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat,tmp_mat,enso_mat) )

# demography plus clim
surv_clim <- left_join(surv, clim_mat) %>%
                subset( !is.na(location) )


# fit climate models -------------------------------------------------
climate_mods <- list(

  # null 
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + (log_area_t0 | year) + (1 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + ppt_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + ppt_tm1 + (log_area_t0 | year) + (1 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + (log_area_t0 | year) + (1 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + oni_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + oni_tm1 + (log_area_t0 | year) + (1 | location)
  
)

# years to do crossvalidation on
years_cv  <- unique(surv$year)
mod_nam   <- c( 'null', 
                'ppt_t0',  'ppt_tm1', 
                'tmp_t0',  'tmp_tm1',
                'oni_t0',  'oni_tm1',
                'spei_to', 'spei_tm1' )

# function to carry out crossvalidation
crossval <- function(yr_ii){ #length(years_cv)
  
  # set up dataframes to train and test year-specific predictions
  train_df <- subset(surv_clim, !(year %in% yr_ii ) )
  test_df  <- subset(surv_clim,   year == yr_ii ) %>% 
                mutate( resp = surv_t1 )

  # fit all models
  mods_clim <- lapply( climate_mods, 
                         function(x) glmer(x, data=train_df, family='binomial') ) 
  
  pred_df <- list()
  for(ii in 1:length(mods_clim)){
    
    pred_df[[ii]] <- test_df %>% 
                        # calculate predictions
                        mutate( pred = predict(mods_clim[[ii]], 
                                               newdata=test_df, 
                                               type = 'response',
                                               re.form=NA) 
                              ) %>% 
                        select( year, pred, resp ) %>% 
                        mutate( model = mod_nam[ii] )
    
  }
  
  pred_df %>% bind_rows 
  
}


# parallel processing model fit -------------------------------------------------------

# detect cores
cores     <- (detectCores() - 1)
# set up as many clusters as detected by 'detectCores()'
cluster   <- parallel::makePSOCKcluster(cores)

# attach packages that will be needed on each cluster
clusterEvalQ(cluster, list(library(lme4), library(boot),  
                           library(dplyr) , library(tidyr)) )

# attach objects that will be needed on each cluster
clusterExport(cluster, c('surv_clim', 'climate_mods', 'mod_nam') )

# Fit 5 alternative models on 200 bootstrapped samples
init_t        <- Sys.time()
boostr_par_l  <- parLapply(cluster, years_cv, crossval)
Sys.time() - init_t

# get scores
score_df <- boostr_par_l %>% 
              bind_rows %>% 
              # calculate negative log likelihood and Squared Error
              mutate( log_lik = dbinom(x    = resp, 
                                       size = 1, 
                                       prob = pred,
                                       log  = T ) * -1 ) %>% 
              mutate( se      = (resp - pred)^2 ) %>% 
              group_by( model ) %>% 
              # sum of neg_log_lik/ root mean squared error  
              summarise( log_lik_sum = sum(log_lik),
                         rmse        = mean(se) %>% sqrt ) %>% 
              ungroup %>% 
              as.data.frame()


# logarithmic score
score_df %>% 
  arrange( log_lik_sum ) %>% 
  select( -rmse ) %>% 
  mutate( dScore = log_lik_sum - log_lik_sum[1] ) %>% 
  write.csv('results/vital_rates/crossval/surv_cv_logscore.csv',
            row.names=F)

# root mean squared error
score_df %>% 
  arrange( rmse ) %>% 
  select( -log_lik_sum ) %>% 
  mutate( dScore = rmse - rmse[1] ) %>% 
  write.csv('results/vital_rates/crossval/surv_cv_rmse.csv',
            row.names=F)

save.image( 'results/vital_rates/crossval/surv_cv.Rdata' )
