rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(brms)
library(lme4)
library(bbmle)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# set rstan options to parallel cores
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# data
lupine_08   <- read.csv("data/lupine_05_08.csv")
lupine_18   <- read.csv("data/lupine_08_18.csv")
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")
lupine_df   <- bind_rows(lupine_08,lupine_18)

# data format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2 )

        
# climate format ----------------------------------------------------------------
years     <- unique(surv$year)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  # set names of climate variables
  clim_names <- paste0( var,c('_t0','_tm1') )
  
  mutate(x, 
         avgt0   = x %>% select(V1:V12) %>% rowSums,
         avgtm1  = x %>% select(V13:V24) %>% rowSums ) %>% 
    select(year, avgt0, avgtm1) %>% 
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


# demography plus clim
surv_clim <- left_join(surv, clim_mat) %>%
                subset( !is.na(location) )


# fit climate models -------------------------------------------------

climate_mods <- list(

  # null 
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_tm1 + (log_area_t0 | year) + (1 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_tm1 + (log_area_t0 | year) + (1 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_tm1 + (log_area_t0 | year) + (1 | location)
  
)

# brier score
brier_score <- function(x){

  sum2 <- test_df %>% 
            mutate( pred = boot::inv.logit(x) ) %>% 
            mutate( sum2 = (resp - pred)^2 ) %>% 
            subset( !is.na(sum2) )
  
  sum2$sum2
  
}

# fit models
mods <- lapply( climate_mods,
                function(x) glmer(x, data=surv_clim, family='binomial') ) %>%
          setNames( c( 'null',
                       'ppt_t0', 'ppt_tm1', 
                       'tmp_t0', 'tmp_tm1',
                       'oni_t0', 'oni_tm1') )

surv_scor <- list()

# fit all models
for(ii in 1:length(years) ){
  
  train_df <- subset(surv_clim, !(year %in% years[ii]) )
  test_df  <- subset(surv_clim,   year == years[ii]) %>% 
                mutate( resp = surv_t1 )
  
  # fit all models
  surv_clim_m <- lapply( climate_mods, 
                         function(x) glmer(x, data=train_df, family='binomial') ) %>% 
                    setNames( c( 'null',
                         'ppt_t0', 'ppt_tm1', 
                         'tmp_t0', 'tmp_tm1',
                         'oni_t0', 'oni_tm1') )
    
  surv_pred       <- lapply(surv_clim_m, predict, newdata=test_df, re.form=NA)
  surv_scor[[ii]] <- lapply(surv_pred, brier_score) %>% as.data.frame
  
}

# get scores
score_df <- Reduce(function(...) rbind(...), surv_scor) %>% 
              colMeans %>% 
              as.data.frame %>% 
              tibble::add_column(model = row.names(.), .before=1) %>% 
              rename( score = "." ) %>% 
              arrange( score )
  
# fit best model with all data
best_mod <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0 + 
                    (log_area_t0 | year) + (1 | location),
                  data = surv_clim, family='binomial')

re_df   <- ranef(best_mod)$year %>% 
              tibble::add_column(.before=1, coef = row.names(.) ) %>% 
              bind_rows( ranef(best_mod)$location %>% 
                         tibble::add_column(.before=1, coef = row.names(.) )
              ) %>% 
              mutate( type_coef = 'ranef' )
  
out_df  <- fixef(best_mod) %>% 
              t %>% t %>% 
              as.data.frame %>% 
              tibble::add_column(.before=1, ranef = row.names(.) ) %>% 
              mutate( type_coef = 'fixef' ) %>% 
              bind_rows( re_df )
     
write.csv(out_df, 
          'results/ml_mod_sel/surv/surv_best_mod_cv.csv',
          row.names=F)
write.csv(score_df, 'results/ml_mod_sel/surv/surv_clim_sel_cv.csv', row.names = F)
