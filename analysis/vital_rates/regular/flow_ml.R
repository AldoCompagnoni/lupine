# simplified script to test for SITE BY CLIMATE interaction
rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(testthat)
library(bbmle)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------
flow <- subset(lupine_df, !is.na(flow_t1) ) %>% 
          subset( area_t1 != 0) %>% 
          mutate( log_area_t1  = log(area_t1),
                  log_area_t12 = log(area_t1)^2,
                  year         = year + 1 ) #%>% 
          #mutate( log_area_t1  = scale(log_area_t1) %>% as.vector,
          #        log_area_t12 = scale(log_area_t12) %>% as.vector )

        
# climate format ----------------------------------------------------------------
years     <- unique(flow$year)
m_obs     <- 5
m_back    <- 36

# format climate - need to select climate predictor first 
tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs)

# demography plus clim
flow_clim  <- left_join(flow, tmp_mat) %>%
                subset( !is.na(location) )


# prepare data
y     <- flow_clim$flow_t1
x1    <- flow_clim %>% 
          select( log_area_t1, V1:V12 ) %>% 
          as.matrix
x2    <- flow_clim %>% 
          select( log_area_t1, V13:V24 ) %>% 
          as.matrix

mod1  <- cv.glmnet( x = x1, y = y,
                    family = 'binomial',n=15,alpha=c(1),keep=T)
mod2  <- cv.glmnet( x = x2, y = y,
                    family = 'binomial',n=15,alpha=c(1),keep=T)


# estimated beta values
beta1 <- coef(mod1, s = mod1$lambda.1se)[,1][-1]
beta2 <- coef(mod2, s = mod2$lambda.1se)[,1][-1]

beta1[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)

beta2[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)


# glmmLasso --------------------------------------------------

# flower data
flow_clim <- flow_clim %>% 
              mutate( year     = as.factor(year),
                      location = as.factor(location) )

# set up lambda values
lam_v <- c(-5:6) %>% exp
aic_l <- list()

# fit models
for( ii in 1:length(lam_v) ){
  
  mod <- 
    glmmLasso( flow_t1 ~ log_area_t1 + 
                         V1 + V2 + V3 + V4 + V5 + V6 + 
                         V7 + V8 + V9 + V10 + V11 + V12, 
               rnd = list( year = ~1 + log_area_t1 ),
               lambda = lam_v[ii],
               family = binomial(link = logit),
               data = flow_clim)
  
  aic_l[[ii]] <- mod$aic
  
}

# final ID
final_i   <- which( (aic_l %>% unlist) == min(aic_l %>% unlist) )

# final model (best model is 9, lambda ~= 20)
final_mod <- glmmLasso( flow_t1 ~ log_area_t1 + 
                          V1 + V2 + V3 + V4 + V5 + V6 + 
                          V7 + V8 + V9 + V10 + V11 + V12, 
                        rnd    = list( year = ~1 + log_area_t1 ),
                        lambda = lam_v[final_i],
                        family = binomial(link = logit),
                        data   = flow_clim )  

# put out the 
tiff( 'results/vital_rates/regular/flow_glmm_lasso.tiff',
      width = 6.3, height = 6.3, unit = 'in', res = 600,
      compression = 'lzw')

coef(final_mod) %>% .[-c(1:4)] %>% plot(type='l', main = 'year one')
abline( h = 0 , lty = 2 )

dev.off()
