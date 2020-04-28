# simplified script to test for SITE BY CLIMATE interaction
rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(mgcv)
library(lme4)
library(bbmle)
library(brms)
library(glmnet)
library(ggplot2)
library(glmmLasso)
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
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3 ) #%>% 
                  # mutate( log_area_t0  = scale(log_area_t0) %>% as.vector,
                  #         log_area_t02 = scale(log_area_t02) %>% as.vector )
  
# climate format ----------------------------------------------------------------
years     <- unique(surv$year)
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

# temperature monthly anomalies
tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) 
  

# demography plus clim
surv_clim <- left_join(surv, tmp_mat) %>%
                subset( !is.na(location) )



# prepare data
y     <- surv_clim$surv_t1
x1    <- surv_clim %>% 
          select( log_area_t0, log_area_t02, log_area_t03, V1:V12 ) %>% 
          as.matrix
x2    <- surv_clim %>% 
          select( log_area_t0, log_area_t02, log_area_t03, V13:V24 ) %>% 
          as.matrix

mod1  <- cv.glmnet( x = x1, y = y,
                    family = 'binomial',n=15,alpha=c(1),keep=T)
mod2  <- cv.glmnet( x = x2, y = y,
                    family = 'binomial',n=15,alpha=c(1),keep=T)


# estimated beta values
beta1 <- coef(mod1, s = mod1$lambda.1se)[,1][-1]
beta2 <- coef(mod2, s = mod2$lambda.1se)[,1][-1]

beta1[-c(1:3)] %>% plot(type='l')
abline(h=0,lty=2)

beta2[-c(1:3)] %>% plot(type='l')
abline(h=0,lty=2)



# glmmLasso --------------------------------------------------

# introduce factors for glmmLasso
surv_clim <- surv_clim %>% 
              mutate( year     = as.factor(year),
                      location = as.factor(location) )

# set up lambda values
lam_v <- c(-5:6) %>% exp
aic_l <- list()

# fit models
for( ii in 1:length(lam_v) ){

  mod <- 
    glmmLasso( surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + 
                 V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, 
               rnd = list( year = ~1 + log_area_t0 ),
               lambda = lam_v[ii],
               family = binomial(link = logit),
               data = surv_clim)
  
  aic_l[[ii]] <- mod$aic

}

# final ID
final_i   <- which( (aic_l %>% unlist) == min(aic_l %>% unlist) )

# final model (best model is 9, lambda ~= 20)
final_mod <- glmmLasso( surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + 
                                   V1 + V2 + V3 + V4 + V5 + V6 + 
                                   V7 + V8 + V9 + V10 + V11 + V12, 
                        rnd    = list( year = ~1 + log_area_t0 ),
                        lambda = lam_v[final_i],
                        family = binomial(link = logit),
                        data   = surv_clim )

# put out the 
tiff( 'results/vital_rates/regular/surv_glmm_lasso.tiff',
      width = 6.3, height = 6.3, unit = 'in', res = 600,
      compression = 'lzw')

coef(final_mod) %>% .[-c(1:4)] %>% plot(type='l', main = 'year one')
abline( h = 0 , lty = 2 )

dev.off()
