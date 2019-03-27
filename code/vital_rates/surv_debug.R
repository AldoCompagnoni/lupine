rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(mgcv)
library(lme4)
library(bbmle)
library(ggplot2)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

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
                          log_area_t03 = log_area_t0^3) #%>% 
                  # mutate( log_area_t0  = scale(log_area_t0) %>% as.vector,
                          # log_area_t02 = scale(log_area_t02) %>% as.vector )
  
# demography plus clim
surv <- left_join(surv, clim_mat) %>%
                subset( !is.na(location) )


# fit models --------------------------------------------


# binned survival probabilities
bins <- 15
h    <- (max(surv$log_area_t0,na.rm=T) - min(surv$log_area_t0,na.rm=T)) / bins
lwr  <- min(surv$log_area_t0,na.rm=T) + (h*c(0:(bins-1)))
upr  <- lwr + h
mid  <- lwr + (1/2*h)

binned_prop <- function(lwr_x, upr_x, response){
  
  tmp <- subset(surv, log_area_t0 > lwr_x & log_area_t0 < upr_x ) 
  
  if( response == 'prob' ){   return( sum(tmp$surv_t1) / nrow(tmp) ) }
  if( response == 'n_size' ){ return( nrow(tmp) ) }
  
}

y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
x_binned <- mid
y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist

mod_smpl <- glm( surv_t1 ~ log_area_t0, data = surv, family='binomial' )
mod_smpl <- glm( surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, data = surv, family='binomial' )
mod_gam  <- gam( surv_t1 ~ s(log_area_t0), data = surv, family='binomial' )
coef_smpl<- coef(mod_smpl) 

par(mfrow=c(1,1))
plot( x_binned,y_binned,ylim=c(0,1) )
newd_df <- data.frame(log_area_t0  = x_binned,
                      log_area_t02 = x_binned^2,
                      log_area_t03 = x_binned^3)
lines(x_binned,
      predict.glm(mod_smpl,newdata=newd_df,type='response'))
lines(x_binned, 
      predict(mod_gam,newdata=newd_df,type='response'),
      col = 'blue', lwd=2)

