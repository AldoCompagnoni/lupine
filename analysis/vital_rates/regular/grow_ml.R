rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(lme4)
library(bbmle)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------
grow        <- lupine_df %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF", "SL")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 )
                  

# climate format ----------------------------------------------------------------
years     <- unique(grow$year)
m_obs     <- 5
m_back    <- 36

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) 

# demography plus clim
grow_clim  <- left_join(grow, tmp_mat) %>%
                subset( !is.na(location) )


# prepare data
y     <- grow_clim$log_area_t1
x1    <- grow_clim %>% 
          select( log_area_t0, V1:V12 ) %>% 
          as.matrix
x2    <- grow_clim %>% 
          select( log_area_t0, V13:V24 ) %>% 
          as.matrix
x3    <- grow_clim %>% 
          select( log_area_t0, V25:V36 ) %>% 
          as.matrix

mod1  <- cv.glmnet( x = x1, y = y,
                    family = 'gaussian',n=15,alpha=c(1),keep=T)
mod2  <- cv.glmnet( x = x2, y = y,
                    family = 'gaussian',n=15,alpha=c(1),keep=T)
mod3  <- cv.glmnet( x = x3, y = y,
                    family = 'gaussian',n=15,alpha=c(1),keep=T)


# estimated beta values
beta1 <- coef(mod1, s = mod1$lambda.1se)[,1][-1]
beta2 <- coef(mod2, s = mod2$lambda.1se)[,1][-1]
beta3 <- coef(mod3, s = mod2$lambda.1se)[,1][-1]

beta1[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)

beta2[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)

beta3[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)


# glmmLasso --------------------------------------------------


# introduce factors for glmmLasso
grow_clim <- grow_clim %>% 
              mutate( year     = as.factor(year),
                      location = as.factor(location) )

# set up lambda values
lam_v <- c(-5:6) %>% exp
aic_l <- list()

# fit models
for( ii in 1:length(lam_v) ){
  
  mod <- 
    glmmLasso( log_area_t1 ~ log_area_t0 + 
                             V1 + V2 + V3 + V4 + V5 + V6 + 
                             V7 + V8 + V9 + V10 + V11 + V12, 
               rnd = list( year = ~1 + log_area_t0 ),
               lambda = lam_v[ii],
               # family = 'gaussian',
               data = grow_clim)
  
  aic_l[[ii]] <- mod$aic
  
}

# final ID
final_i   <- which( (aic_l %>% unlist) == min(aic_l %>% unlist) )

# final model
final_mod <- mod <- glmmLasso( log_area_t1 ~ log_area_t0 + 
                               V1 + V2 + V3 + V4 + V5 + V6 + 
                               V7 + V8 + V9 + V10 + V11 + V12, 
                             rnd = list( year = ~1 + log_area_t0 ),
                             lambda = lam_v[final_i],
                             data = grow_clim)

# put out the 
tiff( 'results/vital_rates/regular/grow_glmm_lasso.tiff',
      width = 6.3, height = 6.3, unit = 'in', res = 600,
      compression = 'lzw')

coef(final_mod) %>% .[-c(1:4)] %>% plot(type='l', main = 'year one')
abline( h = 0 , lty = 2 )

dev.off()
