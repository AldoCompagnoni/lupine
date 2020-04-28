rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(testthat)
library(bbmle)
library(lme4)
library(readxl)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')



# data format --------------------------------------------------------------
fert        <- subset(lupine_df, flow_t1 == 1 ) %>% 
                  subset( area_t1 != 0) %>% 
                  subset( !is.na(numrac_t1) ) %>% 
                  # remove 
                  subset( !(flow_t1 %in% 0) ) %>% 
                  mutate( log_area_t1  = log(area_t1),
                  log_area_t12 = log(area_t1)^2,
                  year         = year + 1 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t1 %in% 0) ) #%>% 
                  # mutate( log_area_t1  = scale(log_area_t1) %>% as.vector,
                  #         log_area_t12 = scale(log_area_t12) %>% as.vector )
        

# climate format ----------------------------------------------------------------
years     <- unique(fert$year)
m_obs     <- 5
m_back    <- 36

# format climate - need to select climate predictor first 
tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) 

# demography plus clim
fert_clim  <- left_join(fert, tmp_mat) %>%
                subset( !is.na(location) )

# prepare data
y     <- fert_clim$numrac_t1
x1    <- fert_clim %>% 
          select( log_area_t1, V1:V12 ) %>% 
          as.matrix
x2    <- fert_clim %>% 
          select( log_area_t1, V13:V24 ) %>% 
          as.matrix

mod1  <- cv.glmnet( x = x1, y = y,
                    family = 'poisson',n=15,alpha=c(1),keep=T)
mod2  <- cv.glmnet( x = x2, y = y,
                    family = 'poisson',n=15,alpha=c(1),keep=T)


# estimated beta values
beta1 <- coef(mod1, s = mod1$lambda.1se)[,1][-1]
beta2 <- coef(mod2, s = mod2$lambda.1se)[,1][-1]

beta1[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)

beta2[-c(1)] %>% plot(type='l')
abline(h=0,lty=2)



# glmmLasso --------------------------------------------------


# introduce factors for glmmLasso
fert_clim <- fert_clim %>% 
              mutate( year     = as.factor(year),
                      location = as.factor(location) )

# set up lambda values
lam_v <- c(-5:6) %>% exp
aic_l <- list()

# fit models
for( ii in 1:length(lam_v) ){
  
  mod <- 
    glmmLasso( numrac_t1 ~ log_area_t1 + 
                 V1 + V2 + V3 + V4 + V5 + V6 + 
                 V7 + V8 + V9 + V10 + V11 + V12, 
               rnd = list( year = ~1 + log_area_t1 ),
               lambda = lam_v[ii],
               family = poisson(link='log'),
               data = fert_clim)
  
  aic_l[[ii]] <- mod$aic
  
}

# final ID
final_i   <- which( (aic_l %>% unlist) == min(aic_l %>% unlist) )

# final model (best model is 9)
final_mod <- glmmLasso( numrac_t1 ~ log_area_t1 + 
                           V1 + V2 + V3 + V4 + V5 + V6 + 
                           V7 + V8 + V9 + V10 + V11 + V12, 
                         rnd = list( year = ~1 + log_area_t1 ),
                         lambda = lam_v[final_i],
                         family = poisson(link='log'),
                         data = fert_clim)

# put out the 
tiff( 'results/vital_rates/regular/fert_glmm_lasso.tiff',
      width = 6.3, height = 6.3, unit = 'in', res = 600,
      compression = 'lzw')

coef(final_mod) %>% .[-c(1:4)] %>% plot(type='l', main = 'year one')
abline( h = 0 , lty = 2 )

dev.off()


select(tmp_mat, V1:V12) %>% 
  cor %>% 
  .[,'V5'] %>% 
  .[c(1:4,6:12)] %>% 
  plot(type='l')

abline(h=0, lty=2)
