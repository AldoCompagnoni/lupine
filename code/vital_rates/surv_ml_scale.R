rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(mgcv)
library(lme4)
library(bbmle)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")

# data format --------------------------------------------------------------
grow        <- lupine_df %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1    = log(area_t1),
                          log_area_t0    = log(area_t0),
                          log_area_t02   = log(area_t0)^2 ) %>% 
                  mutate( log_area_t0_s  = log_area_t0 %>% scale %>% as.vector,
                          log_area_t02_s = log_area_t02 %>% scale %>% as.vector )


# surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
#                   subset( stage_t0 != 'SL' ) %>%
#                   subset( area_t0 != 0) %>%
#                   mutate( log_area_t0    = log(area_t0),
#                           year           = year + 1 ) %>% 
#                   mutate( log_area_t0_s  = log_area_t0 %>% scale %>% as.vector) %>% 
#                   mutate( log_area_t02   = log_area_t0^2,
#                           log_area_t02_s = log_area_t0_s^2 )


best_mod    <- lmer(log_area_t1 ~ log_area_t0 + log_area_t02     + (log_area_t0 | year) + (1 | location),
                     data = grow)
best_mod_s  <- lmer(log_area_t1 ~ log_area_t0_s + log_area_t02_s + (log_area_t0 | year) + (1 | location),
                     data = grow)



fixef(best_mod_s)[2] / sd(grow$log_area_t0)
fixef(best_mod_s)[3] / sd(grow$log_area_t02)
fixef(best_mod_s)[1] + 1 - sum( fixef(best_mod_s)[2:3] * 
                                  c(mean(grow$log_area_t0),mean(grow$log_area_t02)) ) 

# 
cm  <- grow %>% 
          select(log_area_t0,log_area_t02) %>% 
          colMeans
csd <- grow %>% 
          select(log_area_t0,log_area_t02) %>% 
          apply(2,sd)


# function to rescale coefficients
rescale.coefs <- function(beta, mu, sigma) {
    beta2 <- beta ## inherit names etc.
    beta2[-1] <- sigma[1]*beta[-1]/sigma[-1]
    beta2[1]  <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
    beta2
}  

# seems like a big enough gain?
logLik(best_mod)-logLik(best_mod_s)


fixef(best_mod)
(cc    <- rescale.coefs(fixef(best_mod_s), mu=c(0,cm), sigma=c(1,csd) ))



# representation --------------------------------------------------

# unscaled model
x_seq     <- seq( min(surv$log_area_t0, na.rm=T),
                  max(surv$log_area_t0, na.rm=T),
                  length.out = 100 )

coef_u    <- fixef(best_mod)

y_pred_u  <- coef_u[1] + 
             coef_u[2] * x_seq + 
             coef_u[3] * x_seq^2


# scaled model
x_seq_s   <- seq( min(surv$log_area_t0_s, na.rm=T),
                  max(surv$log_area_t0_s, na.rm=T),
                  length.out = 100 )

coef_s    <- fixef(best_mod_s)

y_pred_s  <- coef_s[1] + 
             coef_s[2] * x_seq_s + 
             coef_s[3] * x_seq_s^2



par( mfrow = c(2,1) )

plot(jitter(surv_t1) ~log_area_t0, data=surv)
lines(x_seq, boot::inv.logit(y_pred_u) )

plot(jitter(surv_t1) ~log_area_t0_s, data=surv)
lines(x_seq_s, boot::inv.logit(y_pred_s) )


fixef(best_mod)
fixef(best_mod_s)






library(ggplot2) 
library(dplyr)

# function to rescale variables around 0
zscore <- function(x) (x - mean(x)) / sd(x)

# creates new dataframe with scaled values included
 %>% %>%  <- diamonds %>% 
      select(carat, cut, price, depth) %>% 
      mutate_all(funs(s = zscore(.)))
