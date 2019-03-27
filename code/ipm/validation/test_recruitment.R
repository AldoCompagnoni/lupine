rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)
source('analysis/vital_rates/plot_binned_prop.R')


bs <- data.frame( site_id = 'BS (7)',
                  year    = c(2009:2017) )
dr <- data.frame( site_id = 'DR (3)',
                  year    = c(2009:2014,2016:2017) )
nb <- data.frame( site_id = 'NB (2)',
                  year    = c(2010,2012:2016) )
site_all <- list(bs, dr, nb) %>% bind_rows 

# subset by year and site
sub_s_yr    <- function(x_df){
  x_df %>% 
    subset( (location %in% site_all$site_id) & 
             year     %in% site_all$year       ) 
} 


# read data --------------------------------------------

abort_df  <- read.csv( "data/lupine_all.csv") %>% sub_s_yr
fruit_rac <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr <- read_xlsx('data/seedsperfruit.xlsx')
cosump_df <- read_xlsx('data/consumption.xlsx')
germ      <- read_xlsx('data/seedbaskets.xlsx')
lupine_df <- read.csv( "data/lupine_all.csv") %>% sub_s_yr
  
# vr models --------------------------------------------

# germination. g1, g2 not 100% correct, but by 0.01/0.001 points
germ_p    <- select(germ, g0:g2) %>% colMeans 

# n. of fruits per racemes
fr_rac    <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')
fr_rac_p  <- coef(fr_rac) %>% exp

# n. of seeds per fruit
seed_fr   <- glm(SEEDSPERFRUIT ~ 1, data=mutate(seed_x_fr,  
                 # substitute 0 value with really low value (0.01)
                 SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                         SEEDSPERFRUIT == 0, 
                                         0.01) ),
                 family=Gamma(link = "log"))
seed_fr_p <- coef(seed_fr) %>% exp

# consumption 
cons_p    <- cosump_df %>% 
                subset( Site == 'BS') %>% 
                subset( Year %in% c(2009:2017) ) %>% 
                .$Mean_consumption %>% 
                as.numeric %>% 
                mean(na.rm=T)

# abortion 
abort_m   <- glm(cbind(numab_t0, numrac_t0-numab_t0) ~ 1, data = abort_df, family='binomial')
abort_p   <- coef(abort_m) %>% boot::inv.logit()
  

# calculate raw data --------------------------------------------------

# number of seedlings in each site/year combination
sl_n   <- lupine_df %>% 
            # subset( location == site_id ) %>% 
            subset( !is.na(stage_t0) ) %>% 
            subset( stage_t0 == 'SL' ) %>% 
            count( year, location ) %>% 
            mutate( year = year - 1 ) %>% 
            rename( sl_n = n )
   
# total population number             
pop_n   <- lupine_df %>% 
              subset( !is.na(stage_t0) ) %>% 
              count( year, location ) %>% 
              rename( tot_n = n )

# fertility: number of racemes
rac_n   <- lupine_df %>% 
              subset( flow_t0 == 1 ) %>% 
              subset( area_t0 != 0) %>% 
              subset( !is.na(numrac_t0) ) %>% 
              mutate( log_area_t0  = log(area_t0),
                      log_area_t02 = log(area_t0)^2 ) %>% 
              # remove zero fertility (becase fertility should not be 0)
              # NOTE: in many cases, notab_t1 == 0, 
              # because numab_t1 == 0 also
              # subset( !(numrac_t0 %in% 0) ) %>% 
              group_by( year, location ) %>% 
              summarise( tot_rac = sum(numrac_t0) ) %>% 
              ungroup

# check number of recruits
recr_df <- full_join( pop_n, sl_n ) %>% 
              full_join( rac_n ) %>% 
              subset( !is.na(sl_n) ) %>% 
              mutate( seed_n = tot_rac * 
                               (1-(abort_p + cons_p)) * # abortion + clipped
                               fr_rac_p * 
                               seed_fr_p ) %>% 
              mutate( seed_n = round(seed_n,0) ) %>% 
              # WHY ARE NAs HAPPENING?
              subset( !is.na(seed_n) ) %>% 
              # create failed to germinate
              mutate( fail_g = seed_n - sl_n ) %>% 
              # sometimes you have more seedlings than seeds.
              # in this case, assume germination == 100%
              mutate( fail_g = replace(fail_g,
                                       fail_g < 0,
                                       sl_n[fail_g < 0]) )


mod   <- glm(cbind(sl_n, fail_g) ~ 1, data=recr_df, family='binomial')
x_seq <- seq(0,max(recr_df$seed_n), length.out=24)
y_pre <- x_seq * boot::inv.logit( coef(mod) )

# plot it out
ggplot(recr_df, aes(x=seed_n, y=sl_n)) +
  geom_point( ) +
  xlab( expression('Seed number '['t0']) ) +
  ylab( expression('Seedlings '['t1']) ) +
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size = 20) ) +
  geom_line( aes(x = x_seq, y = y_pre),
             size = 1) + 
  ggsave( 'results/validation/seedlings_vs_seeds.tiff',
          width = 6.3, height = 6.3, compression="lzw" )


# # germination rate is abysmally low!
# mod_fl   <- glmer(flow_t0 ~ log_area_t0 + 
#                             (1|year), data=flow, family='binomial')
# mod_fr   <- glmer(numrac_t0 ~ log_area_t0 +
#                               (1|year), data=fert, family='poisson')
# 
# fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
# seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
#                 data=mutate(seed_x_fr,  
#                             # substitute 0 value with really low value (0.01)
#                             SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
#                                                     SEEDSPERFRUIT == 0, 
#                                                     0.01) ),
#                 family=Gamma(link = "log"))
# germ_coef<- select(germ, g0:g2) %>% colMeans
# 
# # --------------------------------------------------
# 
# 
# pars <- list( mat_siz    = 200,
#               U          = g_lim[2],
#               L          = g_lim[1],
#               
#               abort      = 0.22, # hardcoded for now!
#               clip       = 0.57, # hardcoded for now!
#                        
#               fruit_rac  = fr_rac_p,
#               seed_fruit = seed_fr_p,
#               g0         = germ_coef['g0'],
#               g1         = germ_coef['g1'],
#               g2         = germ_coef['g2']
#               
#             )
# 
# vr_yr <- function(ii){
#     pars$flow_b0 = flow_p[ii,'(Intercept)']
#     pars$flow_b1 = flow_p[ii,'log_area_t0']
#     pars$fert_b0 = fert_p[ii,'(Intercept)']
#     pars$fert_b1 = fert_p[ii,'log_area_t0']
#     pars
# }
# 
# 
# inv_logit <- function(x){ exp(x)/(1+exp(x)) } 
# 
# repr <- function(ii){
#   
#   pars<- vr_yr(ii)
#   
#   n   <- pars$mat_siz
#   L   <- pars$L 
#   U   <- pars$U
#   #these are the upper and lower integration limits
#   h   <- (U-L)/n                   #Bin size
#   b   <- L+c(0:n)*h                #Lower boundaries of bins 
#   x   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
#                        
#   # total racemes prod
#   inv_logit( pars$flow_b0 + pars$flow_b1*x ) * 
#   exp(       pars$fert_b0 + pars$fert_b1*x ) 
#   
# }
#   
# year_trans  <- c(2009:2017)
# 
# seed_n <- function(ii){
#   
#   fert %>% 
#     subset( year == year_trans[ii]) %>% 
#     .$numrac_t0 %>% 
#     sum %>% 
#     `*` (1-pars$abort) * 
#     (1-pars$clip) * 
#     pars$fruit_rac * 
#     pars$seed_fruit
#   
# }
# 
# plot(sapply(1:9,seed_n),
#      sl_n$n[-c(1,2)])
# abline(lm(sl_n$n[-c(1,2)] ~ sapply(1:9,seed_n)))
# 
# lm(sl_n$n[-c(1,2)] ~ sapply(1:9,seed_n)) %>% summary
#         
# # viable racs
# viab_rac <- tot_rac 
# # viable seeds
# viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
# viab_sd
