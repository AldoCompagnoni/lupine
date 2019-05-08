rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)
source('analysis/vital_rates/plot_binned_prop.R')


att<- data.frame( location = 'ATT (8)',
                  year     = c(2008:2012) )
p9 <- data.frame( location = 'POP9 (9)',
                  year     = c(2008:2014) )
bs <- data.frame( location = 'BS (7)',
                  year     = c(2009:2017) )
dr <- data.frame( location = 'DR (3)',
                  year     = c(2009:2014,2016:2017) )
nb <- data.frame( location = 'NB (2)',
                  year     = c(2010,2012:2016) )
site_all <- list(bs, dr, nb, att, p9) %>% bind_rows 

# subset by year and site
sub_s_yr    <- function(x_df){
  x_df %>% 
    subset( (location %in% site_all$location) & 
             year     %in% site_all$year       ) 
}

# read data --------------------------------------------
fruit_rac <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr <- read_xlsx('data/seedsperfruit.xlsx')
germ      <- read_xlsx('data/seedbaskets.xlsx')
lupine_df <- read.csv( "data/lupine_all.csv") %>% sub_s_yr
site_df   <- select(lupine_df,location) %>% 
                unique %>% 
                mutate( Site = gsub(' \\([0-9]\\)','', location) %>% 
                                  toupper ) 
  

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
cons_df    <- read_xlsx('data/consumption.xlsx') %>% 
                mutate( Mean_consumption = Mean_consumption %>% as.numeric) %>% 
                select( Year, Site, Mean_consumption) %>% 
                rename( cons = Mean_consumption,
                        year = Year ) %>% 
                # expand potential "cases"
                complete( Site, year) %>% 
                # update name
                mutate( Site = toupper(Site) ) %>% 
                mutate( cons = replace(cons,
                                       is.na(cons),
                                       mean(cons[Site!='AL'],na.rm=T)
                                       ) ) %>% 
                left_join( site_df ) %>% 
                # remove NA locations
                subset( !is.na(location) ) %>% 
                # remove annoying code
                select( -Site )
                
# abortion 
abor_df  <- subset(lupine_df, !is.na(flow_t0) & flow_t0 == 1 ) %>% 
                subset( !is.na(numrac_t0) ) %>% 
                # remove non-flowering individuals
                subset( !(flow_t0 %in% 0) ) %>% 
                # remove zero fertility (becase fertility should not be 0)
                subset( !(numrac_t0 %in% 0) ) %>% 
                # only years indicated by Tiffany
                subset( year %in% c(2010, 2011, 2013:2017) ) %>% 
                # calculate abortion rates
                mutate( ab_r      = numab_t0 / numrac_t0 ) %>% 
                group_by( location, year ) %>% 
                summarise( ab_r_m = mean(ab_r, na.rm=T) ) %>% 
                ungroup %>% 
                right_join( site_all ) %>% 
                mutate( ab_r_m = replace(ab_r_m,
                                         is.na(ab_r_m),
                                         mean(ab_r_m, 
                                              na.rm=T)) )
  
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
              ungroup %>% 
              left_join( abor_df ) %>% 
              left_join( cons_df ) %>% 
              inner_join( site_all )

# check number of recruits
recr_df <- full_join( pop_n, sl_n ) %>% 
              full_join( rac_n ) %>% 
              subset( !is.na(sl_n) ) %>% 
              # select only year/sites that we need
              inner_join( site_all ) %>% 
              mutate( seed_n = tot_rac * 
                               (1-ab_r_m) * (1-cons) * # abortion + clipped
                               fr_rac_p * 
                               seed_fr_p ) %>% 
              # round to fit a binomial model
              mutate( seed_n = round(seed_n,0) ) %>% 
              # create failed to germinate
              mutate( fail_g = seed_n - sl_n )


# model selection 
mod      <- glm(cbind(sl_n, fail_g) ~ 1, data=recr_df, family='binomial')
mod_loc  <- glm(cbind(sl_n, fail_g) ~ location, data=recr_df, family='binomial')
mod_yr   <- glm(cbind(sl_n, fail_g) ~ year, data=recr_df, family='binomial')
mod_locyr<- glm(cbind(sl_n, fail_g) ~ location:year, data=recr_df, family='binomial')
AIC(mod, mod_loc, mod_yr, mod_locyr)

# prepare prediction data frame
pred_df  <- recr_df %>% select(location) %>% unique
pred_v   <- predict(mod_loc, 
                    type    = 'response',
                    newdata = pred_df )
pred_df  <- mutate( pred_df, 
                    coef = pred_v ) %>% 
            mutate( x = rep(2000,5),
                    y = c(340,320,300,280,260) ) %>% 
            mutate( lab = paste0(location,'=',
                                 round(coef,4)*100,'%')) 


# plot it out
ggplot(recr_df, aes(x     = seed_n, 
                    y     = sl_n,
                    color = location) ) +
  geom_point( ) +
  xlab( expression('Seed number'['t0']) ) +
  ylab( expression('Seedlings'['t1']) ) +
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size = 20) ) +
  geom_abline( data = pred_df,
               aes(intercept = 0, 
                   slope     = coef,
                   color     = location),
               alpha = 0.7,
               size  = 1 ) + 
  # plot percent germination
  geom_text( data= pred_df,
                aes(x=x,
                y=y,
                label = lab),
             vjust = 1) +
  scale_color_viridis_d() +
  ggsave( 'results/ipm/validation/seedlings_vs_seeds.tiff',
          width = 6.3, height = 6.3, compression="lzw" )

# store germination adjustement parameter
pred_df %>% 
  mutate( coef = pred_v ) %>% 
  write.csv('results/ml_mod_sel/germ/germ_adj.csv',row.names=F)


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
