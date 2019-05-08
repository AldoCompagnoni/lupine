# PVA of lupine populations
# 1. reate starting abundances for each population 
  # a. Calculate Stable stage distribution for each population (location)
  # b. Redistribute individuals across stages using observed numbers in 2018
# 2. Estimate and create stochastic IPM
# 3. 
rm(list=ls())
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)
library(lintr)
library(goodpractice)


# data
lupine_df   <- read.csv( "data/lupine_all.csv") 
lupine_cnt  <- read.csv( "data/lupine_counts.csv") 
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
pred_g      <- read_xlsx('data/post predation_lupinus tidestromii.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")
germ        <- read_xlsx('data/seedbaskets.xlsx') %>% 
                 select(g0:g2) %>% 
                 colMeans
germ_adj    <- read.csv('results/ml_mod_sel/germ/germ_adj.csv')

# format climate data ----------------------------------------
years     <- c(2005:2018)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var){
  
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


# vital rates format --------------------------------------------------------------

# first, format site/year combinations
site_df     <- select(lupine_df, year, location) %>% 
                 unique %>% 
                 # create 'Site' column to merge with consumption dat
                 mutate( Site = gsub(' \\([0-9]\\)','',location) ) %>% 
                 subset( year > 2004 ) %>% 
                 arrange( location, year ) %>% 
                 complete(location,year)

surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year ) %>% 
                  left_join( clim_mat ) 

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year ) %>% 
                  left_join( clim_mat ) 

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) 

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
                right_join( select(site_df,-Site) ) %>% 
                mutate( ab_r_m = replace(ab_r_m,
                                         is.na(ab_r_m),
                                         mean(ab_r_m, 
                                              na.rm=T)) )

cons_df    <- read_xlsx('data/consumption.xlsx') %>% 
                mutate( Mean_consumption = Mean_consumption %>% as.numeric) %>% 
                select( Year, Site, Mean_consumption) %>% 
                # update contents/names of variables for merging
                mutate( Site = toupper(Site) ) %>% 
                rename( year = Year ) %>% 
                # expand to all site/year combinations
                right_join( site_df ) %>% 
                mutate( Mean_consumption = replace(Mean_consumption,
                                                   is.na(Mean_consumption),
                                                   mean(Mean_consumption,na.rm=T) 
                                                   ) ) %>% 
                # remove NA locations
                subset( !is.na(location) ) %>% 
                # remove annoying code
                select( -Site ) %>% 
                rename( cons = Mean_consumption ) %>% 
                arrange(location,year)

germ_df     <- site_df %>% 
                 select( location ) %>%
                 unique %>% 
                 arrange( location ) %>% 
                 left_join(germ_adj) %>% 
                 # create germ_obs for AL (1)
                 mutate( germ_obs = replace( germ_obs, 
                                             location == 'AL (1)',
                                             mean(germ_obs[location %in% c('ATT (8)',
                                                                           'POP9 (9)')])) ) %>% 
                 # create germ_obs for BR (6)
                 mutate( germ_obs = replace( germ_obs, 
                                             location == 'BR (6)',
                                             mean(germ_obs[location %in% c('BS (7)',
                                                                           'DR (3)')])) ) %>%   
                 # post-dispersal predation
                 mutate( post_d_p = (germ['g0'] - germ_obs) / germ['g0'] ) %>% 
                 # add germination rats
                 mutate( g0 = germ['g0'],
                         g1 = germ['g1'],
                         g2 = germ['g2'] )


# 1. reate SSD for each population -------------------------------

# fit models without year effect
mod_s    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=surv, family='binomial')
mod_g    <- lmer( log_area_t1 ~ log_area_t0 + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=grow)
mod_g2   <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range( c(grow$log_area_t0, grow$log_area_t1) )
mod_fl   <- glmer(flow_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=flow, family='binomial')
mod_fr   <- glmer(numrac_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | location) + (0 + log_area_t0 | location), 
                  data=fert, family='poisson')

# vital rate models 
surv_p    <- fixef(mod_s)
grow_p    <- fixef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- fixef(mod_fl)
fert_p    <- fixef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ * (1 - 0.43) # Old post-dispersal predation estimate



# models ---------------------------------------------------------
mod_s    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=surv, family='binomial')
mod_g    <- lmer( log_area_t1 ~ log_area_t0 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=grow)
mod_g2   <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range( c(grow$log_area_t0, grow$log_area_t1) )
mod_fl   <- glmer(flow_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) +
                  (1 | location) + (0 + log_area_t0 | location),
                  data=flow, family='binomial')
mod_fr   <- glmer(numrac_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location), 
                  data=fert, family='poisson')
# mod_fr   <- MASS::glm.nb(numrac_t0 ~ log_area_t0 , data=fert )
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))

# vital rate models 
surv_p    <- fixef(mod_s)
grow_p    <- fixef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- fixef(mod_fl)
fert_p    <- fixef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ * (1 - 0.43) # Old post-dispersal predation estimate


# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }

# list of mean IPM parameters. 
pars_mean   <- list( # adults vital rates           
                     surv_b0      = surv_p['(Intercept)'],
                     surv_b1      = surv_p['log_area_t0'],
                     surv_b2      = surv_p['log_area_t02'],
                     surv_b3      = surv_p['log_area_t03'],
                     surv_clim    = surv_p['tmp_tm1'],
                     
                     grow_b0      = grow_p['(Intercept)'],
                     grow_b1      = grow_p['log_area_t0'],
                     grow_sig     = grow_p[3],

                     flow_b0      = flow_p['(Intercept)'],
                     flow_b1      = flow_p['log_area_t0'],
                     flow_clim    = flow_p['tmp_tm1'],
                     
                     fert_b0      = fert_p['(Intercept)'],
                     fert_b1      = fert_p['log_area_t0'],
                     fert_clim    = fert_p['tmp_tm1'],
                     
                     abort        = 0.22, # hardcoded for now!
                     clip         = 0.57, # hardcoded for now!
                     
                     fruit_rac    = fr_rac_p,
                     seed_fruit   = seed_fr_p,
                     g0           = germ_p['g0'],
                     g1           = germ_p['g1'],
                     g2           = germ_p['g2'],
                     
                     recr_sz      = size_sl_p$mean_sl_size,
                     recr_sd      = size_sl_p$sd_sl_size,
                     
                     L_sl         = size_sl_p$min_sl_size,
                     U_sl         = size_sl_p$max_sl_size,
                     
                     L            = g_lim[1],
                     U            = g_lim[2],
                     
                     mat_siz_sl = 100,
                     mat_siz    = 100 )

# test that no NAs present
expect_equal(pars_mean %>% 
               unlist %>% 
               is.na() %>% 
               sum, 0)


# update with yearly parameters
update_par <- function(loc_n){
  
  pars_yr <- pars_mean
  
  # get location-specific parameters
  get_loc <- function(mod_obj,loc_n){
    ranef(mod_obj)$location %>% 
      tibble::add_column(.,
                         .before=1,
                         location=row.names(.)) %>% 
      subset( location == loc_n) %>% 
      select(-location)
  }
  
  # get info on raceme "loss"
  get_rac <- function(cons_df, abor_df, germ_df, loc_n){
    
    # raceme consumption
    out_cons <- cons_df %>% 
                  subset(location == loc_n ) %>% 
                  group_by( location ) %>% 
                  summarise( cons = mean(cons) ) %>% 
                  ungroup %>% 
                  .$cons %>% 
                  as.numeric
    
    # raceme abortion
    out_abor <- abor_df %>% 
                  subset( location == loc_n ) %>%
                  group_by( location ) %>% 
                  summarise( ab_r_m = mean(ab_r_m) ) %>% 
                  ungroup %>% 
                  .$ab_r_m %>% 
                  as.numeric
    
    # germination rates
    out_g    <- germ_df %>% 
                  subset( location == loc_n ) %>% 
                  mutate( g0 = g0 * (1-post_d_p),
                          g1 = g1 * (1-post_d_p),
                          g2 = g2 * (1-post_d_p) ) %>% 
                  select(g0,g1,g2)

    # out it all out
    data.frame( clip = out_cons,
                abor = out_abor) %>% bind_cols( out_g )
    
  }
  
  # location information
  pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_loc(mod_s,loc_n)['(Intercept)']) %>% as.numeric
  pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_loc(mod_s,loc_n)['log_area_t0']) %>% as.numeric

  pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_loc(mod_g,loc_n)['(Intercept)']) %>% as.numeric
  pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_loc(mod_g,loc_n)['log_area_t0']) %>% as.numeric
  
  pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_loc(mod_fl,loc_n)['(Intercept)']) %>% as.numeric
  pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_loc(mod_fl,loc_n)['log_area_t0']) %>% as.numeric

  pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_loc(mod_fr,loc_n)['(Intercept)']) %>% as.numeric
  pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_loc(mod_fr,loc_n)['log_area_t0']) %>% as.numeric
  
  # update data on clipped racemes
  pars_yr$clip       <- get_rac(cons_df, abor_df, germ_df, loc_n)$clip
  pars_yr$abort      <- get_rac(cons_df, abor_df, germ_df, loc_n)$abor
  pars_yr$abort      <- get_rac(cons_df, abor_df, germ_df, loc_n)$abor
  pars_yr$g0         <- get_rac(cons_df, abor_df, germ_df, loc_n)$g0
  pars_yr$g1         <- get_rac(cons_df, abor_df, germ_df, loc_n)$g1
  pars_yr$g2         <- get_rac(cons_df, abor_df, germ_df, loc_n)$g2
  
  pars_yr
  
}

# test function
update_par("POP9 (9)")$clip
update_par("ATT (8)")$clip
update_par("AL (1)")$clip



# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# Survival at size x
sx<-function(x,pars,tmp_anom){
  # survival prob. of each x size class 
  inv_logit(pars$surv_b0 + 
            pars$surv_b1 * x + 
            pars$surv_b2 * x^2 +
            pars$surv_b3 * x^3 +
            pars$surv_clim * tmp_anom) 
}

# update kernel functions
grow_sd <- function(x,pars){
  pars$a*(exp(pars$b*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy <- function(y,x,pars){
  # returns a *probability density distribution* for each x value
  dnorm(y,  mean = pars$grow_b0 + pars$grow_b1*x, 
            sd   = pars$grow_sig)
}

# transition: Survival * growth
pxy<-function(y,x,pars,tmp_anom){
  sx(x,pars,tmp_anom) * gxy(y,x,pars) 
}


# production of seeds from x-sized mothers
fx <-function(x,pars,tmp_anom){

  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + 
                         pars$flow_b1*x +
                         pars$flow_clim*tmp_anom ) * 
              exp(       pars$fert_b0 + 
                         pars$fert_b1*x + 
                         pars$fert_clim*tmp_anom )
              
  # viable racs
  viab_rac <- tot_rac * (1- pars$abort) * (1-pars$clip)
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars){
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd )
}

fxy <- function(y,x,pars,tmp_anom){
  fx(x,pars,tmp_anom) * recs(y,pars) 
}


# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(tmp_anom, pars){
 
  # set up IPM domains --------------------------------------------------------
 
  # plants
  n   <- pars$mat_siz
  L   <- pars$L 
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # populate kernel ------------------------------------------------------------
  
  # seeds mini matrix
  s_mat     <- matrix(0,2,2)

  # seeds that enter 2 yr-old seed bank
  plant_s2   <- fx(y,pars,tmp_anom) * pars$g2
  
  # seeds that enter 1 yr-old seed bank
  plant_s1   <- fx(y,pars,tmp_anom) * pars$g1

  # seeds that go directly to seedlings germinate right away 
  Fmat       <- (outer(y,y, fxy, pars, tmp_anom) * pars$g0 * h) 

  # seeds that enter 2 yr-old seed bank
  s_mat[2,1] <- 1

  # recruits from the 1 yr-old seedbank
  s1_rec     <- h * recs(y, pars) 
  
  # recruits from the 2 yr-old seedbank
  s2_rec     <- h * recs(y, pars)
  
  # survival and growth of adult plants
  Tmat       <- outer(y,y,pxy,pars,tmp_anom) * h
  
  # rotate <- function(x) t(apply(x, 2, rev))
  # outer(y,y, fxy, pars, h) %>% t %>% rotate %>% image
  
  small_K    <- Tmat + Fmat
  
  # Assemble the kernel -------------------------------------------------------------
  
  # top 2 vectors
  from_plant <- rbind( rbind( plant_s2, plant_s1),
                       small_K )

  # leftmost vectors
  from_seed  <- rbind( s_mat,
                       cbind(s1_rec, s2_rec) )

  k_yx       <- cbind( from_seed, from_plant )

  return(k_yx)
   
  # tests "integrating' functions ---------------------------------------------------
  
  # s_sl
  # expect_true( ((outer( rep(1,100), y_s, s_sl, pars, h_s) %>% t %>% colSums) > 0.99) %>% all )

  # gxy 
  # expect_true( ((outer(y,y,gxy,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
  # gxy_s: huge unintentional eviction. Why?
  # expect_true( ((outer(y_s,y,gxy_s,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
}

# average kernels for each population 
loc_v       <- lupine_df$location %>% unique
ker_avg_pop <- function(loc_n) kernel(0, update_par(loc_n) )
ker_avg_l   <- lapply(loc_v, ker_avg_pop) %>% setNames( loc_v )

# stable stage distributions
ssd_ker     <- function(x){
  eK <- eigen(x)
  w  <- Re(eK$vectors[,1])
  w/sum(w)
}
  
# list of stable stage distributions 
ssd_l  <- lapply(ker_avg_l, ssd_ker)

# starting population sizes
obs_n  <- lupine_cnt %>% 
            mutate( AL = replace(AL,
                                 year == 2018,
                                 AL[year == 2017]) ) %>% 
            subset( year == 2018 ) %>% 
            select( -AL.REST ) %>% 
            gather( Site, pop_n, AL:POP9 ) %>% 
            left_join( site_df ) %>% 
            select( -Site )

# starting pulation vectors
start_vec <- function(loc_n){
  
  ssd_x    <- ssd_l[loc_n] %>% unlist %>% as.numeric
  
  # redistribute individuals across non-seeds 
  cnt      <- subset(obs_n, 
                     location == loc_n)$pop_n
  
  # proportions of individuals out of the seedbank
  non_sl_p <- ssd_x[-c(1:2)] %>% sum
  
  # seed bank 1 and 2 sensu Dangremond et al. 2010
  sb2      <- (ssd_x[1] * cnt) / non_sl_p
  sb1      <- (ssd_x[2] * cnt) / non_sl_p
  
  # prop individuals in non seedling stages
  non_sl   <- ssd_x[-c(1:2)] * cnt
  
  # redistributed individuals
  c(sb2,sb1,non_sl)
  
}

# Starting vector!
start_v_l <- lapply(loc_v, start_vec) %>% setNames(loc_v)



# 2. Estimate and create stochastic IPM ----------------------


# models 
mod_s    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=surv, family='binomial')
mod_g    <- lmer( log_area_t1 ~ log_area_t0 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=grow)
mod_g2   <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range( c(grow$log_area_t0, grow$log_area_t1) )
mod_fl   <- glmer(flow_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) +
                  (1 | location) + (0 + log_area_t0 | location),
                  data=flow, family='binomial')
mod_fr   <- glmer(numrac_t0 ~ log_area_t0 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location), 
                  data=fert, family='poisson')
# mod_fr   <- MASS::glm.nb(numrac_t0 ~ log_area_t0 , data=fert )
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))

# vital rate models 
surv_p    <- fixef(mod_s)
grow_p    <- fixef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- fixef(mod_fl)
fert_p    <- fixef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ * (1 - 0.43) # Old post-dispersal predation estimate



# update with yearly parameters
update_par <- function(year_n, loc_n){
  
  pars_sim <- pars_mean
  
  # get year-specific parameters
  get_yr <- function(mod_obj,year_n){
  
    ranef(mod_obj)$year %>% 
      tibble::add_column(.,
                         .before=1,
                         year=row.names(.)) %>% 
      subset( year == year_n) %>% 
      select(-year)
  }
  
  # get location-specific parameters
  get_loc <- function(mod_obj,loc_n){
    ranef(mod_obj)$location %>% 
      tibble::add_column(.,
                         .before=1,
                         location=row.names(.)) %>% 
      subset( location == loc_n) %>% 
      select(-location)
  }
  
  # get info on raceme "loss"
  get_rac <- function(cons_df, abor_df, germ_df, 
                      year_n, loc_n){
    
    # raceme consumption
    out_cons <- cons_df %>% 
                  subset(location == loc_n & year == year_n ) %>% 
                  .$cons %>% 
                  as.numeric
    
    # raceme abortion
    out_abor <- abor_df %>% 
                  subset( location == loc_n & year == year_n ) %>% 
                  .$ab_r_m %>% 
                  as.numeric
    
    # germination rates
    out_g    <- germ_df %>% 
                  subset( location == loc_n ) %>% 
                  mutate( g0 = g0 * (1-post_d_p),
                          g1 = g1 * (1-post_d_p),
                          g2 = g2 * (1-post_d_p) ) %>% 
                  select(g0,g1,g2)

    # out it all out
    data.frame( clip = out_cons,
                abor = out_abor) %>% bind_cols( out_g )
    
  }
  
  # year- and site- specific vital rates
  pars_sim$surv_b0    <- (pars_sim$surv_b0 + get_yr(mod_s,year_n)['(Intercept)']) %>% as.numeric
  pars_sim$surv_b1    <- (pars_sim$surv_b1 + get_yr(mod_s,year_n)['log_area_t0']) %>% as.numeric

  pars_sim$grow_b0    <- (pars_sim$grow_b0 + get_yr(mod_g,year_n)['(Intercept)']) %>% as.numeric
  pars_sim$grow_b1    <- (pars_sim$grow_b1 + get_yr(mod_g,year_n)['log_area_t0']) %>% as.numeric
  
  pars_sim$flow_b0    <- (pars_sim$flow_b0 + get_yr(mod_fl,year_n)['(Intercept)']) %>% as.numeric
  pars_sim$flow_b1    <- (pars_sim$flow_b1 + get_yr(mod_fl,year_n)['log_area_t0']) %>% as.numeric

  pars_sim$fert_b0    <- (pars_sim$fert_b0 + get_yr(mod_fr,year_n)['(Intercept)']) %>% as.numeric
  pars_sim$fert_b1    <- (pars_sim$fert_b1 + get_yr(mod_fr,year_n)['log_area_t0']) %>% as.numeric
  
  # locations
  pars_sim$surv_b0    <- (pars_sim$surv_b0 + get_loc(mod_s,loc_n)['(Intercept)']) %>% as.numeric
  pars_sim$surv_b1    <- (pars_sim$surv_b1 + get_loc(mod_s,loc_n)['log_area_t0']) %>% as.numeric

  pars_sim$grow_b0    <- (pars_sim$grow_b0 + get_loc(mod_g,loc_n)['(Intercept)']) %>% as.numeric
  pars_sim$grow_b1    <- (pars_sim$grow_b1 + get_loc(mod_g,loc_n)['log_area_t0']) %>% as.numeric
  
  pars_sim$flow_b0    <- (pars_sim$flow_b0 + get_loc(mod_fl,loc_n)['(Intercept)']) %>% as.numeric
  pars_sim$flow_b1    <- (pars_sim$flow_b1 + get_loc(mod_fl,loc_n)['log_area_t0']) %>% as.numeric

  pars_sim$fert_b0    <- (pars_sim$fert_b0 + get_loc(mod_fr,loc_n)['(Intercept)']) %>% as.numeric
  pars_sim$fert_b1    <- (pars_sim$fert_b1 + get_loc(mod_fr,loc_n)['log_area_t0']) %>% as.numeric

  # update data on clipped racemes
  pars_sim$clip       <- get_rac(cons_df, abor_df, germ_df, 
                                 year_n, loc_n)$clip
  pars_sim$abort      <- get_rac(cons_df, abor_df, germ_df, 
                                 year_n, loc_n)$abor
  pars_sim$abort      <- get_rac(cons_df, abor_df, germ_df, 
                                 year_n, loc_n)$abor
  pars_sim$g0         <- get_rac(cons_df, abor_df, germ_df, 
                                 year_n, loc_n)$g0
  pars_sim$g1         <- get_rac(cons_df, abor_df, germ_df, 
                                 year_n, loc_n)$g1
  pars_sim$g2         <- get_rac(cons_df, abor_df, germ_df, 
                                 year_n, loc_n)$g2
  
  pars_sim
  
}



# stochastic simulations -------------------------------------


# sequence of it years (always the same), and total IPM size
it        <- 50
yrs       <- 2005:2017
loc_v     <- lupine_df$location %>% unique %>% sort
ker_siz   <- (pars_mean$mat_siz)+2
clim_x    <- seq( min(clim_mat$tmp_tm1),
                  max(clim_mat$tmp_tm1), 
                  length.out = 20 )

# store kernels
ker_df    <- expand.grid( yr     = 2005:2017,
                          loc_v  = loc_v,
                          clim_x = clim_x,
                          stringsAsFactors = F )

# store the kernels in advance!
store_ker <- function(ii){
  
  kernel(ker_df$clim_x[ii],
         update_par( ker_df$yr[ii],
                     ker_df$loc_v[ii]) )

}
  
# store all the kernels
ker_l <- lapply(1:nrow(ker_df), store_ker)

# set up 50-year simulations
sim_df    <- expand.grid( sim_i  = 1:1000,
                          clim_x = clim_x,
                          loc_v  = loc_v,
                          stringsAsFactors = F )

# calculate stochastic lambda for each climate anomaly
project_50 <- function(ii){

  # starting population vector
  start_w <- start_v_l[sim_df$loc_v[ii]] %>% 
                unlist %>% 
                as.numeric
  
  # matrix of pop sizes
  mat     <- matrix(NA,length(start_w),50+1)
  
  # random year sequence
  it      <- 50
  yrs     <- 2005:2017
  set.seed(sim_df$sim_i[ii])
  yr_seq  <- sample(yrs,it,replace=T)
  
  # start the series
  mat[,1] <- start_w

  # stochastic simulations
  for(pp in 1:it){
    
    # select kernel: index is from ker_df
    ker_i    <- which( ker_df$yr     == yr_seq[pp] & 
                       ker_df$loc_v  == sim_df$loc_v[ii] & 
                       ker_df$clim_x == sim_df$clim_x[ii] )
    ker      <- ker_l[[ker_i]]
    
    # project population forward
    mat[,pp+1] <- ker %*% mat[,pp]
  
  }

  mat 
  
}

# all population sizes
proj_l <- lapply(1:nrow(sim_df),project_50)

# # get numbers of visible spp
# get_n_visible <- function(ii){
# 
#   as.matrix(proj_l[[ii]][3:100,] %>% colSums) %>% 
#     t %>% 
#     as.data.frame %>% 
#     bind_rows( sim_df[ii,] )
# 
# }
# n_vis_df <- lapply(1:nrow(sim_df), get_n_visible)


# get numbers of visible spp
get_qe <- function(ii, thresh){

  vec_n <- proj_l[[ii]][3:100,] %>% colSums
  any( vec_n < thresh ) %>% as.numeric
  
}

qe_l  <- lapply(1:nrow(sim_df), get_qe, 5)

# calculate quasi extinction probability 
qe_df <- mutate(sim_df, qe = unlist(qe_l) ) %>% 
            group_by( clim_x, loc_v ) %>% 
            summarise( qe_p = sum(qe) ) %>% 
            ungroup %>% 
            mutate( qe_p = qe_p / 1000 )

# plo
pva_df <- qe_df %>% rename( Location = loc_v ) 
ggplot(pva_df) +
  geom_line( aes( x     = clim_x,
                  y     = qe_p,
                  color = Location ), 
             lineend = 'round',
             alpha   = 1,
             size    = 3 ) +
  scale_color_viridis_d() +
  ylab( 'Probability of quasi-extinction' ) +
  xlab( 'Temperature anomaly' ) +
  ggsave('results/ipm/qe_sims/qe_5.tiff',
         width=6.3, height=6.3, compression='lzw')

# store, just in case
save.image(paste0('results/ipm/qe_sims/eq_',
                  thresh,
                  '.RData'))