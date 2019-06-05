# Stochastic LTRE: lam_s ~ temp. anomaly
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
                     
                     L            = g_lim[1],
                     U            = g_lim[2],
                     
                     mat_siz    = 100 )

# test that no NAs present
expect_equal(pars_mean %>% 
               unlist %>% 
               is.na() %>% 
               sum, 0)


# update with yearly parameters
update_par <- function(year_n, loc_n, vr){
  
  pars_yr <- pars_mean
  
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
  
  if( vr == 'surv'){
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_yr(mod_s,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_yr(mod_s,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_loc(mod_s,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_loc(mod_s,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'grow'){
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_yr(mod_g,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_yr(mod_g,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_loc(mod_g,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_loc(mod_g,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'flow'){
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_yr(mod_fl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_yr(mod_fl,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_loc(mod_fl,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_loc(mod_fl,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'fert'){
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_yr(mod_fr,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_yr(mod_fr,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_loc(mod_fr,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_loc(mod_fr,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'all' ){
    # years
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_yr(mod_s,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_yr(mod_s,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_yr(mod_g,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_yr(mod_g,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_yr(mod_fl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_yr(mod_fl,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_yr(mod_fr,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_yr(mod_fr,year_n)['log_area_t0']) %>% as.numeric
    
    # locations
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_loc(mod_s,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_loc(mod_s,loc_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_loc(mod_g,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_loc(mod_g,loc_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_loc(mod_fl,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_loc(mod_fl,loc_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_loc(mod_fr,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_loc(mod_fr,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  # update data on clipped racemes
  pars_yr$clip       <- get_rac(cons_df, abor_df, germ_df, 
                                year_n, loc_n)$clip
  pars_yr$abort      <- get_rac(cons_df, abor_df, germ_df, 
                                year_n, loc_n)$abor
  pars_yr$abort      <- get_rac(cons_df, abor_df, germ_df, 
                                year_n, loc_n)$abor
  pars_yr$g0         <- get_rac(cons_df, abor_df, germ_df, 
                                year_n, loc_n)$g0
  pars_yr$g1         <- get_rac(cons_df, abor_df, germ_df, 
                                year_n, loc_n)$g1
  pars_yr$g2         <- get_rac(cons_df, abor_df, germ_df, 
                                year_n, loc_n)$g2
  
  pars_yr
  
}

# test function
update_par(2010,"POP9 (9)",'all')$abor
update_par(2010,"ATT (8)",'all')$g0
update_par(2010,"POP9 (9)",'all')$g2

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

ker <- kernel(0,pars_mean)
Re(eigen(ker)$value[1])


# stochastic simulations -------------------------------------

# sequence of it years (always the same), and total IPM size
it        <- 50000
yrs       <- 2005:2017
loc_v     <- lupine_df$location %>% unique %>% sort
set.seed(1776)
yr_seq    <- sample(yrs,it,replace=T)
yr_i      <- yr_seq - 2004
ker_siz   <- (pars_mean$mat_siz)+2


# calculate stochastic lambda for each climate anomaly
lam_stoch <- function(tmp_anom, loc, vr){

  # store kernels
  ker_l     <- list()
  for(ii in seq_along(yrs)){
    ker_l[[ii]] <- kernel(tmp_anom,update_par(yrs[ii],loc, vr))
  }

  # placeholders
  n0        <- rep(1/ker_siz,ker_siz)
  r_v       <- rep(0,it)

  # stochastic simulations
  for(ii in 1:it){
    
    #Store kernel
    ker     <- ker_l[[yr_i[ii]]]
    # ker     <- kernel(tmp_anom,update_par(yrs[yr_i[ii]]))
    
    # expect_true( all.equal(ker1,ker) )
    
    # calculate lambdas
    n0      <- ker %*% n0
    N       <- sum(n0)
    r_v[ii] <- log(N)
    n0      <- n0/N
  
  }

  r_v[-c(1:1000)] %>% mean %>% exp
  
}

# test stochasticity
lamS_noClim <- lam_stoch(0,loc_v[3],'all')

# set up climate covariates
clim_x <- seq( min(clim_mat$tmp_tm1),
               max(clim_mat$tmp_tm1), 
               length.out = 10 )

lam_s_l <- list()

# cycle through the different locations
for(li in 1:length(loc_v)){ #
  lam_s_vec      <- sapply(clim_x, lam_stoch,loc_v[li], 'all')
  lam_s_l[[li]] <- data.frame( climate_anomaly   = clim_x,
                                lam_s            = lam_s_vec,
                                location         = loc_v[li],
                                stringsAsFactors = F)
}

#4,5,7
lam_s_df <- bind_rows( lam_s_l ) %>% 
              mutate( location = as.factor(location) )


# lambdaS_vs_clim
ggplot(lam_s_df, aes(climate_anomaly, lam_s) ) +
  geom_line( aes(color = location ),
             size = 1 ) + 
  viridis::scale_color_viridis( discrete=T ) + 
  geom_hline( yintercept = 1,
              linetype   = 2 ) + 
  ylab( expression(lambda['s']) ) + 
  xlab( 'Climate anomaly' ) +
  theme( axis.title = element_text( size = 30) ) + 
  ggsave('results/ipm/lambdaS_vs_site_clim_noseedl.tiff',
          width = 6.3, height = 5, compression = 'lzw' )



# stochastic both climate and random effects simultaneously -----

# set up a seed that produces a mean ~ 0 (but negative, to be conservative)
set.seed(1835)
clim_v <- rnorm(100,0,6.5)
seq_df <- expand.grid( clim_v = clim_v,
                       yrs    = c(2009:2018) )
it     <- 50000
set.seed(1776)
yr_seq    <- sample(1:nrow(seq_df),it,replace=T)

# simulate stochastic lambda
lam_stoch_clim <- function(vr){

  n0        <- rep(1/ker_siz,ker_siz)
  r_v       <- rep(0,it)

  # store kernels
  ker_l     <- list()
  for(ii in 1:nrow(seq_df)){
    ker_l[[ii]] <- kernel(seq_df[ii,]$clim_v,
                          update_par_spec(seq_df[ii,]$yrs,vr) )
  }
  
  # stochastic simulations
  for(ii in 1:it){
    
    #Store kernel
    ker     <- ker_l[[yr_seq[ii]]]
    
    # calculate lambdas
    n0      <- ker %*% n0
    N       <- sum(n0)
    r_v[ii] <- log(N)
    n0      <- n0/N
  
  }

  r_v[-c(1:1000)] %>% mean %>% exp
  
}

lamS_clim <- lam_stoch_clim('all')


# stochastic autocorrelated ---------------------------------------

# re-format temperature data, and estimate var-covar
years   <- c(1990:2018)
tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              year_anom('tmp') %>% 
              select(year,tmp_tm1)
# calculate empirically observed var-covar matrix
tmp_df <- cbind(tmp_mat[1:28,2],tmp_mat[2:29,2])
Sig    <- cov(tmp_df)

# simulate correlated random process
sim_tmp<- MASS::mvrnorm((50000/2), mu = c(0,0), Sigma = Sig)

# concatenate observations by row 
get_row <- function(ii) sim_tmp[ii,1:2]
vec_xx  <- sapply(1:nrow(sim_tmp), get_row) %>% as.numeric    

# do simulations selecting first 60 
yrs    <- 2009:2018
it     <- 50000
set.seed(1776)
yr_seq <- sample(1:length(yrs),it,replace=T)


# simulate stochastic lambda
lam_s_clim_auto <- function(vr){

  n0        <- rep(1/ker_siz,ker_siz)
  r_v       <- rep(0,it)

  # stochastic simulations
  for(ii in 1:it){
    
    #Store kernel
    ker     <- kernel(vec_xx[ii], 
                      update_par_spec(yrs[yr_seq[ii]],vr) )

    # calculate lambdas
    n0      <- ker %*% n0
    N       <- sum(n0)
    r_v[ii] <- log(N)
    n0      <- n0/N
  
  }

  r_v[-c(1:1000)] %>% mean %>% exp
  
}

lam_s_auto  <- lam_s_clim_auto('all')


# put it all together -----------------------------------------

mean_l <- read.csv('results/lambda_mean.csv')

lam_df <- data.frame( lambda = c('determinitic','no_clim','clim', 'clim_auto',
                                 'clim_surv_sl',
                                 'clim_surv', 'clim_grow', 
                                 'clim_flow', 'clim_fert'),
                      value  = c(mean_l$mean_lambda,lamS_noClim,lamS_clim,lam_s_auto,
                                 lamS_surv_sl,
                                 lamS_surv,lamS_grow,lamS_flow,
                                 lamS_fert) ) %>% 
                      mutate( lambda = factor(lambda, 
                                              levels = c('determinitic','no_clim','clim', 'clim_auto',
                                                         'clim_surv_sl',
                                                         'clim_surv', 'clim_grow', 
                                                         'clim_flow', 'clim_fert'))
                              )

write.csv( lam_df,
           'results/lambda_S.csv',
           row.names=F)

# plot stochastic lambdas
ggplot(lam_df[1:4,]) +
  geom_point(aes(x=lambda,y=value,size=5) ) +
  theme( axis.text.x = element_text( angle=90,
                                     size =20),
         axis.title  = element_text( size =30), 
         legend.position ="none" ) + 
  ylab( expression(lambda[s]) ) +
  xlab( expression('Type of '*lambda) ) + 
  ggsave('results/ipm/lambdas.tiff',
         width=6.3, height=6.3, compression='lzw')

# plot single components
ggplot(lam_df[c(2,5:9),]) +
  geom_point(aes(x=lambda,y=value) ) +
  theme( axis.text.x = element_text( angle=90) ) + 
  ylab( expression(lambda[s]) ) +
  xlab( 'simulation' ) + 
  ggsave('results/ipm/lambda_components.tiff',
         width=6.3, height=6.3, compression='lzw')




# set up a seed that produces a mean ~ 0 (but negative, to be conservative)
clim_means  <- lam_s_df$climate_anomaly
lamS_clim_l <- rep(NA,length(clim_means))

for(ii in 1:length(clim_means)){
  set.seed(1835)
  clim_v <- rnorm(100,clim_means[ii],6.5)
  seq_df <- expand.grid( clim_v = clim_v,
                         yrs    = c(2009:2018) )
  it     <- 50000
  set.seed(1776)
  yr_seq    <- sample(1:nrow(seq_df),it,replace=T)
  
  # simulate stochastic lambda
  lam_stoch_clim <- function(vr){
  
    n0        <- rep(1/ker_siz,ker_siz)
    r_v       <- rep(0,it)
  
    # store kernels
    ker_l     <- list()
    for(ii in 1:nrow(seq_df)){
      ker_l[[ii]] <- kernel(seq_df[ii,]$clim_v,
                            update_par_spec(seq_df[ii,]$yrs,vr) )
    }
    
    # stochastic simulations
    for(ii in 1:it){
      
      #Store kernel
      ker     <- ker_l[[yr_seq[ii]]]
      
      # calculate lambdas
      n0      <- ker %*% n0
      N       <- sum(n0)
      r_v[ii] <- log(N)
      n0      <- n0/N
    
    }
  
    r_v[-c(1:1000)] %>% mean %>% exp
    
  }
  
  lamS_clim_l[ii] <- lam_stoch_clim('all')

}


# Compare two levels of stochastic simulations
lam_s_df %>% 
  mutate( lam_s_ran = lamS_clim_l ) %>% 
  gather(lambda_type, lambda, lam_s:lam_s_ran) %>%
  mutate( lambda_type = replace(lambda_type,
                                lambda_type == 'lam_s',
                                'fixed climate') ) %>% 
  mutate( lambda_type = replace(lambda_type,
                                lambda_type == 'lam_s_ran',
                                'varying climate') ) %>% 
  ggplot( aes(x      = climate_anomaly,
              y      = lambda,
              colour = lambda_type) ) +
  viridis::scale_color_viridis( discrete = T ) +
  geom_line( size = 1) +
  geom_hline( yintercept = 1,
              lty=20 ) +
  theme( axis.title  = element_text(size=20),
         legend.text  = element_text(size=20),
         legend.title = element_text(size=20),
         axis.text    = element_text(size=20) ) +
  ylab( expression(lambda[s])) +
  ggsave('results/ipm/lambdaS_vs_clim_ran.tiff',
         width=7,height=4.5,compression='lzw')



# deterministic climate effects -----------------------------------------
clim_x <- seq( min(clim_mat$tmp_tm1),
               max(clim_mat$tmp_tm1), 
               length.out = 100 )

clim_lambda <- function(xx,pars_mean){
  ker <- kernel(xx,pars_mean)
  Re(eigen(ker)$value[1])
}

lam_vs_clim <- sapply(clim_x, clim_lambda, pars_mean)
plot(clim_x, lam_vs_clim, type='l')
abline(h=1,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
abline(h=1.2,lty=2,lwd=2)

sapply(0, clim_lambda, pars_mean)


df <- data.frame( x=clim_x, y=lam_vs_clim)

ggplot(df, aes( x, y) ) +
  geom_line( size = 1 ) + 
  geom_hline( yintercept = 1,
              linetype   = 2 ) + 
  ylab( expression(lambda) ) + 
  xlab( 'Climate anomaly' ) +
  theme( axis.title = element_text( size = 30) ) + 
  ggsave( 'results/lambda_vs_clim.tiff',
          width = 6.3, height = 6.3, compression="lzw" )
