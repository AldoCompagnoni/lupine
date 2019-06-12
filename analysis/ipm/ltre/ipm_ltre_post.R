# Stochastic LTRE: lam_s_diff ~ temp. anomaly
# If lam_s_0: stochatic lambda at temp. anomaly == 0 
# lam_s_temp: stochatic lambda at temp. anomaly != 0 
# lam_s_diff = lam_s_temp - lam_s_0
# do this for every vital rate
rm(list=ls())
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)


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
fit         <- readRDS('results/vital_rates/bayes/lupine_allvr_annual_anom.rds')


# format climate data ----------------------------------------
years     <- 1990:2018
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(clim_x, clim_var = "ppt", 
                      years, m_back, m_obs ){
  
  # "spread" the 12 months
  clim_m <- select(clim_x, -clim_var )
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(years, clim_back, m_obs, clim_m)
  x_clim    <- year_by_month_mat(clim_x_l, years) %>% 
                  gather(month,t0,V1:V12) %>% 
                  select(year,month,t0) %>% 
                  mutate( month = gsub('V','',month) ) %>%
                  mutate( month = as.numeric(month) )
  
  if( clim_var == 'ppt'){
    raw_df <- x_clim %>% 
                group_by(year) %>% 
                summarise( ppt_t0 = sum(t0) ) %>% 
                ungroup %>% 
                arrange( year ) %>% 
                mutate( ppt_t0  = scale(ppt_t0)[,1] ) %>% 
                mutate( ppt_tm1 = lag(ppt_t0) )
  }
  if( clim_var == 'tmp'){
    raw_df <- x_clim %>% 
                group_by(year) %>% 
                summarise( tmp_t0 = mean(t0) ) %>% 
                ungroup %>% 
                arrange( year ) %>% 
                mutate( tmp_t0  = scale(tmp_t0)[,1] ) %>% 
                mutate( tmp_tm1 = lag(tmp_t0) )
  }
  
  raw_df
  
}

# format climate - need to select climate predictor first 
ppt_mat <- subset(clim, clim_var == "ppt") %>%
              year_anom("ppt", years, m_back, m_obs)

tmp_mat <- subset(clim, clim_var == 'tmean')  %>% 
              year_anom("tmp", years, m_back, m_obs) 
  
# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat, tmp_mat) )


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
all_vr   <- rstan::extract(fit) %>% as.data.frame
mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow)
g_lim    <- range( c(grow$log_area_t0, grow$log_area_t1) )
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))

# vital rate models 
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ * (1 - 0.43) # Old post-dispersal predation estimate



# IPM parameters -------------------------------------------------------------

# factor of all the locations
loc_factors <- lupine_df$location %>% unique %>% as.factor
yr_factors  <- lupine_df$year %>% unique %>% as.factor
thin_seq    <- seq(1,6000, length.out=100) %>% round(0)

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }


# list of mean IPM parameters. 
pars_mean   <- list( # adults vital rates           
                     surv_b0      = mean(all_vr$a_u_s)*2,
                     surv_b1      = mean(all_vr$b_u_s)*2,
                     surv_b2      = mean(all_vr$b_s2),
                     surv_b3      = mean(all_vr$b_s3),
                     surv_clim    = 0,
                     
                     grow_b0      = mean(all_vr$a_u_g)*2,
                     grow_b1      = mean(all_vr$b_u_g)*2,
                     grow_sig     = summary(mod_g)$sigma,

                     flow_b0      = mean(all_vr$a_u_f)*2,
                     flow_b1      = mean(all_vr$b_u_f)*2,
                     flow_clim    = 0,
                     
                     fert_b0      = mean(all_vr$a_u_r)*2,
                     fert_b1      = mean(all_vr$b_u_r)*2,
                     fert_clim    = 0,
                     
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
                     
                     mat_siz    = 100 )

# test that no NAs present
expect_equal(pars_mean %>% 
               unlist %>% 
               is.na() %>% 
               sum, 0)

# FOR posterior SIMULATIONS: update with year- and loc-specific parameters 
update_par_post <- function(year_n, loc_n, ii, vr){
  
  pars_sim <- pars_mean
  
  # get year-specific parameters
  get_yr <- function(vr, year_n, ii){
  
    thin_ii <- thin_seq[ii]
    year_ii <- year_n - 2004
    var_n   <- paste0(vr, '.', year_ii)
    all_vr[var_n][thin_ii,1]
  
  }
  
  # get location-specific parameters
  get_loc <- function(vr, loc_n, ii){
    
    thin_ii <- thin_seq[ii]
    loc_ii  <- which(loc_n == loc_factors)
    var_n   <- paste0(vr, '.', loc_ii)
    all_vr[var_n][thin_ii,1]
    
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
  pars_sim$surv_b0    <- get_yr('a_yr_s',year_n,ii) + get_loc('a_loc_s',loc_n,ii)
  pars_sim$surv_b1    <- get_yr('b_yr_s',year_n,ii) + get_loc('b_loc_s',loc_n,ii)
  
  pars_sim$grow_b0    <- get_yr('a_yr_g',year_n,ii) + get_loc('a_loc_g',loc_n,ii)
  pars_sim$grow_b1    <- get_yr('b_yr_g',year_n,ii) + get_loc('b_loc_g',loc_n,ii)
  
  pars_sim$flow_b0    <- get_yr('a_yr_f',year_n,ii) + get_loc('a_loc_f',loc_n,ii)
  pars_sim$flow_b1    <- get_yr('b_yr_f',year_n,ii) + get_loc('b_loc_f',loc_n,ii)
  
  pars_sim$fert_b0    <- get_yr('a_yr_r',year_n,ii) + get_loc('a_loc_r',loc_n,ii)
  pars_sim$fert_b1    <- get_yr('b_yr_r',year_n,ii) + get_loc('b_loc_r',loc_n,ii)
  
  # 
  if( vr == 'surv' ){
    pars_sim$surv_clim  <- all_vr['b_c_s'][thin_seq[ii],1]
  }
  if( vr == 'flow' ){
    pars_sim$flow_clim  <- all_vr['b_c_f'][thin_seq[ii],1]
  }
  if( vr == 'fert' ){
    pars_sim$fert_clim  <- all_vr['b_c_r'][thin_seq[ii],1]
  }
  
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

# test function
update_par_post(2010,"POP9 (9)",3,'surv')$surv_clim
update_par_post(2010,"POP9 (9)",10,'fert')$fert_clim

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
lam_stoch <- function(sim_ii, loc_n, tmp_anom, vr){

  # store kernels
  ker_l     <- list()
  for(ii in seq_along(yrs)){
    ker_l[[ii]] <- kernel(tmp_anom,update_par_post(yrs[ii], loc_n, sim_ii, vr))
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

# calculate and store stochastic lambdas at climate anomaly == 0
base_lam_s   <- sapply(1:100, lam_stoch, loc_v[5], 0, 'all')
lam_base_df  <- data.frame( sim_ii     = 1:100,
                            lam_s_base = base_lam_s )


# simulate stochastic lambdas by location, vital rate, and climate anomaly

# climate anomalies
clim_x <- seq( min(clim_mat$tmp_t0),
               max(clim_mat$tmp_t0), 
               length.out = 10 )

# Data frame for simulations   
sim_df <- expand.grid( clim_x   = clim_x,
                       sim_ii   = 1:100,
                       loc_n    = loc_v[5],
                       vr       = 'surv', #c('surv','flow','fert'),
                       stringsAsFactors = F ) 

# run the simulations in a smarter way
stoch_sim_post <- function(ii){
  
  lam_s <- lam_stoch(sim_df$sim_ii[ii],
                     sim_df$loc_n[ii],
                     sim_df$clim_x[ii],
                     sim_df$vr[ii])
  
  mutate(sim_df[ii,], 
         lam = lam_s)
  
}

lter_l  <- lapply(1:nrow(sim_df), stoch_sim_post)
lter_df <- bind_rows( lter_l )


# compute sensitivity
sens_df <- full_join( lam_base_df, lter_df ) %>% 
              mutate( sens = lam - lam_s_base )

plot_sens_df <- sens_df %>% 
                  group_by(clim_x) %>% 
                  summarise( lam_min  = min(sens),
                             lam_max  = max(sens),
                             lam_mean = mean(sens) )

# 
ggplot(plot_sens_df) +
  geom_ribbon( aes( x    = clim_x, 
                    ymin = lam_min,
                    ymax = lam_max ),
               alpha = 0.5,
               fill = 'yellow' ) +
  geom_line(   aes( x    = clim_x, 
                    y    = lam_mean ),
               lwd = 1.5 ) +
  ylab( expression(partialdiff*lambda[s]*'/'*partialdiff*'temperature') ) +
  ggtitle( 'Survival sensitivity' ) +
  ggsave( 'results/ipm/ltre/surv_sens.tiff',
          width=6.3,height=6.3,compression='lzw')
  


# set up list for storing lambdas
lam_s_l <- list()

# perform a loop for potential debugging
for(ii in 1:nrow(sim_df)){
  
  lam_s_l[ii] <- lam_stoch(sim_df$location[ii],
                           sim_df$clim_x[ii],
                           sim_df$vr[ii])
  
}

# put LTRE together 
ltre_df  <- sim_df %>% 
              mutate( lam_s    = unlist(lam_s_l) ) %>% 
              left_join(lam_base_df) %>% 
              mutate( lam_diff = lam_s - lam_s_base ) %>% 
              group_by( vr, clim_x ) %>% 
              summarise( lam_diff_mean = mean(lam_diff) ) %>% 
              ungroup %>% 
              rename( vital_rate = vr )


# Plot results 
ggplot(ltre_df) +
  geom_line( aes( x     = clim_x,
                  y     = lam_diff_mean,
                  color = vital_rate ),
             size = 1.5) +
  viridis::scale_color_viridis( discrete=T ) +
  ylab( expression(partialdiff*lambda[s]*'/'*partialdiff*'temperature') ) +
  xlab( 'Temperature anomaly' ) +
  theme( axis.title   = element_text( size = 20),
         legend.title = element_text( size = 15) ) +
  ggsave('results/ipm/ltre/ltre.tiff',
         width=6.3, height=4,compression='lzw')
