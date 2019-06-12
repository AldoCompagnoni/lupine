# PVA of lupine populations WITH DEMOGRAPHIC STOCHASTICITY
# 1. reate starting abundances for each population 
  # a. Calculate Stable stage distribution for each population (location)
  # b. Redistribute individuals across stages using observed numbers in 2018
# 2. Estimate and create stochastic IPM
# 3. IBM (not IPM) stochastic simulations and plot results
rm(list=ls())
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(rstan)
library(testthat)
library(lme4)


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


# 1. reate SSD for each population -------------------------------

# models 
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
                     surv_clim    = mean(all_vr$b_c_s),
                     
                     grow_b0      = mean(all_vr$a_u_g)*2,
                     grow_b1      = mean(all_vr$b_u_g)*2,
                     grow_sig     = mean(all_vr$sigma_y),

                     flow_b0      = mean(all_vr$a_u_f)*2,
                     flow_b1      = mean(all_vr$b_u_f)*2,
                     flow_clim    = mean(all_vr$b_c_f),
                     
                     fert_b0      = mean(all_vr$a_u_r)*2,
                     fert_b1      = mean(all_vr$b_u_r)*2,
                     fert_clim    = mean(all_vr$b_c_r),
                     
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


# update with LOCATION pars only: 
# for site-specific starting N vectors
updt_par_loc <- function(loc_n){
  
  pars_yr <- pars_mean
 
  # get location-specific parameters
  get_loc <- function(vr, loc_n){
  
    loc_ii  <- which(loc_n == loc_factors)
    var_n   <- paste0(vr, '.', loc_ii)
    all_vr[var_n][,1] %>% mean()
    
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
  pars_yr$surv_b0    <- mean(all_vr$a_u_s) + get_loc('a_loc_s',loc_n)
  pars_yr$surv_b1    <- mean(all_vr$b_u_s) + get_loc('b_loc_s',loc_n)
    
  pars_yr$grow_b0    <- mean(all_vr$a_u_g) + get_loc('a_loc_g',loc_n)
  pars_yr$grow_b1    <- mean(all_vr$b_u_g) + get_loc('b_loc_g',loc_n)
  
  pars_yr$flow_b0    <- mean(all_vr$a_u_f) + get_loc('a_loc_f',loc_n)
  pars_yr$flow_b1    <- mean(all_vr$b_u_f) + get_loc('b_loc_f',loc_n)

  pars_yr$fert_b0    <- mean(all_vr$a_u_r) + get_loc('a_loc_r',loc_n)
  pars_yr$fert_b1    <- mean(all_vr$b_u_r) + get_loc('b_loc_r',loc_n)
  
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
updt_par_loc("DR (3)")$surv_b0
updt_par_loc("ATT (8)")$clip
updt_par_loc("AL (1)")$clip



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
ker_avg_pop <- function(loc_n) kernel(0, updt_par_loc(loc_n) )
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



# 3. IBM stochastic simulations -------------------------------------

# FOR SIMULATIONS: update with year- and loc-specific parameters 
update_par_post <- function(year_n, loc_n, ii){
  
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
  pars_sim$surv_clim  <- all_vr['b_c_s'][thin_seq[ii],1]
  
  pars_sim$grow_b0    <- get_yr('a_yr_g',year_n,ii) + get_loc('a_loc_g',loc_n,ii)
  pars_sim$grow_b1    <- get_yr('b_yr_g',year_n,ii) + get_loc('b_loc_g',loc_n,ii)
  
  pars_sim$flow_b0    <- get_yr('a_yr_f',year_n,ii) + get_loc('a_loc_f',loc_n,ii)
  pars_sim$flow_b1    <- get_yr('b_yr_f',year_n,ii) + get_loc('b_loc_f',loc_n,ii)
  pars_sim$flow_clim  <- all_vr['b_c_f'][thin_seq[ii],1]
  
  pars_sim$fert_b0    <- get_yr('a_yr_r',year_n,ii) + get_loc('a_loc_r',loc_n,ii)
  pars_sim$fert_b1    <- get_yr('b_yr_r',year_n,ii) + get_loc('b_loc_r',loc_n,ii)
  pars_sim$fert_clim  <- all_vr['b_c_r'][thin_seq[ii],1]
  
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


# IBM functions
flx <- function(x,pars,tmp_anom){
  inv_logit( pars$flow_b0 + 
             pars$flow_b1*x +
             pars$flow_clim*tmp_anom )  
}

frx <- function(x,pars,tmp_anom){
  exp( pars$fert_b0 + 
       pars$fert_b1*x +
       pars$fert_clim*tmp_anom )  
}
 
# initial population structure

# aboveground plants
n   <- pars_mean$mat_siz
L   <- pars_mean$L 
U   <- pars_mean$U
h   <- (U-L)/n                   #Bin size
b   <- L+c(0:n)*h                #Lower boundaries of bins 
y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints

# n of indiv at stable stag. distrib. (based )
ssd_n     <- round(start_v_l[[5]])

# vector of above-ground sizes
ag_z      <- purrr::map2(ssd_n[-c(1:2)], y, 
                         function(x,y) rep(y, x) ) %>% 
                unlist
tmp_anom  <- 0
n_vec0    <- c( ssd_n[1:2], ag_z )


# project an IBM. Notation from Ellner et al. 2016 book
proj_ibm  <- function(n_vec0, pars, tmp_anom){

  # vector for seed bank 
  sb        <- n_vec0[1:2]
  # vector for above-ground individuals
  z         <- n_vec0[-c(1:2)]
  
  # project survivors
  surv_n1   <- rbinom( n    = length(z), 
                       size = 1, 
                       prob = sx(z, pars, tmp_anom) ) %>% 
                  as.logical
  
  # flowering individuals
  fl0       <- rbinom( n    = length(z), 
                       size = 1, 
                       prob = flx(z, pars, tmp_anom) )
  
  # avg number of racemes per indiv.
  rac0      <- rpois( n     = length(z),
                      lambda= frx(z, pars, tmp_anom) 
                      ) * fl0

  # viable racemes
  via_rac0  <- rbinom( n    = length(rac0),
                       size = rac0,
                       prob = (1 - pars$abort) )
  # via_rac0  <- round(rac0 * (1 - pars$abort), 0)
    
  # intact racemes
  int_rac0  <- rbinom( n    = length(via_rac0),
                       size = via_rac0,
                       prob = (1 - pars$clip) )
  # int_rac0  <- round(via_rac0 * (1 - pars$clip), 0)
  
  seed_n0   <- round(int_rac0 * pars$fruit_rac * pars$seed_fruit, 
                     0) %>% sum
  
  # transition to recruits, 1-year sb (g1), and 2-year sb (g2)
  r0      <- rbinom(1, size = seed_n0, prob = pars$g0 )
  r1      <- rbinom(1, size = seed_n0, prob = pars$g1 )
  r2      <- rbinom(1, size = seed_n0, prob = pars$g2 )
  # r0      <- round(seed_n0 * pars$g0 )
  # r1      <- round(seed_n0 * pars$g1 )
  # r2      <- round(seed_n0 * pars$g2 )
  
  # update seeds
  sb1     <- c(r2, sb[1] + r1)
  
  # new recruits
  z1_r    <- rnorm( n    = r0 + sb[2],
                    mean = pars$recr_sz,
                    sd   = pars$recr_sd )
  
  # growth of surviving individuals
  z1_g    <- rnorm( n    = length(z[surv_n1]),
                    mean = pars$grow_b0 + pars$grow_b1*z[surv_n1], 
                    sd   = pars$grow_sig )
  
  z1      <- c(z1_r, z1_g)
  
  n_vec1  <- c(sb1, z1) 
  
  n_vec1
  
}



# sequence of years (always the same), and total IPM size
it        <- 50
yrs       <- 2005:2017
loc_v     <- lupine_df$location %>% unique %>% sort
ker_siz   <- (pars_mean$mat_siz) + 2
clim_x    <- seq( 0,
                  max(clim_mat$tmp_t0), 
                  length.out = 10 )

# set up 50-year simulations
sim_df    <- expand.grid( post_ii = 1:100,
                          sim_i   = 1:50,
                          clim_x  = 0,
                          # clim_x = clim_x,
                          loc_v   = 'BR (6)',
                          # loc_v  = loc_v[-c(1,2,7)],
                          stringsAsFactors = F )

# calculate stochastic lambda for each climate anomaly
project_50 <- function(ii){

  # create starting population vector state (n_vec0) 
  
  # stabe stage distribution (integers)
  ssd_n    <- start_v_l[sim_df$loc_v[ii]] %>% 
                unlist %>% 
                as.numeric %>% 
                round()
  
  # aboveground plants
  n   <- pars_mean$mat_siz
  L   <- pars_mean$L 
  U   <- pars_mean$U
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  
  # vector of above-ground sizes
  ag_z      <- purrr::map2(ssd_n[-c(1:2)], y, 
                           function(x,y) rep(y, x) ) %>% unlist
  tmp_anom  <- sim_df$clim_x[ii]
  n_vec0    <- c( ssd_n[1:2], ag_z )

  # random year sequence 
  it      <- 50
  yrs     <- 2005:2017
  set.seed(sim_df$sim_i[ii])
  yr_seq  <- sample(yrs,it,replace=T)
  
  n_it    <- list( n_vec0 )
  n_tot   <- rep(NA,it)
  n_tot[1]<- length( n_vec0[-c(1:2)] )
  
  # stochastic simulations
  for(pp in 2:it){
    
    pars_in    <- update_par_post( yr_seq[pp],
                                   sim_df$loc_v[ii],
                                   sim_df$post_ii[ii])
    n_it[[pp]] <- proj_ibm( n_it[[pp-1]], 
                            pars_in, 
                            sim_df$clim_x[ii] )
    n_tot[pp]  <- length(n_it[[pp]][-c(1:2)])
  
  }
   
  n_tot
  
}

# all population sizes
init_t <- Sys.time()
lapply(1:10,project_50)
Sys.time() - init_t

proj_l <- lapply(1:nrow(sim_df),project_50)

# plot a set of stochastic simulations
proj_df <- Reduce(function(...) rbind(...), proj_l ) %>%
              as.data.frame %>%
              setNames( paste0('V',1:50) ) %>%
              bind_cols( sim_df, .) %>% 
              gather( time, n, V1:V50 ) %>%
              mutate( time = gsub('V','',time) ) %>%
              mutate( time = as.numeric(time) ) 



ggplot(proj_df) +
    geom_line( aes(x = time,
                 y = n,
                 group = sim_i),
             alpha=0.3) +
  theme( legend.position = 'none',
         axis.text.x = element_text(angle=90) ) +
  labs( y = 'Population size',
        x = 'Time step') 
#   ggsave( 'results/ipm/qe_sims/DR_clim0_stoch_sim.tiff',
#           width=6.3,height=6.3,compressio='lzw')





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

# # store, just in case
# save.image(paste0('results/ipm/qe_sims/eq_',
#                   thresh,
#                   '.RData'))