# IPM from data
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
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"

all_indiv_sample <- c("BS (7)", 'DR (3)', 'NB (2)')

# data
lupine_df   <- read.csv( "data/lupine_all.csv") #%>% 
                  # subset( location %in% c("BS (7)") ) %>%
                  # subset( year > 2008 )
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')
cons        <- read_xlsx('data/consumption.xlsx') %>% 
                mutate( Mean_consumption = Mean_consumption %>% as.numeric) %>% 
                select( Year, Site, Mean_consumption) %>% 
                # expand potential "cases"
                complete( Site, Year) %>% 
                # update name
                mutate( Site = toupper(Site) )
pred_g      <- read_xlsx('data/post predation_lupinus tidestromii.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")

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
  
subset(clim, clim_var == 'tmean') %>%  
  group_by(year) %>% 
  summarise(mean_t = mean(clim_value) ) %>% 
  as.data.frame

enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_anom('oni')

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat,tmp_mat,enso_mat) )


# vital rates format --------------------------------------------------------------
surv_all    <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 


surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2,
                          log_area_t03= log_area_t0^3 ) %>% 
                  left_join( clim_mat ) 

# write suspect seedlings for Eleanor
seedl %>% 
  subset( log_area_t0 > 4 ) %>% 
  mutate( year = year - 1 ) %>% 
  select(location,year,newid) %>% 
  unique %>% 
  write.csv('results/explorative_graphs/suspect_seedlings.csv',
            row.names=F)


sl_grow     <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2 ) %>% 
                  # only seedling GROWTH
                  subset( surv_t1 == 1 ) %>% 
                  left_join( clim_mat ) 

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset( !(stage_t0 %in%  "SL") ) %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM", "SL")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1) %>% 
                  left_join( clim_mat ) 

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  left_join( clim_mat ) 

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) 

abor        <- subset(lupine_df, !is.na(flow_t0) & flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove 
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t02 = log(area_t0)^2) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  # only years indicated by Tiffany
                  subset( year %in% c(2010, 2011, 2013:2018) ) %>% 


# models ---------------------------------------------------------
mod_sl   <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=seedl, family='binomial')
mod_g_sl <- lm( log_area_t1 ~ log_area_t0, data=sl_grow)
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
mod_ab    <- glmer( cbind(numab_t0, numrac_t0-numab_t0) ~ 
                    (1 | year) + (1 | location),
                    data=abor, family='binomial')
# mod_fr   <- MASS::glm.nb(numrac_t0 ~ log_area_t0 , data=fert )
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))
germ_coef<- select(germ, g0:g2) %>% colMeans


# vital rate models 
surv_sl_p <- fixef(mod_sl)
grow_sl_p <- coef(mod_g_sl)
grow_sl_p <- c(grow_sl_p, summary(mod_g_sl)$sigma)
surv_p    <- fixef(mod_s)
grow_p    <- fixef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- fixef(mod_fl)
fert_p    <- fixef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ_coef * (1 - 0.43)
# abortion
# abor_p    <- list( # set up "grid" of year/location combinations
#                    expand.grid( year = paste0(c(2005:2018)),
#                                 site = rownames(coef(mod_ab)$location),
#                                 stringsAsFactors = F ),
#                    # set up year specific coefficients
#                    data.frame(  year   = rownames(coef(mod_ab)$year),
#                                 year_c = coef(mod_ab)$year[,1] ),
#                    # set up site specific coefficients
#                    data.frame(  site   = rownames(coef(mod_ab)$location),
#                                 site_c = coef(mod_ab)$location[,1] ) ) %>% 
#               Reduce( function(...) full_join(...), . ) %>% 
#               # calculate proportion of aborted racemes
#               mutate( tot_p  = boot::inv.logit(year_c) + boot::inv.logit(site_c) ) %>% 
#               rename( location = site )
abor_p    <-  abor %>% 
                mutate( ab_p = numab_t0 / numrac_t0 ) %>%
                group_by( year, location ) %>% 
                summarise( ab_p_m = mean(ab_p,na.rm=T) ) %>% 
                ungroup  

ggplot(abor_p) +
geom_line( aes(x     = year, 
               y     = ab_p_m, 
               color = location),
           size=1.5) +
viridis::scale_color_viridis(discrete=T) +
# geom_point( data = subset(abor_p, !is.na(tot_p) ),
#             aes(x=year, 
#                 y=tot_p,
#                 shape=location) ) 
  

# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }

# list of mean IPM parameters. 
pars_mean   <- list( # seedlings vital rates
                     surv_sl_b0   = surv_sl_p['(Intercept)'],
                     surv_sl_b1   = surv_sl_p['log_area_t0'],
                     surv_sl_b2   = surv_sl_p['log_area_t02'],
                     surv_sl_b3   = surv_sl_p['log_area_t03'],
                     surv_sl_clim = surv_sl_p['tmp_tm1'],
                     
                     grow_sl_b0   = grow_sl_p["(Intercept)"],
                     grow_sl_b1   = grow_sl_p["log_area_t0"],
                     grow_sl_sig  = grow_sl_p[3], 
                     
                     # adults vital rates           
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
                     
                     one_sz       = sl_grow$log_area_t1 %>% mean(na.rm=T),
                     one_sd       = sl_grow$log_area_t1 %>% sd(na.rm=T),
                     
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
update_par_spec <- function(year_n, vr, loc_n){
  
  pars_yr <- pars_mean
  
  get_yr <- function(mod_obj,year_n){
  
    # # models used in year
    # yr_mod <- attributes(mod_obj) %>%
    #             .$flist %>%
    #             .$year %>%
    #             as.character %>%
    #             unique
    # 
    # if( all(!paste0(2006:2008) %in% yr_mod) ){
    #   year_n <- replace( year_n, year_n== 2006, 2009)
    #   year_n <- replace( year_n, year_n== 2007, 2010)
    #   year_n <- replace( year_n, year_n== 2008, 2011)
    # }
    
    ranef(mod_obj)$year %>% 
      tibble::add_column(.,
                         .before=1,
                         year=row.names(.)) %>% 
      subset( year == year_n) %>% 
      select(-year)
  }
  
  get_loc <- function(mod_obj,loc_n){
    ranef(mod_obj)$location %>% 
      tibble::add_column(.,
                         .before=1,
                         location=row.names(.)) %>% 
      subset( location == loc_n) %>% 
      select(-location)
  }
  
  # get info on raceme "loss"
  get_rac <- function(cons, year_n, loc_n){
    
    # raceme consumption
    locc     <- gsub(' \\([0-9]\\)','',loc_n)
    out_cons <- cons %>% 
                  subset(Site == locc & Year == year_n ) %>% 
                  .$Mean_consumption %>% 
                  as.numeric
    if( is.na(out_cons) ) out_cons <- mean(cons$Mean_consumption,
                                           na.rm=T)
    
    # raceme abortion
    out_abor <- abor_p %>% 
                  subset( site == loc_n & year == year_n ) %>% 
                  .$abor_p
    
    if( is.na(out_abor) ) out_abor <- mean(abor_p$abor_p,na.rm=T)
    
    data.frame( clip = out_cons,
                abor = out_abor)
    
  }
  
  if( vr == 'surv_sl'){
    pars_yr$surv_sl_b0 <- (pars_yr$surv_sl_b0 + get_yr(mod_sl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_sl_b1 <- (pars_yr$surv_sl_b1 + get_yr(mod_sl,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$surv_sl_b0 <- (pars_yr$surv_sl_b0 + get_loc(mod_sl,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_sl_b1 <- (pars_yr$surv_sl_b1 + get_loc(mod_sl,loc_n)['log_area_t0']) %>% as.numeric
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
    pars_yr$surv_sl_b0 <- (pars_yr$surv_sl_b0 + get_yr(mod_sl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_sl_b1 <- (pars_yr$surv_sl_b1 + get_yr(mod_sl,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_yr(mod_s,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_yr(mod_s,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_yr(mod_g,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_yr(mod_g,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_yr(mod_fl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_yr(mod_fl,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_yr(mod_fr,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_yr(mod_fr,year_n)['log_area_t0']) %>% as.numeric
    
    # locations
    pars_yr$surv_sl_b0 <- (pars_yr$surv_sl_b0 + get_loc(mod_sl,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_sl_b1 <- (pars_yr$surv_sl_b1 + get_loc(mod_sl,loc_n)['log_area_t0']) %>% as.numeric
  
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
  pars_yr$clip       <- get_rac(cons, year_n, loc_n)$clip
  pars_yr$abort      <- get_rac(cons, year_n, loc_n)$abor
  
  pars_yr
  
}

# test function
update_par_spec(2006,'all',"BS (7)")$surv_sl_b0
update_par_spec(2010,'all',"POP9 (9)")$abor
update_par_spec(2010,'all',"POP9 (9)")$g0

# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# Transforms all values below/above limits in min/max size
x_range_s <- function(x,pars){
  pmin(pmax(x,pars$L_sl),pars$U_sl)
}

# Survival at size x
sx_s<-function(x,pars,tmp_anom){
  xb <- x_range_s(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_sl_b0 + 
                    pars$surv_sl_b1 * xb + 
                    pars$surv_sl_b2 * (xb^2) + 
                    pars$surv_sl_b3 * (xb^3) + 
                    pars$surv_sl_clim * tmp_anom) )
}

# update kernel functions
grow_sl_sd <- function(x,pars){
  pars$a_sl*(exp(pars$b_sl*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy_s <- function(x,y,pars){
  xb <- x_range_s(x, pars)
  # returns a *probability density distribution* for each x value
  # return( dnorm(y,  mean = pars$grow_sl_b0 + pars$grow_sl_b1*xb, 
  #                   sd   = pars$grow_sl_sig) )
  return( dnorm(y,  mean = pars$one_sz, 
                    sd   = pars$one_sd) )
}

# transition: Survival * growth
pxy_s<-function(x,y,pars,tmp_anom){
  xb <- x_range_s(x, pars)
  return( sx_s(xb,pars,tmp_anom) * gxy_s(xb,y,pars) )
}

x_range <- function(x,pars){
  pmin(pmax(x,pars$L),pars$U)
}

# Survival at size x
sx<-function(x,pars,tmp_anom){
  xb <- x_range(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + 
                    pars$surv_b1 * xb +
                    pars$surv_b2 * xb^2 + 
                    pars$surv_b3 * xb^3 + 
                    pars$surv_clim * tmp_anom) )
}

# update kernel functions
grow_sd <- function(x,pars){
  pars$a*(exp(pars$b*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy <- function(x,y,pars){
  xb <- x_range(x, pars)
  # returns a *probability density distribution* for each x value
  return( dnorm(y,  mean = pars$grow_b0 + pars$grow_b1*xb, 
                    sd   = pars$grow_sig) )
}

# transition: Survival * growth
pxy<-function(x,y,pars,tmp_anom){
  xb <- x_range(x, pars)
  return( sx(xb,pars,tmp_anom) * gxy(xb,y,pars) )
}


# production of seeds from x-sized mothers
fx <-function(x,pars,tmp_anom){
  
  xb       <- x_range(x, pars)
  
  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + 
                         pars$flow_b1   * xb + 
                         pars$flow_clim * tmp_anom ) * 
              exp(       pars$fert_b0 + 
                         pars$fert_b1   * xb + 
                         pars$fert_clim * tmp_anom )
              
  # viable racs
  viab_rac <- tot_rac * ( 1 - (pars$abort+pars$clip) )
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars,h){
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd ) * h
}

# seeds to seedling direct transition
fxy_s <- function(x,y,pars,tmp_anom,h){
  xb <- x_range_s(x, pars)
  fx(xb,pars,tmp_anom) * pars$g0 * recs(y, pars, h)
}


# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(tmp_anom,pars){
 
  # set up IPM domains --------------------------------------------------------
  
  # seedlings
  n   <- pars$mat_siz
  L_s <- pars$L_sl 
  U_s <- pars$U_sl
  #these are the upper and lower integration limits
  h_s <- (U_s-L_s)/n                   #Bin size
  b_s <- L_s+c(0:n)*h_s                #Lower boundaries of bins 
  y_s <- 0.5*(b_s[1:n]+b_s[2:(n+1)])   #Bins' midpoints

  # plants
  n   <- pars$mat_siz
  L   <- pars$L 
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # Kernel "pieces" --------------------------------------------------------
  
  # vectors
  plant_s1  <- numeric(n)
  plant_s2  <- numeric(n)
  rec_s1    <- numeric(n)
  rec_s2    <- numeric(n)
  s1_rec    <- numeric(n)
  s2_rec    <- numeric(n)
  s1_plant  <- numeric(n)
  s2_plant  <- numeric(n)
  
  # seeds 
  s_mat     <- matrix(0,2,2)  
   
  # populate kernel ------------------------------------------------------------
  
  # seeds that enter 1 yr-old seed bank 
  plant_s1  <- fx(y,pars,tmp_anom) * (1 - pars$g0)
  # no seeds go directly to 2 yr-old seed bank!
  plant_s2  <- numeric(n)

  # Impossible: seedlings to seedling
  sl_sl1     <- matrix(0,n,n)
  
  # seeds that go directly to seedings germinate right away 
  pl_sl2     <- (outer(y,y_s, fxy_s, pars, tmp_anom, h_s) ) %>% t 
  # pl_sl2     <- matrix(0,n,n)
  
  # recruits from the 1 yr-old seedbank 
  s1_rec[]   <- recs(y_s, pars, h_s) * pars$g1
  # s1_rec[]   <- rep(0,100)
  
  # seeds that enter 2 yr-old seed bank 
  s_mat[2,1] <- (1 - pars$g1)
  
  # recruits from the 2 yr-old seedbank 
  s2_rec[]   <- recs(y_s, pars, h_s) * pars$g2
  # s2_rec[]   <- rep(0,100)
  
  # seedling survival and growth into plants
  sl_pl3     <- (outer(y_s, y, pxy_s, pars, tmp_anom)*h) %>% t
  # sl_pl3     <- matrix(0,n,n) 
  
  # survival and growth of adult plants
  pl_pl4     <- (outer(y,y,pxy,pars,tmp_anom)*h) %>% t
  
  # utilities for check
  image_k <- function(x) t(apply(x, 2, rev)) %>% image
  # sl_pl3 %>% image_k
  # pl_pl4 %>% image_k
  
  # Assemble the kernel -------------------------------------------------------------
  
  # "big blocks"
  big_mat    <- rbind( cbind(sl_sl1,pl_sl2),
                       cbind(sl_pl3,pl_pl4) )
  
  # top 2 vectors 
  top_vec    <- rbind( c(rec_s1, plant_s1),
                       c(rec_s2, plant_s2) )

  side_vec   <- rbind( s_mat,
                       cbind(s1_rec,   s2_rec  ),
                       cbind(s1_plant, s2_plant) )
  
  
  k_yx       <- cbind( side_vec, 
                       rbind(top_vec, big_mat) 
                       )
     
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
yrs       <- 2009:2018
loc_v     <- lupine_df$location %>% unique %>% sort
set.seed(1776)
yr_seq    <- sample(yrs,it,replace=T)
yr_i      <- yr_seq - 2008
ker_siz   <- (pars_mean$mat_siz*2)+2


# calculate stochastic lambda for each climate anomaly
lam_stoch <- function(tmp_anom, vr, loc){

  # store kernels
  ker_l     <- list()
  for(ii in seq_along(yrs)){
    ker_l[[ii]] <- kernel(tmp_anom,update_par_spec(yrs[ii],vr,loc))
  }

  # placeholders
  n0        <- rep(100/ker_siz,ker_siz)
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
lamS_noClim <- lam_stoch(0,'all', loc_v[6])
lamS_surv_sl<- lam_stoch(0,'surv_sl')
lamS_surv   <- lam_stoch(0,'surv')
lamS_grow   <- lam_stoch(0,'grow')
lamS_flow   <- lam_stoch(0,'flow')
lamS_fert   <- lam_stoch(0,'fert')


# set up climate covariates
clim_x <- seq( min(clim_mat$tmp_tm1),
               max(clim_mat$tmp_tm1), 
               length.out = 10 )

lam_s_l <- list()

# cycle through the different locations
for(li in 1:6){ #length(loc_v)
  lam_s_vec      <- sapply(clim_x, lam_stoch, 'all', loc_v[li])
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
  ggsave('results/ipm/lambdaS_vs_site_clim.tiff',
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
