# 1. Load in all data and Bayesian model results
# 2. Set up data frames for year-by-site plots 
# 3. Compute predictions 
# 4. Plot site/year data and model results
# 5. Plot "average" plots with posteriors
# 6. "ribbon" plots with high/low temperature
rm(list=ls())
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(rstan)
library(ggplot2)
library(testthat)


# data
lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                  mutate( log_area_t0  = log(area_t0) ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3 )
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")
fit_mod     <- readRDS('C:/Users/ac22qawo/lupine/lupine_vr_bayes-5270385_lupine_surv_nc.RDS')
fit_mod     <- readRDS('C:/Users/ac22qawo/lupine/lupine_vr_bayes-5270262_lupine_nc.RDS')

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

surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          surv_t1, tmp_t0, tmp_tm1 )

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          log_area_t1, tmp_t0 )

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          flow_t0, tmp_t0 ) %>% 
                  subset( year != 2018 )

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) %>% 
                  select( location, year, newid,
                          log_area_t0, log_area_t02, log_area_t03,
                          numrac_t0, tmp_t0 ) %>% 
                  subset( year != 2018 )


# All vital rates at once
vr_all <- Reduce( function(...) full_join(...), list( surv, grow,
                                                      flow, fert) )

# 2. Set up data frames for year-by-site plots -----------------------------------------------


# data frame of binned proportions
df_binned_prop <- function(ii, df_in, n_bins, siz_var, rsp_var, grid_y_l){
  
  # make sub-selection of data
  df   <- subset(df_in, year     == grid_y_l$year[ii] & 
                        location == grid_y_l$location[ii] )
  
  if( nrow(df) == 0 ) return( NULL)
  
  size_var <- deparse( substitute(siz_var) )
  resp_var <- deparse( substitute(rsp_var) )
  
  # binned survival probabilities
  h    <- (max(df[,size_var],na.rm=T) - min(df[,size_var],na.rm=T)) / n_bins
  lwr  <- min(df[,size_var],na.rm=T) + (h*c(0:(n_bins-1)))
  upr  <- lwr + h
  mid  <- lwr + (1/2*h)
  
  binned_prop <- function(lwr_x, upr_x, response){
    
    id  <- which(df[,size_var] > lwr_x & df[,size_var] < upr_x) 
    tmp <- df[id,]
    
    if( response == 'prob' ){   return( sum(tmp[,resp_var],na.rm=T) / nrow(tmp) ) }
    if( response == 'n_size' ){ return( nrow(tmp) ) }
    
  }
  
  y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
  x_binned <- mid
  y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist
  
  # output data frame
  data.frame( xx = x_binned, 
              yy  = y_binned,
              nn  = y_n_size) %>% 
    setNames( c(size_var, resp_var, 'n_size') ) %>% 
    mutate( year     = grid_y_l$year[ii], 
            location = grid_y_l$location[ii] )
  
}


# grid of year/location for survival and flowering
grid_y_l    <- expand.grid( year     = surv$year %>% unique %>% sort,
                            location = surv$location %>% unique %>% sort,
                            stringsAsFactors = F)

# produce lists to combine in data frames 
surv_bin_l  <- lapply(1:nrow(grid_y_l), df_binned_prop, surv, 10, 
                                        log_area_t0, surv_t1, grid_y_l)

flow_bin_l  <- lapply(1:nrow(grid_y_l), df_binned_prop, flow, 10, 
                                        log_area_t0, flow_t0, grid_y_l)

# big data frames for "panel plots"
surv_pan_df <- bind_rows( surv_bin_l ) %>% 
                  mutate( loc_lab = location ) %>% 
                  mutate( log_area_t02 = log_area_t0^2, 
                          log_area_t03 = log_area_t0^3,
                          tmp_tm1      = 0,
                          transition   = paste( paste0(year - 1), 
                                                substr(paste0(year),3,4),
                                                sep='-') ) %>% 
                  mutate( year         = as.integer(year - 2004),
                          location     = location %>% as.factor %>% as.integer )

flow_pan_df <- bind_rows( flow_bin_l ) %>% 
                  mutate( loc_lab = location ) %>% 
                  mutate( transition   = paste( paste0(year - 1), 
                                                substr(paste0(year),3,4),
                                                sep='-') ) %>% 
                  mutate( year         = as.integer(year - 2004),
                          location     = location %>% as.factor %>% as.integer )

# growth
grow_pan_df <- grow %>% 
                  mutate( loc_lab = location ) %>% 
                  mutate( transition   = paste( paste0(year - 1), 
                                                substr(paste0(year),3,4),
                                                sep='-') ) %>% 
                  mutate( year         = as.integer(year - 2004),
                          location     = location %>% as.factor %>% as.integer )

# fertility. ONLY IN THIS CASE, two separate data frames

# data frame with data
fert_pan_df <- fert %>% 
                  mutate( loc_lab    = location ) %>% 
                  mutate( transition = paste( paste0(year - 1), 
                                                substr(paste0(year),3,4),
                                                sep='-') ) %>% 
                  mutate( year      = as.integer(year - 2004),
                          location  = location %>% as.factor %>% as.integer )

# data frame of predictions
x_fert      <- seq( min(fert$log_area_t0),
                    max(fert$log_area_t0), 
                    length.out=20 )

# retain only observed year/site combinations
fert_all    <- fert %>% 
                  mutate(   year         = as.integer(year - 2004),
                            location     = location %>% as.factor %>% as.integer ) %>% 
                  select(year, location) %>% 
                  unique
fert_pred_df <- expand.grid( year        = surv$year %>% unique %>% sort,
                             location    = surv$location %>% unique %>% sort,
                             log_area_t0 = x_fert,
                             stringsAsFactors = F ) %>% 
                  mutate( loc_lab    = location ) %>% 
                  mutate( transition   = paste( paste0(year - 1), 
                                                substr(paste0(year),3,4),
                                                sep='-') ) %>% 
                  mutate(   year         = as.integer(year - 2004),
                            location     = location %>% as.factor %>% as.integer ) %>% 
                  right_join( fert_all )



# 3. Compute predictions ----------------------------------------------

# make grid for predictions
pred_grid <- expand.grid( year     = 1:13,
                          location = 1:7,
                          stringsAsFactors = F ) 

# get all parameters for a vital rate
vr_pars_get  <- function(mod_obj, vr){
  
  # extract data frame of mean parameters
  mean_pars <- function(x){
    
    as.matrix(x) %>% 
        apply(2,mean) %>% 
        matrix( nrow=1, ncol=length(.) ) %>% 
        as.data.frame %>% 
        setNames( colnames( as.matrix(x)) ) %>% 
        setNames( gsub('\\[','_',names(.)) ) %>% 
        setNames( gsub('\\]','', names(.)) )  
  
  }

  # extract number from string
  get_num <- function(x){
    regmatches(x, 
           gregexpr("[[:digit:]]{1,2}", 
           x) ) %>% 
      unlist %>% 
      as.integer
  }
  
  # get all parameters from model
  all_pars <- mean_pars(mod_obj)
  
  # vectors of params
  a_yr_v    <- paste0('a_yr_', vr, '_', 1:13)
  b_yr_v    <- paste0('b_yr_', vr, '_', 1:13)
  a_loc_v   <- paste0('a_loc_',vr, '_', 1:7)
  b_loc_v   <- paste0('b_loc_',vr, '_', 1:7)
  
  # random year intercept
  a_yr_df   <- all_pars %>% 
      select( a_yr_v ) %>% 
      stack %>% 
      mutate( ind = as.character(ind) ) %>% 
      mutate( year = get_num( ind ) ) %>% 
      rename( a_yr = values ) %>% 
      select( -ind )

  # random year slope
  b_yr_df   <- all_pars %>% 
      select( b_yr_v ) %>% 
      stack %>% 
      mutate( ind = as.character(ind) ) %>% 
      mutate( year = get_num( ind ) ) %>%  
      rename( b_yr = values ) %>% 
      select( -ind )

  # random location intercept
  a_loc_df  <- all_pars %>% 
      select( a_loc_v ) %>% 
      stack %>% 
      mutate( ind      = as.character(ind) ) %>% 
      mutate( location = get_num( ind ) ) %>%   
      rename( a_loc    = values ) %>% 
      select( -ind )
    
  # random location slope
  b_loc_df  <- all_pars %>% 
      select( b_loc_v ) %>% 
      stack %>% 
      mutate( ind      = as.character(ind) ) %>% 
      mutate( location = get_num( ind ) ) %>%   
      rename( b_loc    = values ) %>% 
      select( -ind )

  # ugly but works
  pred_grid %>% 
    full_join( a_yr_df ) %>% 
    full_join( b_yr_df ) %>% 
    full_join( a_loc_df ) %>% 
    full_join( b_loc_df ) %>% 
    mutate( b_s2 = all_pars$b_s2,
            b_s3 = all_pars$b_s3 )

}

# Extract parameters
s_pars_df <- vr_pars_get(fit_mod, 's')
g_pars_df <- vr_pars_get(fit_mod, 'g')
f_pars_df <- vr_pars_get(fit_mod, 'f')
r_pars_df <- vr_pars_get(fit_mod, 'r')


# Create predictions
surv_df <- left_join(surv_pan_df, s_pars_df) %>% 
              mutate( y_raw = a_yr + a_loc + 
                              ((b_yr + b_loc)* log_area_t0) + 
                              b_s2 * log_area_t02 +
                              b_s3 * log_area_t03 ) %>% 
              mutate( yhat  = boot::inv.logit(y_raw) )

grow_df <- left_join(grow_pan_df, g_pars_df) %>% 
              mutate( alpha = a_yr + a_loc,
                      beta  = b_yr + b_loc )
              

flow_df <- left_join(flow_pan_df, f_pars_df) %>% 
              mutate( y_raw = a_yr + a_loc + 
                              ((b_yr + b_loc)* log_area_t0) ) %>% 
              mutate( yhat  = boot::inv.logit(y_raw) )

fert_df <- left_join(fert_pred_df, r_pars_df) %>% 
              mutate( y_raw = a_yr + a_loc + 
                              ((b_yr + b_loc)* log_area_t0) ) %>% 
              mutate( yhat  = exp(y_raw) )
 

# 4. Plot site/year data and model results ---------------------------


# survival
ggplot(data  = surv_df, 
       aes(x = log_area_t0, 
           y = surv_t1) ) +
  geom_point(alpha = 1,
             pch   = 16,
             size  = 0.5,
             color = 'red') +
  geom_line(aes(x = log_area_t0,
                y = yhat),
            lwd = 1,
            alpha = 0.5)+
  # split in panels
  facet_grid(loc_lab ~ transition) +
  theme_bw() +
  theme( axis.text = element_text( size = 5 ),
         title     = element_text( size = 10 ),
         strip.text.y  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.text.x  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         panel.spacing = unit('0.5',unit='mm') ) +
  ggtitle("Survival of all plants" ) +
  ggsave(filename = "results/vital_rates/surv_all_bayes.tiff",
         dpi = 300, width = 6.3, height = 4, units = "in",
         compression = 'lzw')
                 

# growth
ggplot(data  = grow_df, 
       aes(x = log_area_t0, 
           y = log_area_t1) ) +
  geom_point(alpha = 1,
             pch   = 16,
             size  = 0.5,
             color = 'red') +
  geom_abline( aes(intercept = alpha,
                   slope     = beta),
            lwd = 1,
            alpha = 0.5) +
  # split in panels
  facet_grid(loc_lab ~ transition) +
  theme_bw() +
  theme( axis.text = element_text( size = 5 ),
         title     = element_text( size = 10 ),
         strip.text.y  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.text.x  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         panel.spacing = unit('0.5',unit='mm') ) +
  ggtitle("Growth" ) +
  ggsave(filename = "results/vital_rates/growth_bayes.tiff",
         dpi = 300, width = 6.3, height = 4, units = "in",
         compression = 'lzw')
      

# prob of flowering
ggplot(data  = flow_df, 
       aes(x = log_area_t0, 
           y = flow_t0) ) +
  geom_point(alpha = 1,
             pch   = 16,
             size  = 0.5,
             color = 'red') +
  geom_line(aes(x = log_area_t0,
                y = yhat),
            lwd = 1,
            alpha = 0.5)+
  # split in panels
  facet_grid(loc_lab ~ transition) +
  theme_bw() +
  theme( axis.text = element_text( size = 5 ),
         title     = element_text( size = 10 ),
         strip.text.y  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.text.x  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         panel.spacing = unit('0.5',unit='mm') ) +
  ggtitle("Probability of flowering" ) +
  ggsave(filename = "results/vital_rates/fert_bayes.tiff",
         dpi = 300, width = 6.3, height = 4, units = "in",
         compression = 'lzw')
        
 
# fertility
ggplot(data  = fert_pan_df, 
       aes(x = log_area_t0, 
           y = numrac_t0) ) +
  geom_point(alpha = 1,
             pch   = 16,
             size  = 0.5,
             color = 'red') +
  geom_line( data = fert_df,
              aes(x = log_area_t0,
                  y = yhat),
            lwd = 1,
            alpha = 0.5)+
  # split in panels
  facet_grid(loc_lab ~ transition) +
  theme_bw() +
  theme( axis.text = element_text( size = 5 ),
         title     = element_text( size = 10 ),
         strip.text.y  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.text.x  = element_text( size = 5,
                                       margin = margin(0.5,0.5,0.5,0.5,
                                                       'mm') ),
         strip.switch.pad.wrap = unit('0.5',unit='mm'),
         panel.spacing = unit('0.5',unit='mm') ) +
  ggtitle("Fertility" ) +
  ggsave(filename = "results/vital_rates/fert_bayes.tiff",
         dpi = 300, width = 6.3, height = 4, units = "in",
         compression = 'lzw')


# 5. Plot "average" plots with posteriors -------------------------------


tiff('results/vital_rates/bayes/vr_posterior.tiff',
     unit='in',res=300,height=6.3,width=6.3,compression='lzw')

par(mfrow=c(2,2),mar=c(3,3,0.1,0.1),mgp=c(1.6,0.5,0))

# extract data frame of mean parameters
mean_pars <- function(x){
  
  as.matrix(x) %>% 
      as.data.frame() %>% 
      setNames( colnames( as.matrix(x)) ) %>% 
      setNames( gsub('\\[','_',names(.)) ) %>% 
      setNames( gsub('\\]','', names(.)) )  

}

# survival

# all time
all_pars <- mean_pars(fit_mod)[round(seq(1,6000,length.out = 100),0),]
x_seq    <- seq( min(surv$log_area_t0),
                 max(surv$log_area_t0), length.out=20)

plot_binned_prop(surv, 10, log_area_t0, surv_t1)

# draw 100 posterior samples
draw_posterior <- function(ii){
  y_raw    <- (all_pars$a_u_s[ii]*2) + 
              (all_pars$b_u_s[ii]*2) * x_seq +
               all_pars$b_s2[ii] * (x_seq^2) +
               all_pars$b_s3[ii] * (x_seq^3) +
               all_pars$b_c_s[ii] * mean(surv$tmp_t0)
  lines(x_seq, boot::inv.logit(y_raw), col='grey' )
}

sapply(1:100,draw_posterior)
par(new=TRUE)
plot_binned_prop(surv, 10, log_area_t0, surv_t1)


# Growth
draw_posterior <- function(ii){
  abline( a = (all_pars$a_u_g[ii]*2), 
          b = (all_pars$b_u_g[ii]*2),
          col = 'grey' )
}

plot(log_area_t1 ~ log_area_t0, data=grow)
sapply(1:100,draw_posterior)


# Flowering
draw_posterior <- function(ii){
  y_raw    <- (all_pars$a_u_f[ii]*2) + 
              (all_pars$b_u_f[ii]*2) * x_seq +
               all_pars$b_c_f[ii] * mean(flow$tmp_t0)
  lines(x_seq, boot::inv.logit(y_raw), col='grey' )
}

x_seq    <- seq( min(flow$log_area_t0),
                 max(flow$log_area_t0), length.out=20)

plot_binned_prop(flow, 10, log_area_t0, flow_t0)
sapply(1:100,draw_posterior)
par(new=T)
plot_binned_prop(flow, 10, log_area_t0, flow_t0)



# Fertility
draw_posterior <- function(ii){
  y_raw    <- (all_pars$a_u_r[ii]*2) + 
              (all_pars$b_u_r[ii]*2) * x_seq +
               all_pars$b_c_r[ii] * mean(fert$tmp_t0)
  lines(x_seq, exp(y_raw), col='grey' )
}

x_seq    <- seq( min(fert$log_area_t0),
                 max(fert$log_area_t0), length.out=20)

plot(numrac_t0 ~ log_area_t0, data=fert)
sapply(1:100,draw_posterior)
par(new=T)
plot(numrac_t0 ~ log_area_t0, data=fert)

dev.off()






# extract data frame of mean parameters
mean_pars <- function(x){
  
  as.matrix(x) %>% 
      as.data.frame() %>% 
      setNames( colnames( as.matrix(x)) ) %>% 
      setNames( gsub('\\[','_',names(.)) ) %>% 
      setNames( gsub('\\]','', names(.)) )  

}


# 6. "ribbon" plots with high/low temperature ---------------

# produce dataframe with binned proportions
df_binned_prop <- function(df, n_bins, siz_var, rsp_var){
  
  size_var <- deparse( substitute(siz_var) )
  resp_var <- deparse( substitute(rsp_var) )
  
  # binned survival probabilities
  h    <- (max(df[,size_var],na.rm=T) - min(df[,size_var],na.rm=T)) / n_bins
  lwr  <- min(df[,size_var],na.rm=T) + (h*c(0:(n_bins-1)))
  upr  <- lwr + h
  mid  <- lwr + (1/2*h)
  
  binned_prop <- function(lwr_x, upr_x, response){
    
    id  <- which(df[,size_var] > lwr_x & df[,size_var] < upr_x) 
    tmp <- df[id,]
    
    if( response == 'prob' ){   return( sum(tmp[,resp_var],na.rm=T) / nrow(tmp) ) }
    if( response == 'n_size' ){ return( nrow(tmp) ) }
    
  }
  
  y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
  x_binned <- mid
  y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist
  
  data.frame( x = x_binned, 
              y = y_binned )
       
}



# Surival

# store 100 posterior sample predictions
post_pred <- function(ii, tmp_anom ){
  y_raw    <- (all_pars$a_u_s[ii]*2) + 
              (all_pars$b_u_s[ii]*2) * x_seq +
               all_pars$b_s2[ii] * (x_seq^2) +
               all_pars$b_s3[ii] * (x_seq^3) +
               all_pars$b_c_s[ii] * tmp_anom
  
  data.frame( x     = x_seq,
              y_hat = boot::inv.logit(y_raw),
              id    = ii )
}

# hot and cold temperature
pred_hot <- lapply(1:100,post_pred, 5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              
pred_cold <- lapply(1:100,post_pred, -5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              
# plot it all out
binned_df <- df_binned_prop(surv, 10, log_area_t0, surv_t1) 
  
p1 <- ggplot( pred_hot ) +
  geom_ribbon( aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.5,
               fill = 'red'
              ) + 
  ylim( 0, 1 ) +
  geom_ribbon( data = pred_cold,
               aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.5,
               fill = 'blue'
              ) + 
  geom_point( data = binned_df, 
              aes(x=x,
                  y=y) 
              ) + 
  labs( y = 'Proportion surviving',
        x = expression('log(size)'[t]) ) 




# Growth

# store 100 posterior sample predictions
post_pred <- function(ii, tmp_anom ){
  y_hat    <- (all_pars$a_u_g[ii]*2) + 
              (all_pars$b_u_g[ii]*2) * x_seq
  data.frame( x     = x_seq,
              y_hat = y_hat,
              id    = ii )
}

# hot and cold temperature
pred_hot <- lapply(1:100,post_pred, 5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              
pred_cold <- lapply(1:100,post_pred, -5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              

p2 <- ggplot( grow ) +
  geom_point( aes( x = log_area_t0,
                   y = log_area_t1 ),
              alpha = 0.2 ) +
  geom_ribbon( data = pred_hot,
               aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.8,
               fill = 'grey'
              ) +
  labs( y = expression('log(size)'[t+1]),
        x = expression('log(size)'[t]) ) 



  
# Flowering

# store 100 posterior sample predictions
post_pred <- function(ii, tmp_anom ){
  y_raw    <- (all_pars$a_u_f[ii]*2) + 
              (all_pars$b_u_f[ii]*2) * x_seq +
               all_pars$b_c_f[ii] * tmp_anom
  
  data.frame( x     = x_seq,
              y_hat = boot::inv.logit(y_raw),
              id    = ii )
}

# hot and cold temperature
pred_hot <- lapply(1:100,post_pred, 5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              
pred_cold <- lapply(1:100,post_pred, -5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              
# plot it all out
binned_df <- df_binned_prop(flow, 10, log_area_t0, flow_t0) 
  
p3 <- ggplot( pred_hot ) +
  geom_ribbon( aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.5,
               fill = 'red'
              ) + 
  ylim( 0, 1 ) +
  geom_ribbon( data = pred_cold,
               aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.5,
               fill = 'blue'
              ) + 
  geom_point( data = binned_df, 
              aes(x=x,
                  y=y) 
              ) + 
  labs( y = 'Proportion flowering',
        x = expression('log(size)'[t]) ) 



# Fertility

# store 100 posterior sample predictions
post_pred <- function(ii, tmp_anom ){
  y_raw    <- (all_pars$a_u_r[ii]*2) + 
              (all_pars$b_u_r[ii]*2) * x_seq +
               all_pars$b_c_r[ii] * tmp_anom
  
  data.frame( x     = x_seq,
              y_hat = exp(y_raw),
              id    = ii )
}

# hot and cold temperature
pred_hot <- lapply(1:100,post_pred, 5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              
pred_cold <- lapply(1:100,post_pred, -5 ) %>% 
                bind_rows %>% 
                group_by( x ) %>% 
                summarise( y_min = quantile(y_hat, prob = 0.1),
                           y_max = quantile(y_hat, prob = 0.9) )
              

p4 <- ggplot( fert ) +
  geom_point( aes( x = log_area_t0,
                   y = numrac_t0 ),
              alpha = 0.2 ) + 
  geom_ribbon( data = pred_hot,
               aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.5,
               fill = 'red'
              ) + 
  geom_ribbon( data = pred_cold,
               aes( x    = x, 
                    ymin = y_min,
                    ymax = y_max ),
               alpha = 0.5,
               fill = 'blue'
              ) + 
  labs( y = 'Proportion flowering',
        x = expression('log(size)'[t]) ) 


# save graph!
out <- gridExtra::grid.arrange(p1,p2,p3,p4, ncol=2) 

ggsave('results/vital_rates/bayes/vr_posterior_ribbon.tiff',
        out, width=6.3, height=6.3, compression='lzw')



# mod_s    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_t0 + 
#                   (1 | year) +     (0 + log_area_t0 | year) +
#                   (1 | location) + (0 + log_area_t0 | location), 
#                   data = surv, family='binomial')
# mod_s    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
#                   (1 | year) +     (log_area_t0 | year) +
#                   (1 | location) + (log_area_t0 | location), 
#                   data = surv, family='binomial')
# 
# 
# 
# coefs <- mod_s %>% fixef
# 
# par(mfrow=c(1,1))
# x_seq    <- seq( min(surv$log_area_t0),
#                  max(surv$log_area_t0), length.out=20)
# plot_binned_prop(surv, 10, log_area_t0, surv_t1)
# 
# y_raw    <- coefs[1] + 
#             coefs[2] * (x_seq) +
#             coefs[3] * (x_seq^2) +
#             coefs[4] * (x_seq^3) +
#             coefs[5] * mean(surv$tmp_t0)
# lines(x_seq, boot::inv.logit(y_raw) )


