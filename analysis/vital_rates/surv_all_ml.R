rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(mgcv)
library(lme4)
library(bbmle)
library(ggplot2)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')

# data
lupine_df   <- read.csv("data/lupine_all.csv")
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")

# data format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  # subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0  = log_area_t0_z ) %>% 
                  mutate( log_area_t02 = log_area_t0_z^2,
                          log_area_t03 = log_area_t0_z^3 )

# climate format ----------------------------------------------------------------
years     <- 2003:2018 #unique(surv$year) 
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies from monthly anomalies
year_m_anom <- function(x, var ){

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
              # prism_clim_form("precip", years, m_back, m_obs) %>% 
              # year_anom('ppt')

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              year_anom("tmp", years, m_back, m_obs) 
              # prism_clim_form('tmean', years, m_back, m_obs) %>% 
              # year_anom('tmp')
  
enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_m_anom('oni')

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat,tmp_mat,enso_mat) )

# demography plus clim
surv_clim <- left_join(surv, clim_mat) %>%
                subset( !is.na(location) )


climate_mods <- list(

  # null
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                     function(x) glmer(x, data=surv_clim, family='binomial') ) %>% 
                setNames( c( 'null',
                             'ppt_t0', 'ppt_tm1', 'tmp_t0', 'tmp_tm1', 
                             'oni_t0', 'oni_tm1') )

AICtab(mod_clim, weights=T)

mod <- glmer(climate_mods[[4]], data=surv_clim, family='binomial') 




# fit "structural" models ----------------------------------------------
structure_mods <- list(

  # no quadratic
  surv_t1 ~ log_area_t0 + (1 | year),
  surv_t1 ~ log_area_t0 + (1 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + (1 | year) + (1 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 + (1 | year) + (1 | newid),
  surv_t1 ~ log_area_t0 + (1 | location),
  surv_t1 ~ log_area_t0 + (1 | newid),
  surv_t1 ~ log_area_t0 + (1 | newid) + (1 | location),

  # quadratic
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | newid),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | newid),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | newid) + (1 | location),
  
  # random slope no quadratic
  surv_t1 ~ log_area_t0 + (log_area_t0 | year),
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | newid),

  # random slope quadratic
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year),
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | newid),

  # year + location random slope
  surv_t1 ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 +  (1 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 +  (1 | year) + (log_area_t0 | location) + (1 | newid),
  
  # year + location random slope + quadratic effect
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (log_area_t0 | location) + (1 | newid)
  
)

# fit models
mods <- lapply( structure_mods,
                function(x) glmer(x, data=surv_clim, family='binomial') ) %>% 
          # id: individual random effect
          
          setNames( c( 'yr', 'yr_loc', 'yr_loc_id', 
                       'yr_id', 'loc', 'id', 'id_loc',
                       'yr_2', 'yr_loc_2', 'yr_loc_id_2',
                       'yr_id_2', 'loc_2', 'id_2', 'id_loc_2',
                       'yr_rs', 'yr_loc_rs', 'yr_loc_id_rs', 
                       'yr_id_rs',
                       'yr_rs_2', 'yr_loc_rs_2', 'yr_loc_id_rs_2',
                       'yr_id_rs_2',
                       'yr_rsyl','yr_rsy', 'yr_id_rsyl','yr_id_rsy',
                       'yr_rsyl_2','yr_rsy_2', 'yr_id_rsyl_2','yr_id_rsy_2' ) )

# best model has YEAR, LOCATION, and log_area_t02
AICtab(mods)

climate_mods <- list(
  
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ (1 | year) + (1 | location)
  
)

mods <- lapply( climate_mods,
                function(x) glmer(x, data=surv_clim, family='binomial') )

AICtab(mods, weights=T)

# fit climate models ----------------------------------------------
climate_mods <- list(

  # null
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + ppt_t0_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + ppt_t0_tm2 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + tmp_t0_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + tmp_t0_tm2 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + oni_t0_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + oni_t0_tm2 + (log_area_t0 | year) + (log_area_t0 | location)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                     function(x) glmer(x, data=surv_clim, family='binomial') ) %>% 
                setNames( c( 'null',
                             'ppt_t0', 'ppt_tm1', 'ppt_t0_tm1', 'ppt_t0_tm2',
                             'tmp_t0', 'tmp_tm1', 'tmp_t0_tm1', 'tmp_t0_tm2',
                             'oni_t0', 'oni_tm1', 'oni_t0_tm1', 'oni_t0_tm2') )

AICtab(mod_clim, weights=T)


climate_mods <- list(

  # null
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                     function(x) glmer(x, data=surv_clim, family='binomial') ) %>% 
                setNames( c( 'null',
                             'ppt_t0', 'ppt_tm1', 'tmp_t0', 'tmp_tm1', 
                             'oni_t0', 'oni_tm1') )


# out mod sel
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

out

write.csv(out, 'results/ml_mod_sel/surv/surv_mod_sel_no2.csv', row.names=F)


# fit best model with all data

# oni_t0_tm1
# best_mod <- mod_clim[[out$model[1]]]
best_mod <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + 
                  tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
                  data=surv_clim, family='binomial')
  

# random effects
re_df   <- ranef(best_mod)$year %>% 
              tibble::add_column(.before=1, coef = row.names(.) ) %>% 
              bind_rows( ranef(best_mod)$location %>% 
                         tibble::add_column(.before=1, coef = row.names(.) )
              ) %>% 
              mutate( type_coef = 'ranef' )

# out (binding fixed effects)
out_df  <- fixef(best_mod) %>% 
              t %>% t %>% 
              as.data.frame %>% 
              tibble::add_column(.before=1, ranef = row.names(.) ) %>% 
              mutate( type_coef = 'fixef' ) %>% 
              bind_rows( re_df )
     
write.csv(out_df, 
          'results/ml_mod_sel/surv/surv_best_mod.csv',
          row.names=F)


# plot ---------------------------------------------------------------

best_mod <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (log_area_t0 | year) + (log_area_t0 | location),
                  data=surv_clim, family='binomial')
  
mod1 <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (1 | year) + (1 | location) + 
                  (0 + log_area_t0 | year) + (0 + log_area_t0 | location),
                  data=surv_clim, family='binomial')

mod2 <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (log_area_t0 | year) + (log_area_t0 | location) + 
                  (0 + log_area_t0 | year) + (0 + log_area_t0 | location) +
                  (0 + log_area_t02 | year) + (0 + log_area_t02 | location),
                  data=surv_clim, family='binomial')

mod3 <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1 + 
                  (log_area_t0 | year) + (log_area_t0 | location) + 
                  (0 + log_area_t0 | year) + (0 + log_area_t0 | location) +
                  (0 + log_area_t02 | year) + (0 + log_area_t02 | location) + 
                  (0 + log_area_t03 | year) + (0 + log_area_t03 | location),
                  data=surv_clim, family='binomial')

AIC(best_mod, mod1, mod2, mod3)


coefs  <- fixef(best_mod)

# quantiles of climate predictor
clim_quant  <- surv_clim$tmp_tm1 %>% unique %>% quantile

x_seq  <- seq( min(surv_clim$log_area_t0, na.rm=T),
               max(surv_clim$log_area_t0, na.rm=T),
               length.out = 100 )

y_low  <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * (x_seq^2) + 
          coefs[4] * (x_seq^3) + 
          coefs[5] * clim_quant['25%']
y_high <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * (x_seq^2) +
          coefs[4] * (x_seq^3) +
          coefs[5] * clim_quant['75%']



tiff('results/ml_mod_sel/surv_no2.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mfrow=c(1,1), mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot_binned_prop(surv_clim, 8, log_area_t0, surv_t1)
lines(x_seq, boot::inv.logit(y_low), lwd=2, col = 'blue')
lines(x_seq, boot::inv.logit(y_high), lwd=2, col = 'red')

legend(2.2,0.6, 
       c('low tmp tm1','high tmp tm1'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()


# year-by-site plots ----------------------------------------------------------------


# data frame of binned proportions
df_binned_prop <- function(ii, df_in, n_bins, siz_var, rsp_var, s_y_df){
  
  # make sub-selection of data
  df   <- subset(df_in, year     == s_y_df$year[ii] & 
                        location == s_y_df$location[ii] )
  
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
    mutate( year     = s_y_df$year[ii], 
            location = s_y_df$location[ii] )
  
}

# grid of year/location
s_y_df      <- expand.grid( year     = surv_clim$year %>% unique %>% sort,
                            location = surv_clim$location %>% unique %>% sort,
                            stringsAsFactors = F)

# produce data frames 
surv_bin_l  <- lapply(1:nrow(s_y_df), df_binned_prop, surv_clim, 10, 
                                      log_area_t0, surv_t1, s_y_df)

# big data frame
surv_bin_df <- bind_rows( surv_bin_l ) %>% 
                  mutate( log_area_t02 = log_area_t0^2, 
                          log_area_t03 = log_area_t0^3,
                          tmp_tm1      = 0,
                          transition   = paste( paste0(year - 1), 
                                                substr(paste0(year),3,4),
                                                sep='-') )

# include predicted values
surv_bin_df <- surv_bin_df %>% 
                  mutate( yhat         = predict(best_mod, 
                                                 newdata=surv_bin_df, 
                                                 type='response') )


# plot it all out
ggplot(data  = surv_bin_df, 
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
  facet_grid(location ~ transition) +
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
  ggsave(filename = "results/vital_rates/surv_all.tiff",
         dpi = 300, width = 6.3, height = 4, units = "in",
         compression = 'lzw')


# # exploratory graphs --------------------------------------------------------------------
# 
# # sample sizes
# rep_df  <- surv_clim %>% 
#               group_by( year, location ) %>% 
#               summarise( succ  = sum(surv_t1, na.rm=T),
#                          fail  = sum(surv_t1 == 0, na.rm=T) ) %>% 
#               ungroup %>% 
#               mutate( tot = succ + fail,
#                       n_i = round(succ/(succ+fail),1),
#                       y = 0.8) %>% 
#               as.data.frame %>% 
#               select(year,location, n_i, tot, y ) 
# 
# 
# # Make the multi-panel plot 
# ggplot(data  = surv_clim, 
#        aes(x = log_area_t0, 
#            y = surv_t1) ) +
#   geom_point(alpha = 0.5,
#              pch = 1) +
#   # Fit linear Regression
#   geom_smooth(method = "glm", 
#               method.args = list(family = "binomial"),
#               size = 0.5) +
#   # split in panels
#   facet_grid(location ~ year) +
#   theme_bw() +
#   ggtitle("Survival (stage_t0 <> SL, surv_t1 <> NA, area_t0 <> 0)") +
#   geom_text(
#     data    = rep_df,
#     mapping = aes(x = -Inf, y = -Inf, label = n_i),
#     hjust   = -0.1,
#     vjust   = -1
#   ) +
#    geom_text(
#     data    = rep_df,
#     mapping = aes(x = -Inf, y = y, label = tot),
#     hjust   = -0.1,
#     vjust   = -1
#   ) +
#   ggsave(filename = "results/ml_mod_sel/surv_checks.png",
#          dpi = 300, width = 50, height = 30, units = "cm")
  

# year pooled data checks -----------------------------------------------------
pooled_df  <- surv_clim %>% 
                  group_by( year ) %>% 
                  summarise( succ  = sum(surv_t1, na.rm=T),
                             fail  = sum(surv_t1 == 0, na.rm=T) ) %>% 
                  ungroup %>% 
                  mutate( tot  = succ + fail,
                          prop = round(succ/(succ+fail),3) ) %>% 
                  as.data.frame %>% 
                  select(year, prop )

# pool climate and survival info together
pool_df <- left_join(pooled_df,clim_mat) 
  

    
# plot by variable
plot_by_var <- function(var){
  
  # nonstandard evaluation
  var2 <- rlang::enquo(var)
  
  ggplot(data  = pool_df, 
         aes(x = !! var2, 
             y = prop,
             size = 2) ) +
    ylab('Plant survival')+
    theme(axis.title=element_text(size=20),
          legend.position="none") +
    geom_point(alpha = 1,
               pch   = 16)
}
   
# plot it all out 
library(gridExtra)
p <- grid.arrange(plot_by_var(ppt_t0), plot_by_var(ppt_tm1), 
                  plot_by_var(ppt_t0_tm1), plot_by_var(ppt_t0_tm2), 
                  plot_by_var(tmp_t0), plot_by_var(tmp_tm1), 
                  plot_by_var(tmp_t0_tm1), plot_by_var(tmp_t0_tm2), 
                  plot_by_var(oni_t0), plot_by_var(oni_tm1), 
                  plot_by_var(oni_t0_tm1), plot_by_var(oni_t0_tm2),
                  nrow = 3, ncol=4)

ggsave(filename = "results/ml_mod_sel/surv_tab.png",
       p, 
       dpi = 300, width = 50, height = 30, units = "cm")
