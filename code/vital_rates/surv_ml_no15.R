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
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2 ) %>% 
                  mutate( log_area_t0  = scale(log_area_t0) %>% as.vector,
                          log_area_t02 = scale(log_area_t02) %>% as.vector ) %>% 
                  subset( !(year %in% 2015) )
  
# climate format ----------------------------------------------------------------
years     <- unique(surv$year)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
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


# demography plus clim
surv_clim <- left_join(surv, clim_mat) %>%
                subset( !is.na(location) )


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


# fit climate models ----------------------------------------------

climate_mods <- list(

  # null
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0_tm2 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0_tm2 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location)
  
)

# fit models
mods_clim <- lapply( climate_mods,
                     function(x) glmer(x, data=surv_clim, family='binomial') ) %>% 
                setNames( c( 'null',
                             'ppt_t0', 'ppt_tm1', 
                             'tmp_t0', 'tmp_tm1', 
                             'oni_t0', 'oni_tm1') )

AICtab(mods_clim, weights=T)

# out mod sel
out <- data.frame( model  = AICtab(mods, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mods, weights=T) %>% .$dAIC,
                   df     = AICtab(mods, weights=T) %>% .$df,
                   weight = AICtab(mods, weights=T) %>% .$weight )

write.csv(out, 'results/ml_mod_sel/surv/surv_mod_sel.csv', row.names=F)


# fit best model with all data
best_mod <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + tmp_tm1 + 
                    (log_area_t0 | year) + (1 | location),
                  data = surv_clim, family='binomial')

best_mod1 <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + tmp_tm1 + 
                    (log_area_t0 | year) + (log_area_t0 | location),
                  data = surv_clim, family='binomial')


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

# restrict size structure
# small_df <- subset(surv_clim, log_area_t0 < 4)
# 
# mod1 <- glm( surv_t1 ~ 1, data = small_df,
#              family='binomial')
# mod  <- glm( surv_t1 ~ log_area_t0, data = small_df,
#              family='binomial')
# 
# # 
# boot::inv.logit(coef(mod1))
# 
# coef(mod)
# 
# plot(jitter(surv_t1) ~ log_area_t0, data = small_df)
# x_seq <- seq( min(small_df$log_area_t0, na.rm=T),
#               max(small_df$log_area_t0, na.rm=T),
#               length.out = 100 )
# y_low  <- coef(mod)[1] + coef(mod)[2] * x_seq
# lines(x_seq, boot::inv.logit(y_low))
# abline( h= boot::inv.logit(coef(mod1)) )




# binned survival probabilities
h    <- (max(surv_clim$log_area_t0,na.rm=T) - min(surv_clim$log_area_t0,na.rm=T)) / 15
lwr  <- min(surv_clim$log_area_t0,na.rm=T) + (h*c(0:14))
upr  <- max(surv_clim$log_area_t0,na.rm=T) + (h*c(0:14))
mid  <- lwr + (1/2*h)

binned_prop <- function(lwr_x, upr_x, response){
  
  tmp <- subset(surv_clim, log_area_t0 > lwr_x & log_area_t0 < upr_x ) 
  
  if( response == 'prob' ){   return( sum(tmp$surv_t1) / nrow(tmp) ) }
  if( response == 'n_size' ){ return( nrow(tmp) ) }
  
}

y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
x_binned <- mid
y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist


mod_smpl <- glm( surv_t1 ~ log_area_t0, data = surv_clim, family='binomial' )
coef_smpl<- coef(mod_smpl) 

coefs       <- best_mod %>% fixef
# quantiles of climate predictor
clim_quant  <- surv_clim$tmp_tm1 %>% unique %>% quantile

x_seq <- seq( min(surv_clim$log_area_t0, na.rm=T),
              max(surv_clim$log_area_t0, na.rm=T),
              length.out = 100 )

y_low  <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * clim_quant['25%']

y_high <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * clim_quant['75%']


y_smpl <- boot::inv.logit(coef_smpl[1] + coef_smpl[2] * x_seq)


tiff('results/ml_mod_sel/surv.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mfrow=c(2,1), mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
# plot(jitter(surv_t1) ~ log_area_t0, data = surv_clim,
#      ylab = 'survival probability')

plot(x_binned, y_binned,
     ylab = 'survival probability',
     ylim = c(0,1) )
lines(x_seq, boot::inv.logit(y_low), lwd=2, col = 'blue')
lines(x_seq, boot::inv.logit(y_high), lwd=2, col = 'red')

plot(x_binned, y_binned,
     ylab = 'survival probability',
     ylim = c(0,1) )
lines(x_seq, y_smpl)

legend(2.2,0.6, 
       c('low tmp tm1','high tmp tm1'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()



ggplot(data = surv_clim, 
       aes(x = log_area_t0, 
           y = surv_t1)) +
geom_jitter(height = 0.05, 
            alpha = 0.5, 
            pch = 1) +
# Fit Logistic Regression - glm with family = "binomial"
geom_smooth(method = "glm", 
            method.args = list(family = "binomial")) +
# split in panels
theme_bw() +
ggtitle("Survival (stage_t0 <> SL, surv_t1 <> NA, area_t0 <> 0)")



# Spline
mod <- gam( surv_t1 ~ s(log_area_t0, k = 10), 
            family = binomial,  data = surv_clim )
plot(mod)
