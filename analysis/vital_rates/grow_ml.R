rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(lme4)
library(bbmle)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------
grow        <- lupine_df %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF", "SL")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2 )
                  

# climate format ----------------------------------------------------------------
years     <- unique(grow$year)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  # set names of climate variables
  clim_names <- paste0( var,c('_t0','_tm1') )
  
  mutate(x, 
         avgt0   = x %>% select(V1:V12) %>% rowSums,
         avgtm1  = x %>% select(V13:V24) %>% rowSums ) %>% 
    select(year, avgt0, avgtm1) %>% 
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
grow_clim  <- left_join(grow, clim_mat) %>%
                subset( !is.na(location) )


# fit "structural" models ----------------------------------------------
structure_mods <- list(

  # no quadratic
  log_area_t1 ~ log_area_t0 + (1 | year),
  log_area_t1 ~ log_area_t0 + (1 | year) + (1 | location),
  log_area_t1 ~ log_area_t0 + (1 | year) + (1 | location) + (1 | newid),
  log_area_t1 ~ log_area_t0 + (1 | year) + (1 | newid),
  log_area_t1 ~ log_area_t0 + (1 | location),
  log_area_t1 ~ log_area_t0 + (1 | newid),
  log_area_t1 ~ log_area_t0 + (1 | newid) + (1 | location),

  # quadratic
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | year),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | location),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | location) + (1 | newid),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | newid),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | location),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | newid),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | newid) + (1 | location),

  # random slope no quadratic
  log_area_t1 ~ log_area_t0 + (log_area_t0 | year),
  log_area_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | location),
  log_area_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  log_area_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | newid),

  # random slope quadratic
  log_area_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | newid),
  
  # year + location random slope
  log_area_t1 ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 +  (1 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  log_area_t1 ~ log_area_t0 +  (1 | year) + (log_area_t0 | location) + (1 | newid),
  
  # year + location random slope + quadratic effect
  log_area_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  log_area_t1 ~ log_area_t0 + log_area_t02 + (1 | year) + (log_area_t0 | location) + (1 | newid)


)

# fit models
mods <- lapply( structure_mods,
                function(x) lmer(x, data=grow_clim ) ) %>%
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
AICtab(mods, weights=T)

# fit climate models ----------------------------------------------
climate_mods <- list(

  # null
  log_area_t1 ~ log_area_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  log_area_t1 ~ log_area_t0 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  log_area_t1 ~ log_area_t0 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  log_area_t1 ~ log_area_t0 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  log_area_t1 ~ log_area_t0 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                function(x) lmer(x, data=grow_clim ) ) %>% 
              setNames( c( 'null',
                           'ppt_t0', 'ppt_tm1', 
                           'tmp_t0', 'tmp_tm1', 
                           'oni_t0', 'oni_tm1') )
                           
AICtab(mod_clim, weights=T)

# out mod sel
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

write.csv(out, 'results/ml_mod_sel/grow/grow_mod_sel_loc.csv', row.names=F)


# fit best model with all data
best_mod <- lmer(log_area_t1 ~ log_area_t0 + (log_area_t0 | year) + (log_area_t0 | location),
                  data = grow_clim )

# Growth variance model?
x       <- fitted(best_mod)
y       <- resid(best_mod)^2
gr_var  <- nls(y~a*exp(b*x),start=list(a=1,b=0))
y_pred  <- coef(gr_var)['a'] * exp(coef(gr_var)['b']*x)
plot(x,y)
lines(x,y_pred)
plot(x,y_pred)

# upper and lower observed sizes
siz_vec <- c( min(grow_clim$log_area_t0, grow_clim$log_area_t1, na.rm=T),
              max(grow_clim$log_area_t0, grow_clim$log_area_t1, na.rm=T) ) %>% 
                setNames( c('min_size','max_size') )

# random effects
re_df   <- ranef(best_mod)$year %>% 
              tibble::add_column(.before=1, coef = row.names(.) ) %>% 
              bind_rows( ranef(best_mod)$location %>% 
                         tibble::add_column(.before=1, coef = row.names(.) )
              ) %>% 
              mutate( type_coef = 'ranef' )

# out (binding fixed effects)
out_df  <- fixef(best_mod) %>% 
              # store also model's SD
              c( setNames(summary(best_mod)$sigma,'sigma') ) %>% 
              # store unequal variance model
              c( coef(gr_var) ) %>% 
              # store min/max sizes
              c( siz_vec ) %>% 
              t %>% t %>% 
              as.data.frame %>% 
              tibble::add_column(.before=1, ranef = row.names(.) ) %>% 
              mutate( type_coef = 'fixef' ) %>% 
              bind_rows( re_df )
     
write.csv(out_df, 
          'results/ml_mod_sel/grow/grow_best_mod.csv',
          row.names=F)


# plot ----------------------------------------------------------------------

coefs <- best_mod %>% fixef

# quantiles of climate predictor
x_seq <- seq( min(grow_clim$log_area_t0, na.rm=T),
              max(grow_clim$log_area_t0, na.rm=T),
              length.out = 100 )
y_prd <- coefs[1] + coefs[2] * x_seq 


tiff('results/ml_mod_sel/grow.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mar = c(3.5,3.5,2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot(log_area_t1 ~ log_area_t0, data=grow_clim,
     ylab = 'log_size_t1',
     main = 'NOTE: best model has no climate..')
lines(x_seq, y_prd, lwd=2)

dev.off()

