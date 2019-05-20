rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(bbmle)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------
abor        <- subset(lupine_df, !is.na(flow_t0) & flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove 
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t02 = log(area_t0)^2) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  # only years indicated by Tiffany
                  subset( year %in% c(2010, 2011, 2013:2018) )


# climate format ----------------------------------------------------------------
years     <- unique(abor$year)
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
abor_clim  <- left_join(abor, clim_mat) %>%
                subset( !is.na(location) )

# may/june temp
ppt_month <- subset(clim, clim_var == "ppt") %>%
                prism_clim_form("precip", years, 1, m_obs) %>% 
                setNames( c('year','may_ppt') )
tmp_month <- subset(clim, clim_var == "tmean") %>%
                prism_clim_form("tmean", years, 1, m_obs) %>% 
                setNames( c('year','may_tmp') )
enso_month<- subset(enso, clim_var == 'oni' ) %>%
                month_clim_form('oni', years, 1, m_obs) %>% 
                setNames( c('year','may_oni') )

# put together all climate
clim_mon  <- Reduce( function(...) full_join(...),
                     list(ppt_month, tmp_month, enso_month) ) 
                
abor_clim_mon  <- left_join(abor, clim_mon) %>%
                       subset( !is.na(location) )


# fit "structural" models ----------------------------------------------
structure_mods <- list(

  # no quadratic
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | year),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | year) + (1 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | year) + (1 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | year) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (1 | newid) + (1 | location),

  # quadratic
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | year),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | year) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | newid) + (1 | location),
  
  # random slope no quadratic
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (log_area_t0 | year),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (log_area_t0 | year) + (1 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (log_area_t0 | year) + (1 | newid),

  # random slope quadratic
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (log_area_t0 | year),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | newid),
  
  # year + location random slope
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 +  (1 | year) + (log_area_t0 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 +  (1 | year) + (log_area_t0 | location) + (1 | newid),
  
  # year + location random slope + quadratic effect
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | year) + (log_area_t0 | location),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + log_area_t02 + (1 | year) + (log_area_t0 | location) + (1 | newid)
  
)

# fit models
mods <- lapply( structure_mods,
                function(x) glmer(x, data=abor_clim, family='binomial') ) %>% 
          setNames( c( 'yr', 'yr_loc', 'yr_loc_id', 
                       'yr_id', 'loc', 'id', 'id_loc',
                       'yr_2', 'yr_loc_2', 'yr_loc_id_2',
                       'yr_id_2', 'loc_2', 'id_2', 'id_loc_2',
                       'yr_rs', 'yr_loc_rs', 'yr_loc_id_rs', 
                       'yr_id_rs',
                       'yr_rs_2',        'yr_loc_rs_2', 
                       'yr_loc_id_rs_2', 'yr_id_rs_2',
                       'yr_rsyl',        'yr_rsy', 
                       'yr_id_rsyl',     'yr_id_rsy',
                       'yr_rsyl_2',      'yr_rsy_2', 
                       'yr_id_rsyl_2',   'yr_id_rsy_2' ) )

# best model has YEAR, LOCATION, and log_area_t02
AICtab(mods, weights=T)


# fit climate models ----------------------------------------------
climate_mods <- list(

  # null
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  
  # precipitation
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  
  # temperature
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  
  # enso
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  cbind(numab_t0, numrac_t0-numab_t0) ~ log_area_t0 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data=abor_clim, family='binomial') ) %>% 
              setNames( c( 'null',
                           'ppt_t0', 'ppt_tm1', 
                           'tmp_t0', 'tmp_tm1', 
                           'oni_t0', 'oni_tm1') )
                           # 'ppt_t0_loc', 'ppt_tm1_loc', 
                           # 'tmp_t0_loc', 'tmp_tm1_loc', 
                           # 'oni_t0_loc', 'oni_tm1_loc') )
no_ppt_tm1 <- mod_clim[c(1:2,4:7)]
AICtab(mod_clim, weights=T)
AICtab(no_ppt_tm1, weights=T)
AICtab(no_ppt_tm1, weights=T)

# out mod sel
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

write.csv(out, 'results/ml_mod_sel/abor/abor_mod_sel_loc.csv', row.names=F)

# fit best model with all data
best_mod <- glmer(cbind(numab_t1, numall_t1-numab_t1) ~ log_area_t1 + ppt_tm1 
                  + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
                  data = abor_clim, family='binomial' )

best_mod <- glm(cbind(numab_t1, numrac_t1-numab_t1) ~ 1,
                data = subset(abor_clim, year %in% c(2010,2011, 2013:2018) ),
                family='binomial' )
boot::inv.logit(coef(best_mod))

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
          'results/ml_mod_sel/abor/abor_best_mod.csv',
          row.names=F)


# plot -----------------------------------------------------------

# coefficients
coefs    <- best_mod %>% fixef

# quantiles of climate predictor
clim_quant <- abor_clim$ppt_tm1 %>% unique %>% quantile

x_seq <- seq( min(abor_clim$numall_t1, na.rm=T),
              max(abor_clim$numall_t1, na.rm=T),
              length.out = 100 )

y_low  <- boot::inv.logit(coefs[1] + coefs[3] * clim_quant['25%'])
y_high <- boot::inv.logit(coefs[1] + coefs[3] * clim_quant['75%'])


tiff('results/ml_mod_sel/abor.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot(numab_t1 ~ numall_t1,data = abor_clim,
     ylab = 'n. of flowers aborted',
     xlab = 'n. of flowers')
lines(x_seq, x_seq*y_low, lwd=2, col = 'blue')
lines(x_seq, x_seq*y_high, lwd=2, col = 'red')

legend(0,70, 
       c('low ppt tm1','high ppt tm1'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()



# fit climate models ----------------------------------------------
climate_mods <- list(

  # null
  cbind(numab_t1, numall_t1-numab_t1) ~ log_area_t1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # precipitation
  cbind(numab_t1, numall_t1-numab_t1) ~ log_area_t1 + may_ppt + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # temperature
  cbind(numab_t1, numall_t1-numab_t1) ~ log_area_t1 + may_tmp + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # enso
  cbind(numab_t1, numall_t1-numab_t1) ~ log_area_t1 + may_oni + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data=abor_clim, family='binomial') ) %>% 
              setNames( c( 'null',
                           'may_ppt', 'may_tmp', 
                           'may_oni') )

AICtab(mod_clim, weights=T)
