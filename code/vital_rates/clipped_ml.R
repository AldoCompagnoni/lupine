rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
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
clip        <- subset(lupine_df, !is.na(flow_t1) & flow_t1 == 1 ) %>% 
                  subset( area_t1 != 0) %>% 
                  subset( !is.na(numrac_t1) ) %>% 
                  # remove 
                  subset( !(flow_t1 %in% 0) ) %>% 
                  mutate( log_area_t1  = log(area_t1),
                  log_area_t12 = log(area_t1)^2,
                  year         = year + 1 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  subset( !(numrac_t1 %in% 0) )


# climate format ----------------------------------------------------------------
years     <- unique(clip$year)
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
clip_clim  <- left_join(clip, clim_mat) %>%
                subset( !is.na(location) )


# fit "structural" models ----------------------------------------------
structure_mods <- list(

  # NULLs
  cbind(numcl_t1, numint_t1) ~ (1 | year),
  cbind(numcl_t1, numint_t1) ~ (1 | year) + (1 | location),
  cbind(numcl_t1, numint_t1) ~ (1 | year) + (1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ (1 | year) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ (1 | location),
  cbind(numcl_t1, numint_t1) ~ (1 | newid),
  cbind(numcl_t1, numint_t1) ~ (1 | newid) + (1 | location),

  # no quadratic
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | year),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | year) + (1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | year) + (1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | year) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (1 | newid) + (1 | location),

  # quadratic
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | year),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | year) + (1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | year) + (1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | year) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | newid) + (1 | location),
  
  # random slope no quadratic
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (log_area_t1 | year),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (log_area_t1 | year) + (1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (log_area_t1 | year) + (1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (log_area_t1 | year) + (1 | newid),

  # random slope quadratic
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (log_area_t1 | year),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (1 | newid),
  
  # year + location random slope
  cbind(numcl_t1, numint_t1) ~ log_area_t1 +  (log_area_t1 | year) + (log_area_t1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 +  (1 | year) + (log_area_t1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 +  (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 +  (1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # year + location random slope + quadratic effect
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (log_area_t1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | year) + (log_area_t1 | location),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + log_area_t12 + (1 | year) + (log_area_t1 | location) + (1 | newid)
  
)

# fit models
mods <- lapply( structure_mods,
                function(x) glmer(x, data=clip_clim, family='binomial') ) %>% 
          setNames( c( 'yr0', 'yr_loc0', 'yr_loc_id0',
                       'yr_id0', 'loc0', 'id0', 'id_loc0',
                       'yr', 'yr_loc', 'yr_loc_id', 
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
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # precipitation
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + ppt_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + ppt_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # temperature
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + tmp_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + tmp_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # enso
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + oni_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  cbind(numcl_t1, numint_t1) ~ log_area_t1 + oni_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid)

)

# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data=clip_clim, family='binomial') ) %>% 
              setNames( c( 'null',
                           'ppt_t0', 'ppt_tm1', 
                           'tmp_t0', 'tmp_tm1', 
                           'oni_t0', 'oni_tm1') )
                           # 'ppt_t0_loc', 'ppt_tm1_loc', 
                           # 'tmp_t0_loc', 'tmp_tm1_loc', 
                           # 'oni_t0_loc', 'oni_tm1_loc') )

AICtab(mod_clim, weights=T)

# out mod sel
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

write.csv(out, 'results/ml_mod_sel/clip/clip_mod_sel_loc.csv', row.names=F)

# fit best model with all data
best_mod <- glmer(cbind(numcl_t1, numall_t1-numcl_t1) ~ log_area_t1 
                  + tmp_t0 + (log_area_t1 | year) + (log_area_t1 | location) 
                  + (1 | newid),
                  data = clip_clim, family='binomial' )

best_mod <- glm(cbind(numcl_t1, notab_t1-numcl_t1) ~ 1,
                data = subset(clip_clim, !(location %in% 'AL (1)')),
                family='binomial')

coef(best_mod) %>% boot::inv.logit()

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
          'results/ml_mod_sel/clip/clip_best_mod.csv',
          row.names=F)


# plot -----------------------------------------------------------

coefs      <- best_mod %>% fixef

# quantiles of climate predictor
clim_quant <- clip_clim$tmp_t0 %>% unique %>% quantile

x_seq <- seq( min(clip_clim$numall_t1, na.rm=T),
              max(clip_clim$numall_t1, na.rm=T),
              length.out = 100 )

y_low  <- boot::inv.logit(coefs[1] + coefs[3] * clim_quant['25%'])
y_high <- boot::inv.logit(coefs[1] + coefs[3] * clim_quant['75%'])


tiff('results/ml_mod_sel/clip.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot(numcl_t1 ~ numall_t1, data=clip_clim,
     ylab = '% clipped')
lines(x_seq, x_seq*y_low, lwd=2, col = 'blue')
lines(x_seq, x_seq*y_high, lwd=2, col = 'red')

legend(0,90, 
       c('low tmp t0','high tmp t0'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()

