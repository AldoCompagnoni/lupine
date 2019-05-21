rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(testthat)
library(bbmle)
library(lme4)
library(readxl)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')



# data format --------------------------------------------------------------
fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log_area_t0_z,
                          log_area_t02 = log_area_t0_z^2 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) 

# climate format ----------------------------------------------------------------
years     <- 1990:2018 #unique(fert$year)
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
              year_anom("ppt", years, m_back, m_obs) #%>% 
              # prism_clim_form("precip", years, m_back, m_obs) %>% 
              # year_anom('ppt')

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              year_anom("tmp", years, m_back, m_obs)
              # prism_clim_form('tmean', years, m_back, m_obs) %>% 
              # year_anom('tmp')
  
enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_m_anom('oni')

spei_mat <- subset(clim, clim_var == 'spei' ) %>%
              spei_clim_form(years, m_back, m_obs) %>% 
              year_m_anom('spei')

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat, tmp_mat, enso_mat, spei_mat) )

# demography plus clim
fert_clim  <- left_join(fert, clim_mat) %>%
                subset( !is.na(location) )


# climate model selection ---------------------------------------
climate_mods <- list(

  # null
  numrac_t0 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  numrac_t0 ~ log_area_t0 + log_area_t02 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  numrac_t0 ~ log_area_t0 + log_area_t02 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  numrac_t0 ~ log_area_t0 + log_area_t02 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location),

  # spei
  numrac_t0 ~ log_area_t0 + log_area_t02 + spei_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + spei_tm1 + (log_area_t0 | year) + (log_area_t0 | location)

)

# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data=fert_clim, family='poisson') ) %>% 
              setNames( c( 'null',
                           'ppt_t0',  'ppt_tm1',  
                           'tmp_t0',  'tmp_tm1',  
                           'oni_t0',  'oni_tm1',  
                           'spei_t0', 'spei_tm1') )

AICtab(mod_clim, weights=T)








# n flowers/individual -----------------------------------------------------------------
f_i_df  <- fert_clim %>% 
              group_by( year,location ) %>% 
              summarise( f_n  = sum(numrac_t1),
                         tot  = n(),
                         oni_t0_tm2 = mean(oni_t0_tm2, na.rm=T) ) %>% 
              ungroup %>% 
              mutate( n_i = f_n/tot ) %>% 
              as.data.frame

par(mfrow=c(1,1),mar=c(2.5,2.5,0.5,0.5))
plot(n_i ~ year, data=f_i_df)
plot(n_i ~ oni_t0_tm2, data=f_i_df,  pch = 16)


# fit "structural" models ----------------------------------------------
structure_mods <- list(

  # no quadratic
  numrac_t1 ~ log_area_t1 + (1 | year),
  numrac_t1 ~ log_area_t1 + (1 | year) + (1 | location),
  numrac_t1 ~ log_area_t1 + (1 | year) + (1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + (1 | year) + (1 | newid),
  numrac_t1 ~ log_area_t1 + (1 | location),
  numrac_t1 ~ log_area_t1 + (1 | newid),
  numrac_t1 ~ log_area_t1 + (1 | newid) + (1 | location),

  # quadratic
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | year),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | year) + (1 | location),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | year) + (1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | year) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | location),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | newid) + (1 | location),
  
  # random slope no quadratic
  numrac_t1 ~ log_area_t1 + (log_area_t1 | year),
  numrac_t1 ~ log_area_t1 + (log_area_t1 | year) + (1 | location),
  numrac_t1 ~ log_area_t1 + (log_area_t1 | year) + (1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + (log_area_t1 | year) + (1 | newid),

  # random slope quadratic
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (1 | location),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (1 | newid),
  
  # year + location random slope
  numrac_t1 ~ log_area_t1 +  (log_area_t1 | year) + (log_area_t1 | location),
  numrac_t1 ~ log_area_t1 +  (1 | year) + (log_area_t1 | location),
  numrac_t1 ~ log_area_t1 +  (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 +  (1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # year + location random slope + quadratic effect
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (log_area_t1 | location),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | year) + (log_area_t1 | location),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + (1 | year) + (log_area_t1 | location) + (1 | newid)
  
)

# fit models
mods <- lapply( structure_mods,
                function(x) glmer(x, data=fert_clim, family='poisson') ) %>% 
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
  numrac_t1 ~ log_area_t1 + log_area_t12 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # precipitation
  numrac_t1 ~ log_area_t1 + log_area_t12 + ppt_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + ppt_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + ppt_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + ppt_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # temperature
  numrac_t1 ~ log_area_t1 + log_area_t12 + tmp_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + tmp_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + tmp_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + tmp_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  
  # enso
  numrac_t1 ~ log_area_t1 + log_area_t12 + oni_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + oni_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + oni_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + oni_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),

  # spei
  numrac_t1 ~ log_area_t1 + log_area_t12 + spei_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + spei_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + spei_t0_tm1 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid),
  numrac_t1 ~ log_area_t1 + log_area_t12 + spei_t0_tm2 + (log_area_t1 | year) + (log_area_t1 | location) + (1 | newid)

)

climate_mods <- list(

  # null
  numrac_t0 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # precipitation
  numrac_t0 ~ log_area_t0 + log_area_t02 + ppt_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + ppt_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # temperature
  numrac_t0 ~ log_area_t0 + log_area_t02 + tmp_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + tmp_tm1 + (log_area_t0 | year) + (log_area_t0 | location),
  
  # enso
  numrac_t0 ~ log_area_t0 + log_area_t02 + oni_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + oni_tm1 + (log_area_t0 | year) + (log_area_t0 | location),

  # spei
  numrac_t0 ~ log_area_t0 + log_area_t02 + spei_t0 + (log_area_t0 | year) + (log_area_t0 | location),
  numrac_t0 ~ log_area_t0 + log_area_t02 + spei_tm1 + (log_area_t0 | year) + (log_area_t0 | location)

)


# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data=fert_clim, family='poisson') ) %>% 
              setNames( c( 'null',
                           'ppt_t0',  'ppt_tm1',  
                           'tmp_t0',  'tmp_tm1',  
                           'oni_t0',  'oni_tm1',  
                           'spei_t0', 'spei_tm1') )

AICtab(mod_clim, weights=T)

# out mod sel
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

out

write.csv(out, 'results/ml_mod_sel/fert/fert_mod_sel.csv', row.names=F)


best_mod <- mod_clim[[out$model[1]]]
best_mod <- glmer(numrac_t1 ~ log_area_t1 + log_area_t12 + tmp_t0 + (log_area_t1 | year) + (log_area_t1 | location) + (1|newid), 
                  data=fert_clim, family='poisson')
          
# fruits per raceme
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        

# seeds per fruit
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))

# germination
germ_coef<- select(germ, g0:g2) %>% colMeans

# random effects
re_df    <- ranef(best_mod)$year %>% 
              tibble::add_column(.before=1, coef = row.names(.) ) %>% 
              bind_rows( ranef(best_mod)$location %>% 
                         tibble::add_column(.before=1, coef = row.names(.) )
              ) %>% 
              mutate( type_coef = 'ranef' )

# out (binding fixed effects)
out_df  <- fixef(best_mod) %>% 
              # include fruit_per_raceme
              c( setNames( coef(fr_rac) %>% exp, 'fruit_per_raceme') ) %>% 
              # include seed_per_fruit
              c( setNames( coef(fr_rac) %>% exp, 'seed_per_fruit') ) %>% 
              c( germ_coef ) %>% 
              t %>% t %>% 
              as.data.frame %>% 
              tibble::add_column(.before=1, ranef = row.names(.) ) %>% 
              mutate( type_coef = 'fixef' ) %>% 
              bind_rows( re_df )
     
write.csv(out_df, 
          'results/ml_mod_sel/fert/fert_best_mod.csv',
          row.names=F)


# plot -----------------------------------------------------------

mod_glm <- glm(numrac_t1 ~ log_area_t1 + log_area_t12 + tmp_, 
               data=fert_clim, family='poisson')
coef_glm<- coef(mod_glm)
x_seq <- seq( min(fert_clim$log_area_t1, na.rm=T),
              max(fert_clim$log_area_t1, na.rm=T),
              length.out = 100 )
y_pred<- coef_glm[1] + 
          coef_glm[2] * x_seq + 
          coef_glm[3] * x_seq^2

par( mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot(numrac_t1 ~ log_area_t1, data=fert_clim,
     ylab = 'n. of flowers')
lines(x_seq, exp(y_pred), lwd=2, col = 'blue')

y_low  <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 

lines(x_seq, exp(y_low), lwd=2, col = 'red')

ggplot(data  = fert_clim, 
       aes(x = log_area_t1, 
           y = numrac_t1) ) +
  geom_point(alpha = 0.5,
             pch = 1) +
  # Fit linear Regression
  geom_smooth(method = "glm", size = 0.5, 
              method.args = list(family = "poisson"))


# 
coefs    <- best_mod %>% fixef

# quantiles of climate predictor
clim_quant <- fert_clim$tmp_t0 %>% unique %>% quantile

x_seq <- seq( min(fert_clim$log_area_t1, na.rm=T),
              max(fert_clim$log_area_t1, na.rm=T),
              length.out = 100 )

y_low  <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * clim_quant['25%']

y_high <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * clim_quant['75%']


tiff('results/ml_mod_sel/fert.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot(numrac_t1 ~ log_area_t1, data=fert_clim,
     ylab = 'n. of flowers')
lines(x_seq, exp(y_low), lwd=2, col = 'blue')
lines(x_seq, exp(y_high), lwd=2, col = 'red')

legend(0,150, 
       c('low tmp t0','high tmp t0'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()


# exploratory graphs --------------------------------------------------------------------

# sample sizes
rep_df  <- fert_clim %>% 
              group_by( year,location ) %>% 
              summarise( f_n  = sum(numrac_t1),
                         tot  = n() ) %>% 
              ungroup %>% 
              mutate( n_i = round(f_n/tot,1),
                      y   = 5 ) %>% 
              as.data.frame %>% 
              select(year,location, n_i, tot, y )
              

fert_plot <- fert_clim %>% 
                mutate( log_numrac_t1 = log(numrac_t1) )

# Make the multi-panel plot 
ggplot(data  = fert_plot, 
       aes(x = log_area_t1, 
           y = log_numrac_t1) ) +
  geom_point(alpha = 0.5,
             pch = 1) +
  # Fit linear Regression
  geom_smooth(method = "lm", size = 0.5) +
  # split in panels
  facet_grid(location ~ year) +
  theme_bw() +
  ggtitle("N. of racemes (stage_t0 <> SL, area_t0 <> 0, numrac_t0 <>0, flow_t0 = 1)") + 
  geom_text(
    data    = rep_df,
    mapping = aes(x = -Inf, y = -Inf, label = n_i),
    hjust   = -0.1,
    vjust   = -1
  ) +
   geom_text(
    data    = rep_df,
    mapping = aes(x = -Inf, y = y, label = tot),
    hjust   = -0.1,
    vjust   = -1
  ) +
  ggsave(filename = "results/ml_mod_sel/fertility_checks.png",
         dpi = 300, width = 50, height = 30, units = "cm")


# year pooled data checks -----------------------------------------------------
pooled_df  <- fert_clim %>% 
                group_by( year ) %>% 
                summarise( f_n  = sum(numrac_t1),
                           tot  = n() ) %>% 
                ungroup %>% 
                mutate( n_i = round(f_n/tot,1) )
  
# pool climate and survival info together
pool_df <- left_join(pooled_df,clim_mat) 

    
# plot by variable
plot_by_var <- function(var){
  
  # nonstandard evaluation
  var2 <- rlang::enquo(var)
  
  ggplot(data  = pool_df, 
         aes(x = !! var2, 
             y = n_i,
             size = 2) ) +
    ylab('Flowers per plant')+
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

ggsave(filename = "results/ml_mod_sel/fert_tab.png",
       p, 
       dpi = 300, width = 50, height = 30, units = "cm")
