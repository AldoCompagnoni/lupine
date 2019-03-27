rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(lme4)
library(bbmle)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------
seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2,
                          log_area_t03= log_area_t0^3 )

# write out seedling size 
sl_size <- data.frame( mean_sl_size = mean(seedl$log_area_t0),
                       sd_sl_size   = sd(seedl$log_area_t0),
                       max_sl_size  = max(seedl$log_area_t0),
                       min_sl_size  = min(seedl$log_area_t0),
                       stringsAsFactors = F)

write.csv(sl_size, 'results/ml_mod_sel/size_sl/seedl_size.csv', row.names=F)

# climate format ----------------------------------------------------------------
years     <- unique(seedl$year)
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
sl_clim  <- left_join(seedl, clim_mat) %>%
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

  # random slope no quadratic
  surv_t1 ~ log_area_t0 + (log_area_t0 | year),
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 + (log_area_t0 | year) + (1 | newid),

  # year + location random slope
  surv_t1 ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 +  (1 | year) + (log_area_t0 | location),
  surv_t1 ~ log_area_t0 +  (log_area_t0 | year) + (log_area_t0 | location) + (1 | newid),
  surv_t1 ~ log_area_t0 +  (1 | year) + (log_area_t0 | location) + (1 | newid)
  
)

# fit models
mods <- lapply( structure_mods,
                function(x) glmer(x, data=sl_clim, family='binomial') ) %>% 
          setNames( c( 'yr', 'yr_loc', 'yr_loc_id', 
                       'yr_id', 'loc', 'id', 'id_loc',
                       'yr_rs', 'yr_loc_rs', 'yr_loc_id_rs', 
                       'yr_id_rs',
                       'yr_rsyl','yr_rsy', 'yr_id_rsyl','yr_id_rsy') )

# best model has YEAR, LOCATION, and log_area_t02
AICtab(mods, weights=T)


# fit climate models -------------------------------------------------
climate_mods <- list(

  # null 
  surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year) + (1 | location),
  
  # precipitation
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_tm1 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0_tm1 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + ppt_t0_tm2 + (log_area_t0 | year) + (1 | location),
  
  # temperature
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_tm1 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0_tm1 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + tmp_t0_tm2 + (log_area_t0 | year) + (1 | location),
  
  # enso
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_t0 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_tm1 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_t0_tm1 + (log_area_t0 | year) + (1 | location),
  surv_t1 ~ log_area_t0 + log_area_t02 + oni_t0_tm2 + (log_area_t0 | year) + (1 | location)
  
)

# fit models
mod_clim <- lapply( climate_mods,
                function(x) glmer(x, data = sl_clim, family='binomial') ) %>% 
              setNames( c( 'null',
                           'ppt_t0', 'ppt_tm1', 'ppt_t0_tm1', 'ppt_t0_tm2',
                           'tmp_t0', 'tmp_tm1', 'tmp_t0_tm1', 'tmp_t0_tm2',
                           'oni_t0', 'oni_tm1', 'oni_t0_tm1', 'oni_t0_tm2') )

# out data frame
out <- data.frame( model  = AICtab(mod_clim, weights=T) %>% attributes %>% .$row.names,
                   dAIC   = AICtab(mod_clim, weights=T) %>% .$dAIC,
                   df     = AICtab(mod_clim, weights=T) %>% .$df,
                   weight = AICtab(mod_clim, weights=T) %>% .$weight )

write.csv(out, 'results/ml_mod_sel/surv_sl/surv_sl_mod_sel.csv', row.names = F)

# fit best model with all data

# tmp_tm1
best_mod <- mod_clim[[out$model[1]]]
best_mod <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + tmp_tm1 +
                 (log_area_t0 | year) + (1 | location),
                data = sl_clim, family='binomial')

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
          'results/ml_mod_sel/surv_sl/surv_sl_best_mod.csv',
          row.names=F)

# plot -----------------------------------------------------------
coefs       <- best_mod %>% fixef
coefs       <- glm(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + tmp_tm1,
                     data = sl_clim, family='binomial') %>% coef

# quantiles of climate predictor
clim_quant  <- sl_clim$tmp_tm1 %>% unique %>% quantile

x_seq <- seq( min(sl_clim$log_area_t0, na.rm=T),
              max(sl_clim$log_area_t0, na.rm=T),
              length.out = 100 )

y_low  <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * x_seq^3 +
          coefs[5] * clim_quant['25%']

y_high <- coefs[1] + 
          coefs[2] * x_seq + 
          coefs[3] * x_seq^2 +
          coefs[4] * x_seq^3 +
          coefs[5] * clim_quant['75%']


tiff('results/ml_mod_sel/surv_sl.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par( mar = c(3.5,3.5,0.2,0.2), mgp = c(2.1,0.7,0), cex.lab = 2 )
plot_binned_prop(sl_clim, 10, log_area_t0, surv_t1)
lines(x_seq, boot::inv.logit(y_low), lwd=2, col = 'blue')
lines(x_seq, boot::inv.logit(y_high), lwd=2, col = 'red')

legend(2.2,0.6, 
       c('low tmp tm1','high tmp tm1'), lwd = 2,
       col = c('blue','red'), bty = 'n', cex = 2)

dev.off()

# exploratory graphs --------------------------------------------------------------------

# sample sizes
rep_df  <- sl_clim %>% 
              group_by( year, location ) %>% 
              summarise( succ  = sum(surv_t1, na.rm=T),
                         fail  = sum(surv_t1 == 0, na.rm=T) ) %>% 
              ungroup %>% 
              mutate( tot = succ + fail,
                      n_i = round(succ/(succ+fail),1),
                      y = 0.8) %>% 
              as.data.frame %>% 
              select(year,location, n_i, tot, y ) 


# Make the multi-panel plot 
ggplot(data  = sl_clim, 
       aes(x = log_area_t0, 
           y = surv_t1) ) +
  geom_point(alpha = 0.5,
             pch = 1) +
  # Fit linear Regression
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"),
              size = 0.5) +
  # split in panels
  facet_grid(location ~ year) +
  theme_bw() +
  ggtitle("Survival (stage_t0 <> SL, surv_t1 <> NA, area_t0 <> 0)") +
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
  ggsave(filename = "results/ml_mod_sel/surv_sl_checks.png",
         dpi = 300, width = 50, height = 30, units = "cm")



# year pooled data checks -----------------------------------------------------
pooled_df  <- sl_clim %>% 
                  group_by( year ) %>% 
                  summarise( succ  = sum(surv_t1, na.rm=T),
                             fail  = sum(surv_t1 == 0, na.rm=T) ) %>% 
                  ungroup %>% 
                  mutate( tot  = succ + fail,
                          prop = round(succ/(succ+fail),3) ) %>% 
                  as.data.frame %>% 
                  select(year, prop )



pool_df <- left_join(pooled_df,clim_mat) 
  
# plot by variable
plot_by_var <- function(var){
  
  # nonstandard evaluation
  var2 <- rlang::enquo(var)
  
  ggplot(data  = pool_df, 
         aes(x = !! var2, 
             y = prop,
             size = 2) ) +
    ylab('Seedling survival')+
    theme(axis.title=element_text(size=20),
          legend.position="none") +
    geom_point(alpha = 1,
               pch   = 16)
}
   


# install.packages("gridExtra")
library(gridExtra)
p <- grid.arrange(plot_by_var(ppt_t0), plot_by_var(ppt_tm1), 
                  plot_by_var(ppt_t0_tm1), plot_by_var(ppt_t0_tm2), 
                  plot_by_var(tmp_t0), plot_by_var(tmp_tm1), 
                  plot_by_var(tmp_t0_tm1), plot_by_var(tmp_t0_tm2), 
                  plot_by_var(oni_t0), plot_by_var(oni_tm1), 
                  plot_by_var(oni_t0_tm1), plot_by_var(oni_t0_tm2),
                  nrow = 3, ncol=4)

ggsave(filename = "results/ml_mod_sel/surv_sl_tab.png",
       p, 
       dpi = 300, width = 50, height = 30, units = "cm")

  