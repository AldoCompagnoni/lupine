rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(lme4)
library(bbmle)
library(gridExtra)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# climate format ----------------------------------------------------------------
lupine_df$transition
# 4.19.2018, 2019 data not yet in
years     <- 1990:2018
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
                    list(ppt_mat,tmp_mat,enso_mat) ) %>% 
              gather(measure, anomaly, ppt_t0:oni_t0_tm2) 


# climate anomalies by year ------------------------
clim_mat %>% 
  subset( grepl('_t0$',measure) ) %>% 
  ggplot( aes(x=year, 
                     y=anomaly) ) +
  geom_line( aes(color=measure),
             lwd = 1 ) + 
  scale_color_viridis_d() +
  ylab( 'Anomaly' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2 ) + 
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size  = 20,
                                    angle = 90) ) + 
  ggsave('results/climate/anomalies_vs_year.tiff',
         width=6.3,height=4,compression='lzw')
  

# histograms of anomalies -------------------------------------
p1 <- clim_mat %>%
        subset( grepl('_t0$',measure) ) %>% 
        ggplot( ) + 
        geom_histogram( aes(anomaly),
                        binwidth = 3 ) + 
        facet_grid( ~ measure) + 
        xlab( '' ) +
        ylab( '' ) +
        theme(plot.margin = unit(c(0,0,0,0), "cm"),
              panel.spacing = unit('0.5',unit='mm'),
              strip.text.x  = element_text( size = 10,
                                       margin = margin(0,0,0,0,
                                                       'mm') ) )
  
p2 <- clim_mat %>%
        subset( grepl('_t0$',measure) ) %>%
        subset( year > 2005 ) %>% 
        ggplot( ) + 
        geom_histogram( aes(anomaly),
                        binwidth = 3 ) + 
        facet_grid( ~ measure) + 
        xlab( '' ) + 
        theme(plot.margin = unit(c(0,0,0,0), "cm"),
              panel.spacing = unit('0.5',unit='mm'),
              strip.text.x  = element_text( size = 10,
                                       margin = margin(0,0,0,0,
                                                       'mm') ) )

p3 <- clim_mat %>%
        subset( grepl('_t0$',measure) ) %>%
        subset( year < 2006 ) %>% 
        ggplot( ) + 
        geom_histogram( aes(anomaly),
                        binwidth = 3 ) + 
        facet_grid( ~ measure) +
        ylab( '' ) +
        theme(plot.margin = unit(c(0,0,0,0), "cm"),
              panel.spacing = unit('0.5',unit='mm'),
              strip.text.x  = element_text( size = 10,
                                       margin = margin(0,0,0,0,
                                                       'mm') ) )

g <- arrangeGrob(p1,p2,p3,nrow=3,ncol=1)
ggsave(filename = 'results/climate/anom_hist.tiff',
       plot = g,
       dpi = 300, width = 6.3, height = 4, units = "in",
       compression = 'lzw')
