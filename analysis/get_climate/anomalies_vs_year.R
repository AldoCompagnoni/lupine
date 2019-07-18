rm(list=ls())
library(dplyr)
library(tidyr)
library(testthat)
library(rstan)
library(rstanarm)
library(lme4)
library(ggplot2)
library(ggthemes)
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
m_back    <- 12


# calculate yearly anomalies summing up  monthly anomalies --------
month_anom <- function(x, var){
  
  # set names of climate variables
  clim_names <- paste0( var,'_t0' )
  
  x %>% 
    mutate(sumt0 = x %>% select(V1:V12) %>% rowSums) %>% 
    select(year, sumt0) %>% 
    setNames( c('year',clim_names) )
  
}

ppt_mo <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs) %>% 
              month_anom('ppt')

tmp_mo <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              month_anom('tmp')
  
enso_mo <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              month_anom('oni')

spei_mo <- subset(clim, clim_var == 'spei' ) %>%
              spei_clim_form(years, m_back, m_obs) %>% 
              month_anom('spei')

# put together all climate
clim_mo <- Reduce( function(...) full_join(...),
                   list(ppt_mo,tmp_mo,enso_mo,spei_mo) ) %>% 
              gather(measure, value, ppt_t0:spei_t0)


# (summed) monthly weather anomalies ------------------------------------
yr_anom <- function(clim_x, clim_var = "ppt", 
                    yearz, m_back, m_obs, anom){
  
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
  clim_x_l  <- lapply(yearz, clim_back, m_obs, clim_m)
  x_clim    <- year_by_month_mat(clim_x_l, yearz) %>% 
                  gather(month,value,V1:V12)
  
  if( clim_var == 'ppt'){
    out_raw <- x_clim %>% 
                group_by(year) %>% 
                summarise( ppt_raw = sum(value) )
  }
  if( clim_var == 'tmp'){
    out_raw <- x_clim %>% 
                group_by(year) %>% 
                summarise( tmp_raw = mean(value) )
  }
  
  if( anom ){
    out_raw[,2] <- scale(out_raw[,2], center = T, scale = T)[,1]
    out_raw <- out_raw %>%
                setNames( c('year', 
                            paste0(clim_var,'_anom')) )
  }
  
  out_raw
  
}

# format climate - need to select climate predictor first 
ppt_raw <- subset(clim, clim_var == "ppt") %>%
              yr_anom("ppt", years, m_back, m_obs, F)

tmp_raw <- subset(clim, clim_var == 'tmean') %>% 
              yr_anom("tmp", years, m_back, m_obs, F) 

# put together all climate
clim_raw <- Reduce( function(...) full_join(...),
                    list(ppt_raw,tmp_raw) ) %>% 
              gather(measure, value, ppt_raw:tmp_raw)


# yearly anomalies  -------------------------------------------

# format climate 
ppt_yr <- subset(clim, clim_var == "ppt") %>%
              yr_anom("ppt", years, m_back, m_obs, T)

tmp_yr <- subset(clim, clim_var == 'tmean') %>% 
              yr_anom("tmp", years, m_back, m_obs, T) 

oni_yr <- enso_mo %>% 
              mutate( oni_anom  = oni_t0/12 )

spei_yr<- spei_mo %>% 
              mutate( spei_anom = spei_t0/12 )

# put together all climate
clim_yr <- Reduce( function(...) full_join(...),
                   list(ppt_yr, tmp_yr, 
                        oni_yr, spei_yr) ) %>% 
              dplyr::select(-oni_t0, -spei_t0 ) %>% 
              gather(measure, value, ppt_anom:spei_anom)

# Plots ----------------------------------------------


# all monthly climate anomalies by year
clim_mo %>% 
  subset( grepl('_t0$',measure) ) %>% 
  ggplot( aes(x=year, 
              y=value) ) +
  geom_line( aes(color=measure),
             lwd = 2 ) + 
  scale_color_viridis_d() +
  ylab( 'Anomaly' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2) + 
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size  = 20,
                                    angle = 0) ) + 
  ggsave('results/climate/anomalies_vs_year.tiff',
         width=12,height=5,compression='lzw')


clim_mo %>% 
  subset( measure != 'oni_t0' ) %>% 
  subset( grepl('_t0$',measure) ) %>%
  mutate( measure = replace(measure, measure=='tmp_t0',
                            'Temperature') ) %>% 
  mutate( measure = replace(measure, measure=='ppt_t0',
                            'Precipitation') ) %>% 
  mutate( measure = replace(measure, measure=='spei_t0',
                            'Aridity I.') ) %>% 
  ggplot( aes(x=year, 
              y=value) ) +
  geom_line( aes(color=measure),
             lwd = 2 ) + 
  scale_color_viridis_d() +
  ylab( 'Anomaly' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2) + 
  theme( axis.title = element_text( size  = 20),
         axis.text  = element_text( size  = 20,
                                    angle = 0) ) + 
  ggsave('results/climate/anomalies_vs_year_noOni.tiff',
         width=12,height=5,compression='lzw')
  


# Raw time series versus anomalies ---------------------
p1 <- clim_mo %>% 
  subset( grepl('ppt_|tmp_',measure) ) %>% 
  ggplot( aes(x=year, 
              y=value) ) +
  geom_line( aes(color=measure),
             lwd = 1 ) + 
  scale_color_viridis_d() +
  ylab( 'Monthly anomalies' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2 ) + 
  theme( axis.title = element_text( size = 10),
         axis.text  = element_text( size  = 10,
                                    angle = 90) )

p2 <- clim_yr %>% 
        ggplot( aes(x=year, 
                    y=value) ) +
        geom_line( aes(color=measure),
                   lwd = 1 ) + 
        scale_color_viridis_d() +
        ylab( 'Yearly anomaly' ) + 
        xlab( 'Year' ) +
        geom_vline( xintercept = 2006,
                    lty  = 2 ) + 
        theme( axis.title = element_text( size = 10),
               axis.text  = element_text( size  = 10,
                                          angle = 70) )

p3 <- clim_raw %>% 
  subset( measure == 'ppt_raw') %>% 
  ggplot( aes(x=year, 
              y=value) ) +
  geom_line( aes(color=measure),
             lwd = 1 ) + 
  scale_color_viridis_d() +
  ylab( 'Tot. year precip. (mm)' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2 ) + 
  theme( axis.title = element_text( size = 10),
         axis.text  = element_text( size  = 10,
                                    angle = 70),
         legend.position = 'none' ) 
  
  
p4 <- clim_raw %>% 
  subset( measure == 'tmp_raw') %>% 
  ggplot( aes(x=year, 
              y=value) ) +
  geom_line( aes(color=measure),
             lwd = 1 ) + 
  scale_color_viridis_d() +
  ylab( 'Mean temperature (mm)' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2 ) + 
  theme( axis.title = element_text( size = 10),
         axis.text  = element_text( size  = 10,
                                    angle = 70),
         legend.position = 'none' )
    

g <- arrangeGrob(p1,p3,p2,p4,nrow=2,ncol=2)
ggsave(filename = 'results/climate/raw_yr_mo_anom.tiff',
       plot = g,
       dpi = 300, width = 8, height = 6, units = "in",
       compression = 'lzw')

 
# raw time series (all) ---------------------------------------


clim_raw %>% names
clim_yr %>% names

clim_raw_all <- clim_yr %>% 
                  subset( grepl('oni|spei', measure) ) %>% 
                  bind_rows( clim_raw )
  


p1 <- clim_raw_all %>% 
        subset( measure == 'tmp_raw' ) %>% 
        ggplot( aes(x = year, 
                    y = value) ) +
        geom_line( lwd = 2,
                   color = '#E69F00' ) + 
        scale_colour_colorblind() +
        ylab( 'Temperature (°C)' ) + 
        xlab( 'Year' ) +
        geom_vline( xintercept = 2006,
                    lty  = 2 ) + 
        theme( axis.title      = element_text( size = 10),
               legend.position = 'none' )  

p2 <- clim_raw_all %>% 
        subset( measure == 'ppt_raw' ) %>% 
        ggplot( aes(x = year, 
                    y = value) ) +
        geom_line( lwd = 2,
                   color = "#56B4E9" ) + 
        scale_colour_colorblind() +
        ylab( 'Precipitation (mm)' ) + 
        xlab( 'Year' ) +
        geom_vline( xintercept = 2006,
                    lty  = 2 ) + 
        theme( axis.title      = element_text( size = 10),
               legend.position = 'none' ) 


p3 <- clim_raw_all %>% 
        subset( measure == 'oni_anom' ) %>% 
        ggplot( aes(x = year, 
                    y = value) ) +
        geom_line( lwd = 2,
                   color = "#009E73" ) + 
        scale_colour_colorblind() +
        ylab( 'ONI' ) + 
        xlab( 'Year' ) +
        geom_vline( xintercept = 2006,
                    lty  = 2 ) + 
        theme( axis.title      = element_text( size = 10),
               legend.position = 'none' ) 

p4 <- clim_raw_all %>% 
        subset( measure == 'spei_anom' ) %>% 
        ggplot( aes(x = year, 
                    y = value) ) +
        geom_line( lwd = 2,
                   color = "#F0E442") + 
        scale_color_colorblind() +
        ylab( 'SPEI' ) + 
        xlab( 'Year' ) +
        geom_vline( xintercept = 2006,
                    lty  = 2 ) + 
        theme( axis.title      = element_text( size = 10),
               legend.position = 'none' )  

g <- arrangeGrob(p1,p3,p2,p4,nrow=2,ncol=2)
ggsave(filename = 'results/climate/raw_yr_clim.tiff',
       plot = g,
       dpi = 300, width = 8, height = 6, units = "in",
       compression = 'lzw')


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
       dpi = 300, width = 8, height = 6, units = "in",
       compression = 'lzw')


# Final figure with only annual anomalies -----------------------------------

# convert 
clim_yr %>% 
  mutate( measure = replace(measure, measure=='tmp_anom',
                            'Temperature') ) %>% 
  mutate( measure = replace(measure, measure=='ppt_anom',
                            'Precipitation') ) %>% 
  mutate( measure = replace(measure, measure=='spei_anom',
                            'Aridity I.') ) %>% 
  mutate( measure = replace(measure, measure=='oni_anom',
                            'ONI') ) %>% 
  ggplot( aes(x=year, 
              y=value) ) +
  geom_line( aes(color=measure),
             lwd = 2 ) + 
  scale_colour_colorblind() + 
  ylab( 'Yearly anomaly' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2 ) + 
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size  = 20,
                                    angle = 0) ) + 
  ggsave('results/climate/yr_anom_vs_year.tiff', 
         width=12, height=5, compression = 'lzw')


# Dual two axes -------------


clim_yr %>% head
spread(clim_raw, measure, value) %>% 
  .$ppt_raw %>% sd

735 + (244.5898*2)

# attempt two axes
# spread(clim_raw, measure, value) %>% 
clim_yr %>% 
  ggplot( aes(x = year) ) + 
  geom_line( aes(y=ppt_anom),
             lwd = 2,
             color = 'black') +
  scale_y_continuous(  ) + 
  ylab( "Precipitation (mm)") +
  geom_line( aes(y=tmp_anom/2),
             lwd = 2,
             color = 'red') + 
  scale_y_continuous(sec.axis = sec_axis(~.*2,
                                         name = "Temperature (°C)")
                    )

                                         breaks = c(-1, 0, 2),
                                         labels = c('11','12','13')

  scale_colour_colorblind() + 
  ylab( 'Yearly anomaly' ) + 
  xlab( 'Year' ) +
  geom_vline( xintercept = 2006,
              lty  = 2 ) + 
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size  = 20,
                                    angle = 0) )


# load library for axf data
# devtools::install_github('MarkusLoew/ReadAxfBOM')
library(ReadAxfBOM) 

# import data
obs <- ReadAxfBOM("http://www.bom.gov.au/fwo/IDV60901/IDV60901.94866.axf")

# show first observations
head(obs)


# two scales
ggplot(obs, aes(x = Timestamp)) +
  geom_line(aes(y = air_temp, colour = "Temperature")) +
  geom_line(aes(y = rel_hum/5, colour = "Humidity")) +
  scale_y_continuous(sec.axis = sec_axis(~.*5, 
                                                name = "Relative humidity [%]")
                            )
# now adding the secondary axis, following the example in the help file ?scale_y_continuous
# and, very important, reverting the above transformation


# modifying colours and theme options
p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Air temperature [°C]",
              x = "Date and time",
              colour = "Parameter")
p <- p + theme(legend.position = c(0.8, 0.9))
p