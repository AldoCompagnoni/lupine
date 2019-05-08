# Hyper-simplified IPM for validation
# I compare observed lambda to deterministic lambda
# only for sites BS, NB, and DR
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
library(bbmle)
source('analysis/vital_rates/plot_binned_prop.R')

# read in and format data
lupine_df  <- read.csv('data/lupine_all.csv')
site_yr_df <- select(lupine_df, year, location) %>% 
                unique %>% 
                mutate( Site = gsub(' \\([0-9]\\)','',location) ) %>% 
                mutate( Site = toupper(Site) ) %>% 
                subset( year > 2004 )

# consumption 
cons_df    <- read_xlsx('data/consumption.xlsx') %>%
                rename( year = Year ) %>% 
  
                mutate( Mean_consumption = Mean_consumption %>% as.numeric) %>% 
                select( year, Site, Mean_consumption) %>% 
                mutate( Site = toupper(Site) ) %>% 
                # expand potential "cases"
                right_join( site_yr_df ) %>% 
                select( -Site )

# abortion 
abor_df  <- read.csv('data/lupine_all.csv') %>% 
                subset( !is.na(flow_t0) & flow_t0 == 1 ) %>% 
                subset( !is.na(numrac_t0) ) %>% 
                # remove non-flowering individuals
                subset( !(flow_t0 %in% 0) ) %>% 
                # remove zero fertility (becase fertility should not be 0)
                subset( !(numrac_t0 %in% 0) ) %>% 
                # only years indicated by Tiffany
                subset( year %in% c(2010, 2011, 2013:2017) ) %>% 
                # calculate abortion rates
                mutate( ab_r      = numab_t0 / numrac_t0 ) %>% 
                group_by( location, year ) %>% 
                summarise( ab_r_m = mean(ab_r, na.rm=T) ) %>% 
                ungroup %>% 
                right_join( site_yr_df ) %>% 
                select( -Site )

# data frame to plot everything together
ab_cons_df <- full_join(abor_df, cons_df) %>% 
                rename( Abortion    = ab_r_m,
                        Consumption = Mean_consumption ) %>% 
                gather(measure, value, Abortion, Consumption)

# abortion
ggplot(ab_cons_df) + 
  geom_line( aes(x     = year,
                 y     = value,
                 color = location ),
             size=2,
             alpha=0.8) +
  scale_x_continuous(breaks = 0:2100) +
  scale_color_viridis_d() + 
  ylab( 'Proportion of racemes' ) +
  theme( axis.text.x = element_text(angle=70) ) + 
  facet_grid( measure ~ 1 ) +
  ggsave('results/vital_rates/abort_consumption.tiff',
         height=6.3,width=6.3,compression='lzw')
    