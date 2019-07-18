# Hyper-simplified IPM for validation
# I compare observed lambda to deterministic lambda
# only for sites BS, NB, and DR
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
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
                summarise( ab_r_m  = mean(ab_r, na.rm=T),
                           rep_ind = n() ) %>% 
                ungroup %>% 
                right_join( site_yr_df ) %>% 
                select( -Site )

# data frame to plot everything together
ab_cons_df <- full_join(abor_df, cons_df) %>% 
                rename( Abortion    = ab_r_m,
                        Consumption = Mean_consumption ) %>% 
                gather(measure, value, Abortion, Consumption)


# plots ------------------------------------------------------

# abortion/consumption graph
ggplot(ab_cons_df) + 
  geom_line( aes(x     = year,
                 y     = value,
                 color = location ),
             size=2,
             alpha=0.8) +
  scale_x_continuous(breaks = 0:2100) +
  scale_colour_colorblind() + 
  ylab( 'Proportion of racemes' ) +
  theme( axis.text.x = element_text(angle=70) ) + 
  facet_grid( measure ~ 1 ) +
  ggsave('results/vital_rates/abort_consumption.tiff',
         height=6.3,width=6.3,compression='lzw')
    

# abortion 
abor_df %>% 
  subset( year > 2009 & year < 2018 ) %>% 
  rename( pop = location ) %>% 
  mutate( pop = gsub('\\(|\\)|[0-8]| ','',pop)) %>% 
  mutate( pop = gsub('99','9',pop) ) %>% 
  ggplot() +
  geom_line( aes(x     = year, 
                 y     = ab_r_m, 
                 color = pop ),
             lwd = 1.5) +
  scale_colour_colorblind() +
  theme( axis.title = element_text( size = 20), 
         axis.text = element_text( size = 15), 
         legend.text = element_text( size = 10), 
         legend.title = element_text( size = 15), 
         legend.margin =  margin(t = 0, r = 10, b = 0, l = 2),
         legend.box.margin = margin(-10,-10,-10,-10) ) +
  labs( x     = 'Year',
        y     = 'Average abortion rate',
        color = 'Population' ) +
  ggsave('results/vital_rates/abortion.tiff',
         height=5,width=6.3,compression='lzw')

# consumption 
cons_df %>% 
  # subset( year > 2009 & year < 2018 ) %>% 
  rename( pop = location ) %>% 
  mutate( pop = gsub('\\(|\\)|[0-8]| ','',pop)) %>% 
  mutate( pop = gsub('99','9',pop) ) %>% 
  ggplot() +
  geom_line( aes(x     = year, 
                 y     = Mean_consumption, 
                 color = pop ),
             lwd = 1.5) +
  scale_colour_colorblind() +
  theme( axis.title = element_text( size = 20), 
         axis.text = element_text( size = 15), 
         legend.text = element_text( size = 10), 
         legend.title = element_text( size = 15), 
         legend.margin =  margin(t = 0, r = 10, b = 0, l = 2),
         legend.box.margin = margin(-10,-10,-10,-10) ) +
  labs( x     = 'Year',
        y     = 'Average consumption rate',
        color = 'Population' ) +
  ggsave('results/vital_rates/consumption.tiff',
         height=5,width=6.3,compression='lzw')


# abortion 
abor_tot_df  <- read.csv('data/lupine_all.csv') %>% 
                  subset( !is.na(flow_t0) & flow_t0 == 1 ) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  # only years indicated by Tiffany
                  subset( year %in% c(2010, 2011, 2013:2017) ) %>% 
                  
                  group_by( location, year ) %>% 
                  # calculate abortion rates
                  summarise( tot_ab_t0  = sum(numab_t0,  na.rm=T),
                             tot_rac_t0 = sum(numrac_t0, na.rm=T),
                             rep_tot    = n() ) %>% 
                  ungroup %>% 
                  mutate( ab_r_tot = tot_ab_t0 / tot_rac_t0 )
                

abor_ctrl <- full_join( abor_df, abor_tot_df ) %>% 
                select(-tot_ab_t0, -tot_rac_t0, -rep_tot, -rep_ind ) %>% 
                rename( ab_r_ind = ab_r_m ) %>% 
                gather( measure, value, ab_r_ind:ab_r_tot ) %>% 
                subset( year > 2009 & year < 2018 )


# compare two methods to calculate abortion rates
ggplot(abor_ctrl) +
  geom_line( aes(x     = year, 
                 y     = value, 
                 lty   = measure,
                 color = location),
             lwd = 1.5,
             alpha = 0.8) +
  scale_colour_colorblind() +
  ggsave( 'results/vital_rates/abort_indiv_vs_tot.tiff', 
          height=3,width=6.3,compression='lzw' )


# store abortion rates, with sample sizes, calculated with "total mean"
select(abor_tot_df, year, location, ab_r_tot, rep_tot) %>% 
  mutate( ab_r_tot = round(ab_r_tot, 3) ) %>% 
  mutate( ab_lab   = paste0(ab_r_tot,' (',rep_tot,')') ) %>% 
  select( year, location, ab_lab ) %>% 
  spread( year, ab_lab ) %>% 
  write.csv( 'results/vital_rates/abor_rates_tot_rep.csv',
             row.names=F)

# store abortion rates, with sample sizes, calculated with "individual mean"
select(abor_df, year, location, ab_r_m, rep_ind) %>%
  subset( year %in% c(2010,2011,2013,2014,2015,2016,2017) ) %>% 
  mutate( ab_r_m = round(ab_r_m, 3) ) %>% 
  mutate( ab_lab = paste0(ab_r_m,' (',rep_ind,')') ) %>% 
  select( year, location, ab_lab ) %>% 
  spread( year, ab_lab ) %>% 
  write.csv( 'results/vital_rates/abor_rates_ind_rep.csv',
             row.names=F)

# store consumption rates, with sample sizes
read_xlsx('data/consumption.xlsx') %>% 
  select( Site, Year, nplants, Mean_consumption) %>% 
  mutate( Mean_consumption = as.numeric(Mean_consumption) ) %>% 
  mutate( cons = round(Mean_consumption, 3) ) %>% 
  mutate( cons = paste0(cons,' (',nplants,')') ) %>% 
  select( Year, Site, cons ) %>% 
  mutate( cons = replace(cons, cons == 'NA (NA)', NA) ) %>% 
  spread( Year, cons ) %>% 
  write.csv( 'results/vital_rates/cons_rates_rep.csv',
             row.names=F)
