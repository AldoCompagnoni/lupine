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
                  

# calculate yearly anomalies as in analyses ---------------------------
years   <- 2005:2018
m_obs   <- 5
m_back  <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  mutate(x, 
         avgt0   = x %>% select(V1:V12) %>% rowSums ) %>% 
    select(year, avgt0) %>% 
    setNames( c('year', paste0(var,'_anom')) )
  
}

# format climate - need to select climate predictor first 
ppt_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs) %>% 
              year_anom('ppt')

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              year_anom('tmp')
  
enso_mat <- enso %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_anom('oni')

clim_anom <- Reduce(function(...) full_join(...), 
                    list(ppt_mat, tmp_mat, enso_mat) )

# set up "demographic year"
demo_year <- function(x_df){
  mutate(x_df,
         year = replace(year, 
                        month_num > 5, 
                        year[month_num > 5] + 1) )
}
clim   <- demo_year( clim )
enso   <- demo_year( enso ) %>% 
              subset( clim_var == 'oni' ) %>%
              select( -mon, -clim_var ) %>%
              rename( oni = clim_value )


# Format monthly data ----------------------------------------------------

# climate monthly wide (by measure) form
clim_mon_w <- spread(clim, clim_var, clim_value) %>% 
                left_join( enso ) %>% 
                subset( !is.na(oni) )

# climate monthly wide-wide (columns represent month-by-measure) form
clim_mon_w_w <- clim %>% 
                  mutate( clim_m = paste(clim_var,month_num,sep='_') ) %>% 
                  select( -month_num, -clim_var ) %>% 
                  spread(clim_m, clim_value)

# Format yearly data ----------------------------------------------------

# Aggregate 
clim_ppt   <- clim %>% 
                select( -month_num ) %>% 
                subset( clim_var == 'ppt' ) %>% 
                group_by( year ) %>% 
                summarise( ppt = sum(clim_value) )

# temperature 
clim_tmp   <- clim %>% 
                select( -month_num ) %>% 
                subset( clim_var == 'tmean' ) %>% 
                group_by( year ) %>% 
                summarise( tmp = mean(clim_value) )

# ONI 
clim_oni   <- enso %>% 
                select( -month_num ) %>% 
                group_by( year ) %>% 
                summarise( oni = mean(oni) )


# yearly climatic values
clim_yr    <- full_join( clim_ppt, clim_tmp ) %>% 
                left_join( clim_oni ) %>% 
                subset( !is.na(oni) ) %>% 
                left_join( clim_anom )

# store demographic yearly climate
write.csv( subset(clim_yr, year > 2004),
           'climate_demographic_year_lupine.csv', row.names=F)

# plots ------------------------------------------------------------
tiff('results/climate/monthly_tmp_vs_ppt.tiff', 
     unit="in", width=5, height=6.3, res=400, compression="lzw")

par( mfrow = c(2,1), mar=c(3,3,0.2,0.2), 
     mgp = c(2,1,0) )
plot(tmean ~ ppt,data=clim_mon_w)
plot(tmp ~ ppt,data=clim_yr)

dev.off()

tiff('results/climate/monthly_pairs_tmp_ppt_oni.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

par(mar=c(0,0,0.2,0.2),oma=c(0,0,0,0))
clim_mon_w %>% 
  subset( year > 2004 ) %>% 
  select( -year, -month_num ) %>% 
  pairs

dev.off()

tiff('results/climate/yearly_pairs_tmp_ppt_oni.tiff', 
     unit="in", width=6.3, height=6.3, res=400, compression="lzw")

clim_yr_study <- clim_yr %>% 
  subset( year > 2004 ) %>% 
  mutate( col = 'black') %>% 
  mutate( col = replace(col, year == 2015, 'red'))

cor( select(clim_yr_study, -year, -col) )

par(mfrow=c(2,2), mar=c(3.5,3.5,0.2,0.2),
    mgp=c(2,0.8,0), oma=c(0,0,0,0))
plot(ppt ~ tmp, data=select(clim_yr_study, -year), col = col, pch =16,
     xlab = "Annual temperature", ylab = 'Annual precipitation')
plot(tmp ~ oni, data=select(clim_yr_study, -year), col = col, pch =16,
     xlab = "Annual ONI index", ylab = 'Annual temperature')
plot(ppt ~ oni, data=select(clim_yr_study, -year), col = col, pch =16,
     xlab = "Annual ONI index", ylab = 'Annual precipitation')
# pairs(select(clim_yr, -year), col = col )

dev.off()

cor(select(clim_mon_w, -year, -month_num) )
cor( select(clim_yr, -year) )


# monthly ppt~temp info
tiff('results/climate/month_vs_month_tmp_ppt_corr.tiff', 
     unit="in", width=6.3, height=8, res=400, compression="lzw")

par( mfrow = c(4,3), mar=c(3,3,0.2,0.2), 
     mgp = c(1.7,0.6,0), cex.lab=1 )

clim_mon_w_w <- subset(clim_mon_w_w, year > 2004)

plot(ppt_1  ~ tmean_1, data=clim_mon_w_w)
abline(lm(ppt_1  ~ tmean_1, data=clim_mon_w_w))
plot(ppt_2  ~ tmean_2, data=clim_mon_w_w)
abline(lm(ppt_2  ~ tmean_2, data=clim_mon_w_w))
plot(ppt_3  ~ tmean_3, data=clim_mon_w_w)
abline(lm(ppt_3  ~ tmean_3, data=clim_mon_w_w))
plot(ppt_4  ~ tmean_4, data=clim_mon_w_w)
abline(lm(ppt_4  ~ tmean_4, data=clim_mon_w_w))

plot(ppt_5  ~ tmean_5, data=clim_mon_w_w)
abline(lm(ppt_5  ~ tmean_5, data=clim_mon_w_w))
plot(ppt_6  ~ tmean_6, data=clim_mon_w_w)
abline(lm(ppt_6  ~ tmean_6, data=clim_mon_w_w))
plot(ppt_7  ~ tmean_7, data=clim_mon_w_w)
abline(lm(ppt_7  ~ tmean_7, data=clim_mon_w_w))
plot(ppt_8  ~ tmean_8, data=clim_mon_w_w)
abline(lm(ppt_8  ~ tmean_8, data=clim_mon_w_w))

plot(ppt_9  ~ tmean_9, data=clim_mon_w_w)
abline(lm(ppt_9  ~ tmean_9, data=clim_mon_w_w))
plot(ppt_10 ~ tmean_10, data=clim_mon_w_w)
abline(lm(ppt_10  ~ tmean_10, data=clim_mon_w_w))
plot(ppt_11 ~ tmean_11, data=clim_mon_w_w)
abline(lm(ppt_11  ~ tmean_11, data=clim_mon_w_w))
plot(ppt_12 ~ tmean_12, data=clim_mon_w_w)
abline(lm(ppt_12  ~ tmean_12, data=clim_mon_w_w))

dev.off()


# tmp_cia <- clim_mon_w_w %>% .[-1,] %>% .[-32,] %>% 
#   .[,grep('tmean_', names(clim_mon_w_w), value=T )] %>% 
#   cor
# 
# ppt_cia <- clim_mon_w_w %>% .[-1,] %>% .[-32,] %>% 
#   .[,grep('ppt_', names(clim_mon_w_w), value=T )] %>% 
#   cor 
# 
# ppt_cia[upper.tri(ppt_cia)] %>% mean
# tmp_cia[upper.tri(tmp_cia)] %>% mean
