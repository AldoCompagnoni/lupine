setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
clim      <- read.csv("data/lupine_fc_vars.csv")
enso      <- read.csv("data/enso_data.csv")
years     <- 2002:2017
m_obs     <- 6
m_back    <- 12

# ENSO data ----------------------------------------------------------------------

oni_df   <- subset(enso, clim_var == "oni") %>% 
              dplyr::select(-clim_var,-mon) %>% 
              rename( oni = clim_value )
soi_df   <- subset(enso, clim_var == "soi") 


# fetch climate ------------------------------------------------------------------

prec_df  <- subset(clim, clim_var == "prate") %>%
                subset( year > 2001 ) %>%
                mutate( value = replace(value, value < 0, 0) ) %>% 
                subset( year > 2001 )
airt_df  <- subset(clim, clim_var == "airt") %>% 
                subset( year > 2001 )

# format day one
day_one   <- as.Date( paste0("1/1/", 2002 ), 
                        format="%d/%m/%Y" ) 
  
# climate data
prec_fc    <- as.Date(1:nrow(prec_df), day_one-1) %>%
                as.character %>%
                as.data.frame(stringsAsFactors=F) %>%
                separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                bind_cols(prec_df) %>%
                dplyr::select(-year,-day) %>%
                setNames( c("year", "month", "day", "clim_var", "value") ) %>% 
                group_by(year,month) %>% 
                summarise( ppt_fc = sum(value) ) %>% 
                ungroup %>% 
                rename( month_num = month ) %>% 
                mutate( month_num = as.numeric(month_num),
                        year      = as.numeric(year) )

airt_fc    <- as.Date(1:nrow(airt_df), day_one-1) %>%
                as.character %>%
                as.data.frame(stringsAsFactors=F) %>%
                separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                bind_cols(airt_df) %>%
                dplyr::select(-year,-day) %>%
                setNames( c("year", "month", "day", "clim_var", "value") ) %>% 
                group_by(year,month) %>% 
                summarise( meant_fc = mean(value) ) %>% 
                ungroup %>% 
                rename( month_num = month ) %>% 
                mutate( month_num = as.numeric(month_num),
                        year      = as.numeric(year) )

# format prism data -------------------------------------------------------------------------
prism_f   <- paste0('data/weather_station/prism/', list.files('data/weather_station/prism') )
prism_l   <- lapply(prism_f, read.csv) %>% bind_rows
yr_prism  <- prism_l$year %>% unique %>% .[-c(1:4,32)]
prism_df  <- prism_l %>% 
                bind_rows %>% 
                rename( month_num  = month,
                        clim_value = value, 
                        clim_var   = variable ) %>% 
                spread(clim_var, clim_value) %>% 
                rename( ppt_prism   = ppt,
                        tmean_prism = tmean )

# format Point Reyes data -------------------------------------------------------------------
point_r_f <- paste0('data/weather_station/point_reyes/', 
                    list.files('data/weather_station/point_reyes/') )

point_r_df<- lapply(point_r_f[-c(1:3)], read.csv) %>% 
                bind_rows %>% 
                mutate( ppt   = ppt * 25.4,
                        meant = (meant - 32) * (5/9) ) %>% 
                group_by(yr,mon) %>% 
                summarise( ppt_pr   = sum(ppt,na.rm=T),
                           meant_pr = mean(meant,na.rm=T)) %>% 
                ungroup %>% 
                left_join( 
                  data.frame( month_num = 1:12,
                              mon       = month.name) ) %>%
                rename( year = yr )


# format BODEGA data ----------------------------------------------------------------
dir   <- 'C:/cloud/Dropbox/lupine/data/weather_station/UCC_ghcn_USW00093245_2018_06_11_1528730025/'
fil_n <- 'UCC_ghcn_USW00093245_2018_06_11_1528730025.csv'
bod   <- read.csv( paste0(dir, fil_n), skip=15 )[-c(1:14),] %>% 
            dplyr::select(-Snow.Depth,-Snow.Fall,-Ref.Evapotranspiration) %>% 
            setNames( c('Day','prec','min_temp','max_temp') ) %>% 
            mutate( prec      = replace(prec,     prec     == "M", NA),
                    min_temp  = replace(min_temp, min_temp == "M", NA),
                    max_temp  = replace(max_temp, max_temp == "M", NA) ) %>%
            mutate( prec      = as.numeric(prec),
                    min_temp  = as.numeric(min_temp),
                    max_temp  = as.numeric(max_temp) ) %>% 
            mutate( mean_temp = (min_temp + max_temp)/2 ) %>% 
            dplyr::select(Day, prec, mean_temp) %>% 
            separate( Day, into = c('year','month','day'), sep="-" ) %>%
            # exclude june in 2006 (not complete)
            # subset( !(year == '2008' & month == '06') ) %>% 
            # get monthly means
            group_by(year, month) %>% 
            summarise( prec_bod = sum(prec,na.rm=T),
                       airt_bod = mean(mean_temp, na.rm=T) ) %>% 
            arrange(year, month) %>%
            ungroup %>% 
            mutate( year  = as.numeric(year),
                    month = as.numeric(month) ) %>% 
            rename( month_num = month )
  
# pair with point reyes
mont_df <- data.frame(month_num = 1:12,
                      mon = month.name, 
                      stringsAsFactors = F)


uno = bod %>% 
  rename( prec_bod = prec,
          airt_bod = airt )
prova = point_r_df %>% 
  rename(year = yr) %>% 
  full_join( mont_df) %>% 
  full_join( uno)

prova %>% 
  ungroup %>% 
  select(prec_bod, ppt_pr, meant_pr, airt_bod) %>% 
  pairs

# years bodega
yr_bod  <- bod %>% .$year %>% unique %>% .[-c(1,11)]

# bodega formatted files
bod_airt<- bod %>% 
              select(-prec_bod ) %>% 
              rename( clim_value = airt_bod) %>%
              as.data.frame %>% 
              month_clim_form("soi", yr_bod, m_back, m_obs)

bod_prec<- bod %>% 
              select(-airt_bod ) %>% 
              rename( clim_value = prec_bod ) %>% 
              as.data.frame %>% 
              month_clim_form("soi", yr_bod, m_back, m_obs)

# bod_df  <- bind_cols(2,agg_df) %>% 
#               group_by(yr) %>% 
#               summarise( prec_bod = sum(prec),
#                          airt_bod = mean(airt) ) %>% 
#               ungroup %>% 
#               rename( year = yr) %>% 
#               right_join(clim_df) %>% 
#               mutate( prec_bod = scale(prec_bod),
#                       airt_bod = scale(airt_bod) )


# Plot Jan-to-Dec anomalies -----------------------------------------------------

# plot  
jan_dec_df <- Reduce(function(...) full_join(...), 
                     list(bod, prism_df, point_r_df, airt_fc, prec_fc, oni_df) )

# interactive graphs
tiff(paste0('results/climate/prec_anomaly.tiff'),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

jan_dec_df %>%
  subset( year < 2017) %>% 
  mutate( prec_bod    = scale(prec_bod),
          airt_bod    = scale(airt_bod),
          ppt_prism   = scale(ppt_prism),
          tmean_prism = scale(tmean_prism),
          ppt_pr      = scale(ppt_pr),
          meant_pr    = scale(meant_pr),
          ppt_fc      = scale(ppt_fc),
          meant_fc    = scale(meant_fc) ) %>% 
  dplyr::select(prec_bod, ppt_pr, ppt_prism, ppt_fc) %>%
  pairs

dev.off()


tiff(paste0('results/climate/airt_anomaly.tiff'),
     unit="in", width=6.3, height=6.3, res=600,compression="lzw")

jan_dec_df %>%
  subset( year < 2017) %>% 
  mutate( prec_bod    = scale(prec_bod),
          airt_bod    = scale(airt_bod),
          ppt_prism   = scale(ppt_prism),
          tmean_prism = scale(tmean_prism),
          ppt_pr      = scale(ppt_pr),
          meant_pr    = scale(meant_pr),
          ppt_fc      = scale(ppt_fc),
          meant_fc    = scale(meant_fc) ) %>% 
  dplyr::select(airt_bod, tmean_prism, meant_pr, meant_fc, oni) %>%
  pairs

dev.off()

