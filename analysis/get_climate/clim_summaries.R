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

# fetch climate ---------------------------------------------------------------------

# format the enso data
fc_month <- function(fc_df,var){
  
  # format day one
  day_one   <- as.Date( paste0("1/1/", first(fc_df$year) ), 
                        format="%d/%m/%Y" ) 
  
  # format to get means
  clim_d    <- as.Date(1:nrow(fc_df), day_one-1) %>%
                  as.character %>%
                  as.data.frame(stringsAsFactors=F) %>%
                  separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                  bind_cols(fc_df) %>% 
                  dplyr::select(-year,-day,-day1,-clim_var) %>%
                  setNames( c("year", "month", "value") ) %>% 
                  mutate( year  = as.numeric(year),
                          month = as.numeric(month) ) %>% 
                  group_by(year,month)
            
  # calculate monthly means
  if( var == 'prec'){
    clim_out  <- clim_d %>% summarise( ppt_fc  = sum(value) ) 
  }
  
  if( var == 'airt'){
    clim_out  <- clim_d %>% summarise( airt_fc = mean(value) )
  }
 
  clim_out %>% 
    ungroup %>% 
    as.data.frame %>% 
    return()

}

# make summaries by year
anom_by_yr <- function(yrz, data_df,clim_var){
  
  fc_df   <- as.data.frame(fc_df)
  
  ids     <- (which(data_df$year == yrz & data_df$month == 6 )-11):
              which(data_df$year == yrz & data_df$month == 6 )
  
  if( grepl('airt',clim_var) ) summ <- mean(data_df[ids,clim_var], na.rm=T)
  if( grepl('ppt', clim_var) ) summ <- sum(data_df[ids,clim_var], na.rm=T)

  data.frame(year  = data_df[ids,'year'][1],
             month = summ,
             stringsAsFactors = F) %>% 
   setNames( c('year',clim_var) ) %>% 
   return()

}

# monthly means
prec_df  <- subset(clim, clim_var == "prate") %>%
              mutate( value = replace(value, value < 0, 0) ) %>% 
              fc_month('prec')
airt_df  <- subset(clim, clim_var == "airt") %>% fc_month('airt')

# yearly means
fc_ppt_form <- lapply(years, anom_by_yr, prec_df,'ppt_fc') %>% bind_rows
fc_ait_form <- lapply(years, anom_by_yr, airt_df,'airt_fc') %>% bind_rows


# format Point Reyes data -------------------------------------------------------------------
point_r_f <- paste0('data/weather_station/point_reyes/', 
                    list.files('data/weather_station/point_reyes/') )
point_r_df<- lapply(point_r_f, read.csv) %>% 
                bind_rows %>% 
                mutate( ppt   = ppt * 25.4,
                        meant = (meant - 32) * (5/9) ) %>% 
                group_by(yr, mon) %>% 
                summarise( ppt_pr   = sum(ppt,na.rm=T),
                           airt_pr  = mean(meant,na.rm=T)) %>% 
                ungroup %>% 
                left_join( 
                  data.frame( month = 1:12,
                              mon   = month.name,
                              stringsAsFactors = F) ) %>%
                rename( year        = yr ) %>% 
                as.data.frame %>% 
                arrange( year, month) %>% 
                mutate( airt_pr = replace(airt_pr, is.nan(airt_pr),NA) )
                


pr_airt_form <- lapply(years[7:length(years)], anom_by_yr, point_r_df, 'airt_pr') %>% 
                    bind_rows
pr_ppt_form  <- lapply(years[7:length(years)], anom_by_yr, point_r_df,  'ppt_pr') %>% 
                    bind_rows


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
            summarise( ppt_bod  = sum(prec,na.rm=T),
                       airt_bod = mean(mean_temp, na.rm=T) ) %>% 
            arrange(year, month) %>%
            ungroup %>% 
            mutate( year  = as.numeric(year),
                    month = as.numeric(month) ) %>% 
            as.data.frame
  

bod_airt_form <- lapply(years[8:length(years)], anom_by_yr, bod, 'airt_bod') %>% bind_rows
bod_ppt_form  <- lapply(years[8:length(years)], anom_by_yr, bod, 'ppt_bod') %>% bind_rows


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
                        airt_prism = tmean )

# calculate anomalies
prism_airt_form <- lapply(years, anom_by_yr, prism_df, 'airt_prism') %>% bind_rows
prism_ppt_form  <- lapply(years, anom_by_yr, prism_df, 'ppt_prism') %>% bind_rows


# ENSO data ----------------------------------------------------------------------

oni_df   <- subset(enso, clim_var == "oni") 
soi_df   <- subset(enso, clim_var == "soi") 

# format the enso data
format_enso <- function(yrz, enso_df,var){
  
  yr_mon <- enso_df %>%
              dplyr::select(year, month_num,clim_value) %>% 
              setNames( c('year','month',var) ) 
  
  ids     <- (which(yr_mon$year == yrz & yr_mon$month == 6 )-11):
              which(yr_mon$year == yrz & yr_mon$month == 6 )
  
  data.frame(year  = yr_mon[ids,'year'][1],
                     month = mean(yr_mon[ids,var],na.rm=T),
                     stringsAsFactors = F) %>%
    setNames( c('year',var) ) %>% 
    return()
    
}

# data formatted 
oni_form <- lapply(years, format_enso, oni_df, 'oni') %>% bind_rows
soi_form <- lapply(years, format_enso, soi_df, 'soi') %>% bind_rows


# climate summaries --------------------------------------------------------------
clim_raw <- Reduce(function(...) full_join(...), 
                     list(pr_airt_form, bod_airt_form, 
                          prism_airt_form, fc_ait_form,
                          pr_ppt_form, bod_ppt_form,
                          prism_ppt_form, fc_ppt_form, 
                          oni_form, soi_form) ) %>% 
              arrange(year)
        

# scaled values
scal     <- apply(select(clim_raw,-year), 2, scale) %>% 
              as.data.frame 

trans    <- paste0(clim_raw$year,'-',clim_raw$year+1)

# scaled climate summaries
clim_scl <- bind_cols(select(clim_raw,year), 
                      apply(select(clim_raw,-year), 2, scale) %>% as.data.frame ) %>% 
              tibble::add_column( transition = trans, .before = 2)

write.csv(clim_scl, 'results/climate/year_anomalies.csv',row.names=F )
