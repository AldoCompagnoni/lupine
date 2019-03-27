# Produce a dataframe of climate predictors
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

# format prism data -------------------------------------------------------------------------
clim      <- read.csv("data/lupine_fc_vars.csv")
years     <- 2002:2018
m_obs     <- 6
m_back    <- 12

prec_df  <- subset(clim, clim_var == "prate") %>%
                subset( year > 2001 ) %>% 
                mutate( value = replace(value, value < 0, 0) ) 
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

fc_df     <- full_join(prec_fc, airt_fc) #%>% mutate( source = 'fc' )

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
                        meant_prism = tmean ) #%>% 
                #mutate(source = 'prism' )

# CHELSA -----------------------------------------------------------------------
chelsa_df <- read.csv('data/weather_station/chelsa/climate_chelsa_2010_2011.csv') %>% 
                rename( meant_chelsa = value,
                        month_num    = month ) %>% 
                select(  year, month_num, meant_chelsa )

# format Point Reyes data -------------------------------------------------------------------
point_r_f <- paste0('data/weather_station/point_reyes/', 
                    list.files('data/weather_station/point_reyes/') )

pr_df     <- lapply(point_r_f[-c(1:3)], read.csv) %>% 
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
                rename( year = yr ) %>% 
                # mutate( source = 'pr' ) %>% 
                dplyr::select(- mon )

# format BODEGA data ----------------------------------------------------------------
dir     <- 'C:/cloud/Dropbox/lupine/data/weather_station/UCC_ghcn_USW00093245_2018_06_11_1528730025/'
fil_n   <- 'UCC_ghcn_USW00093245_2018_06_11_1528730025.csv'
bod_df  <- read.csv( paste0(dir, fil_n), skip=15 )[-c(1:14),] %>% 
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
                       meant_bod = mean(mean_temp, na.rm=T) ) %>% 
            arrange(year, month) %>%
            ungroup %>% 
            mutate( year  = as.numeric(year),
                    month = as.numeric(month) ) %>% 
            rename( month_num = month ) #%>% 
            #mutate( source = 'bod' )
  
# FT ROSS --------------------------------------------------------------------------------
dir     <- 'C:/cloud/Dropbox/lupine/data/weather_station/UCC_ghcn_USC00043191_2018_09_01_1535816841/'
fil_n   <- 'UCC_ghcn_USC00043191_2018_09_01_1535816841.csv'
ftr_df  <- read.csv( paste0(dir, fil_n), skip=15 )[-c(1:14),] %>% 
            dplyr::select(-Snow.Depth,-Snow.Fall) %>% 
            setNames( c('Day','prec','min_temp','max_temp') ) %>% 
            mutate( prec      = replace(prec,     prec     == "M" | prec == "S" | prec == "T", NA),
                    min_temp  = replace(min_temp, min_temp == "M" | min_temp == "S", NA),
                    max_temp  = replace(max_temp, max_temp == "M" | max_temp == "S", NA) ) %>%
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
            summarise( ppt_ftr  = sum(prec,na.rm=T),
                       meant_ftr = mean(mean_temp, na.rm=T) ) %>% 
            arrange(year, month) %>%
            ungroup %>% 
            mutate( year  = as.numeric(year),
                    month = as.numeric(month) ) %>% 
            rename( month_num = month ) #%>% 
            #mutate( source = 'bod' )
  
names(ftr_df)
names(bod_df)
names(pr_df)
names(prism_df)
names(fc_df)


# all climate variables put together
clim_df <- Reduce(function(...) inner_join(...), 
                  list(pr_df, bod_df, prism_df, 
                       fc_df, ftr_df) ) %>% #chelsa_df, 
              arrange(year,month_num) %>% 
              mutate( time = paste(month_num,year,sep='/') %>% as.character ) %>% 
              mutate( time = paste0('15/',time)) %>% 
              mutate( time_date = as.Date(time, '%d/%m/%Y'))
dates_x   <- subset(clim_df, month_num == 1)

# Bodega/Point Reyes/Prism
tiff('results/climate/ppt_compare_pr_bod_prism.tiff',
     unit="in", width=9, height=6.3, res=600,compression="lzw")

par(mar=c(3.5,4,0.2,0.2),mgp=c(2,0.8,0))
matplot(clim_df$time_date,select(clim_df,ppt_bod,ppt_pr, ppt_prism),
        type='l', lty=1, xaxt="n",
        ylab = 'monthly precipitation (mm)',
        xlab = 'Year')
axis(1,at=dates_x$time_date,labels=dates_x$year)
legend('topleft', c('Bodega','Point Reyes','Prism'), 
       col=1:3, lty=1, bty='n')

dev.off()

# Bodega/Point Reyes/Prism BEFORE 2017
tiff('results/climate/ppt_compare_pr_bod_prism_pre_2017.tiff',
     unit="in", width=9, height=6.3, res=600,compression="lzw")

clim_pre_17 <- subset(clim_df, year < 2017)
date_pre_17 <- subset(clim_pre_17, month_num == 1)

par(mar=c(3.5,4,0.2,0.2),mgp=c(2,0.8,0))
matplot(clim_pre_17$time_date,
        select(clim_pre_17,ppt_bod,ppt_pr, ppt_prism),
        type='l', lty=1, xaxt="n",
        ylab = 'monthly precipitation (mm)',
        xlab = 'Year')
axis(1,at=date_pre_17$time_date,labels=date_pre_17$year)
legend('topleft', c('Bodega','Point Reyes','Prism'), 
       col=1:3, lty=1, bty='n')

dev.off()

# Bodega vs. prism
tiff('results/climate/ppt_compare_bod_prism.tiff',
     unit="in", width=9, height=6.3, res=600,compression="lzw")

par(mar=c(3.5,4,0.2,0.2),mgp=c(2,0.8,0))
matplot(clim_df$time_date,select(clim_df,ppt_bod, ppt_prism),
        type='l', lty=1:2, xaxt="n", col=c(1,3),
        ylab = "monthly precipitation (mm)",
        xlab = 'Year')
axis(1,at=dates_x$time_date,labels=dates_x$year)
legend('topleft', c('Bodega','Prism'), 
       col=c(1,3), lty=c(1,2), bty='n')

dev.off()


# temperature -------------------------------------------

# Bodega/Point Reyes/Prism/Fetch Climate
tiff('results/climate/airt_compare_pr_bod_prism_ftr_fc.tiff',
     unit="in", width=9, height=6.3, res=600,compression="lzw")

par(mar=c(3.5,4,0.2,0.2),mgp=c(2,0.8,0))
matplot(clim_df$time_date,
        select(clim_df,
               meant_bod,
               meant_pr, 
               meant_prism,
               meant_ftr,
               meant_fc), # meant_chelsa
        type='l', lty=1, xaxt="n",
        ylim = c(8,23),
        ylab = 'Mean Monthly Temperature (Celsius)',
        xlab = 'Year')
axis(1,at=dates_x$time_date,labels=dates_x$year)
legend('topleft', c('Bodega','Point Reyes','Prism','ftr','Fetch climate'), 
       col=1:5, lty=1, bty='n')

dev.off()

# Bodega/Point Reyes/Fetch Climate
tiff('results/climate/airt_compare_pr_bod_prism.tiff',
     unit="in", width=9, height=6.3, res=600,compression="lzw")

par(mar=c(3.5,4,0.2,0.2),mgp=c(2,0.8,0))
matplot(clim_df$time_date,
        select(clim_df,
               meant_bod,
               meant_pr, 
               meant_prism),
        type='l', lty=1, xaxt="n",
        ylim = c(8,23),
        col=c(1:3), 
        ylab = 'Mean Monthly Temperature (Celsius)',
        xlab = 'Year')
axis(1,at=dates_x$time_date,labels=dates_x$year)
legend('topleft', c('Bodega','Point Reyes','Prism'), 
       col=c(1:3), lty=1, bty='n')

dev.off()
