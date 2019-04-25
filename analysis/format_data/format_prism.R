# Produce a dataframe of climate predictors
# 1. format prism data 
# 2. calculate SPEI (IMPORTANT: pay attention to missing months,
#    at end of climate time series)
library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

#1. format prism data -------------------------------------------------------------------------
prism_f   <- paste0('data/weather_station/prism/', list.files('data/weather_station/prism') )
prism_l   <- lapply(prism_f, read.csv) %>% bind_rows
yr_prism  <- prism_l$year %>% unique %>% .[-c(1:4,32)]
prism_df  <- prism_l %>% 
                bind_rows %>% 
                rename( month_num  = month,
                        clim_value = value, 
                        clim_var   = variable ) %>% 
                select( year,month_num, clim_var, clim_value)  

# add 12-month SPEI to data ----------------------------------

# create placeholder for calculations
pr_clim     <- prism_df %>% 
                 rename( MONTH = month_num ) %>% 
                 spread( clim_var, clim_value)

# Compute potential evapotranspiration (PET) and climatic water balance (BAL)
pr_clim$PET <- thornthwaite(pr_clim$tmean, 38.06667)
pr_clim$BAL <- pr_clim$ppt - pr_clim$PET

# transform in 
pr_ts <- ts(pr_clim[,-c(1,2)], 
            end=c(2018,6), 
            frequency=12)

# calculate SPEI
pr_spei12 <- spei(pr_ts[,'BAL'], 12)

# format only months with complete data (1987-2017)
spei_df   <- matrix(pr_spei12$fitted[1:(12*31)],
               nrow=31,ncol=12,
               byrow=T) %>% 
              # IMPORTANT: pay attention to missing months,
              # at end of climate time series
              rbind( c(pr_spei12$fitted[((12*31)+1):((12*31)+6)],
                       rep(NA,6) ) ) %>% 
              as.data.frame %>% 
              # add column of years (link to prism_df)
              tibble::add_column( year = prism_df$year %>% unique %>% sort,
                                  .before=1 ) %>% 
              gather( month_num, clim_value, V1:V12 ) %>% 
              mutate( month_num = gsub('V','',month_num) ) %>% 
              mutate( month_num = as.numeric(month_num) ) %>% 
              tibble::add_column(clim_var = 'spei',
                                 .before  = 3) %>% 
              arrange( year, month_num) %>% 
              # remove NAs
              subset( !is.na(clim_value ) )
        
# write out formatted file
write.csv(bind_rows(prism_df, spei_df),
          'data/prism_point_reyes_87_18.csv',
          row.names=F)
