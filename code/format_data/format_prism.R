# Produce a dataframe of climate predictors
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
options(stringsAsFactors = F)

# format prism data -------------------------------------------------------------------------
prism_f   <- paste0('data/weather_station/prism/', list.files('data/weather_station/prism') )
prism_l   <- lapply(prism_f, read.csv) %>% bind_rows
yr_prism  <- prism_l$year %>% unique %>% .[-c(1:4,32)]
prism_df  <- prism_l %>% 
                bind_rows %>% 
                rename( month_num  = month,
                        clim_value = value, 
                        clim_var   = variable ) %>% 
                select( year,month_num, clim_var, clim_value)  

# write out formatted file
write.csv(prism_df,
          'C:/cloud/Dropbox/lupine/data/prism_point_reyes_87_18.csv',
          row.names=F)
