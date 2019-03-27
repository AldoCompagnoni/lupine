setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)

f_name <- 'UCC_ghcn_USW00093245_2018_06_11_1528730025.csv'
dir    <- 'data/weather_station/UCC_ghcn_USW00093245_2018_06_11_1528730025/'

list.files('data/weather_station')

bod    <- read.csv(paste0(dir, f_name), skip = 15) %>% 
            .[-c(1:14),] %>% 
            dplyr::select(-Ref.Evapotranspiration,
                          -Snow.Depth, -Snow.Fall) %>% 
            separate(Day, into = c('day','month','year'), 
                     sep = '-')
bod[bod=='M'] = NA
         

