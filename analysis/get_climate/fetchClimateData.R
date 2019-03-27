rm(list=ls(all=TRUE)); graphics.off();
setwd("C:/cloud/MEGA/Projects/SIDE/lupine")
require(RFc)
require(dplyr)
require(magrittr)
# Precipitation rate variable in FC is prate, units mm/month
# Temperature variable is airt, units degrees C


# fetch lupine climate data -----------------------------------------------------

# fect yearly climate data 
fetch_yr_clim <- function(yr,var_name){
  
  fcTimeSeriesDaily(variable  = var_name,
                    latitude  = 38.085, 
                    longitude = -122.967,
                    firstYear = yr, 
                    lastYear  = yr)
  
}

#  format fetched climate
form_fetched <- function(raw_df, var_name){
  
  yrs      <- 1967:2017
  yr_by_yr <- function(x, y, var_name){
    
    data.frame( year     = x,
                day      = 1:length(y$value),
                clim_var = var_name,
                value    = as.numeric(y$value) )
    
  }
  
  form_df_l <- Map(yr_by_yr, yrs, raw_df, var_name)
  
  Reduce(function(...) rbind(...), form_df_l)
  
}

# wrapper for two functions above
fetch_var_clim <- function(var_name){
  
  # fetch yearly climate climate
  raw_clim <- lapply(1967:2017, fetch_yr_clim, var_name)
  
  print(paste0("done downloading ", var_name))
  
  # format RFc in a data frame
  form_fetched(raw_clim, var_name)
  
}

clim_vars <- c("prate","airt", "soilmoist", "frs", "wvp", "wvsp", 
               "windspeed", "sunp", "relhum_land", "abshum")

clim_var_l1  <- lapply( clim_vars[1:2], fetch_var_clim)
clim_var_l3  <- lapply( clim_vars[3], fetch_var_clim)
#clim_var_l4  <- lapply( clim_vars[4], fetch_var_clim)
clim_var_l5  <- lapply( clim_vars[5], fetch_var_clim)
clim_var_l67 <- lapply( clim_vars[6:7], fetch_var_clim)

clim_var_l8 <- lapply( clim_vars[8], fetch_var_clim)

clim_var_l910 <- lapply( clim_vars[9:10], fetch_var_clim)

clim_var_l   <- list(clim_var_l1[[1]],clim_var_l1[[2]],
                     clim_var_l3[[1]],clim_var_l5[[1]])


clim_var_l   <- list(clim_var_l5[[1]],clim_var_l67[[1]],clim_var_l67[[2]],
                     clim_var_l8[[1]],clim_var_l910[[1]],clim_var_l910[[2]])



clim_var_df <- Reduce(function(...) rbind(...), clim_var_l) 

write.csv(clim_var_df, "data/lupine_fc_vars2.csv", row.names=F)
write.csv(clim_var_df, "data/lupine_fc_vars.csv", row.names=F)
cia <- read.csv("data/lupine_fc_vars.csv")
