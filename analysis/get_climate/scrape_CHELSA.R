# PRELIMINARY PRELIMINARY
# Just started script. working on downloading this data
setwd("C:/cloud/Dropbox/lupine/")
library(dplyr)
library(tidyr)
library(archive)
library(raster)
library(measurements)
library(leaflet)
sA_dir      <- 'C:/cloud/Dropbox/lupine/data/weather_station/chelsa'


# convert lat/lon in decimal form
conv_plot_coord <- function(lat_in, lon_in, from_unit){
  
  coord_df <- data.frame( lat = conv_unit(lat_in,  
                                          from = from_unit, 
                                          to = 'dec_deg'),
                          lon = conv_unit(lon_in, 
                                          from = from_unit, 
                                          to = 'dec_deg'),
                          stringsAsFactors = F) %>% 
                mutate(   lat = as.numeric(lat),
                          lon = as.numeric(lon) )
  
  return(coord_df)
  
}

# from degrees+minutes+seconds to decimal degrees
leaflet(data = conv_plot_coord('38 05 39', '-122 57 00', 'deg_min_sec') ) %>% 
  addTiles() %>% 
  addCircleMarkers(~lon, ~lat)


# site coordinates 
site_coord<- conv_plot_coord('38 05 39', '-122 57 00', 'deg_min_sec') %>% 
                setNames( c('Lat', 'Long') ) %>% 
                dplyr::select(Long, Lat)


# CHELSA set up -------------------------------------------------------------------------------

# what do I need from CHELSA?
chelsa_df <- expand.grid( variable = c('tmean'),
                          year     = c(2010:2013),
                          month    = c(1:12),
                          stringsAsFactors = F) %>% 
                arrange(variable,year,month)

# set up reading
read_dir  <- 'https://www.wsl.ch/lud/chelsa/data/timeseries/tmean/'
# read_dir  <- 'https://www.wsl.ch/lud/chelsa/data/timeseries/tmin/'
# CHELSA_tmin_2_2010_V1.2.1.tif

# produce file name based on index 'ii'
produce_file_name <- function(ii){
  
  file_n    <- paste0('CHELSA_tmean_',
  # file_n    <- paste0('CHELSA_tmin_',
                       chelsa_df$month[ii],
                       '_',
                       chelsa_df$year[ii],
                       '_V1.2.1.tif')
  
  return(file_n)
  
}

# get all file links (from file name)
file_names <- lapply(1:nrow(chelsa_df), produce_file_name) %>% unlist
file_links <- paste0(read_dir, file_names)
file_dest  <- paste0('C:/temp_dir/', file_names)

# Extract CHELSA DATA ------------------------------------------------------------------------

# extract year and monthly data
extract_year_month <- function(ii){
  
  Sys.sleep(15)
  # extac with archive 
  # devtools::install_github('jimhester/archive')
  file_path <- file_links[ii]
  download.file( file_path, destfile = file_dest[ii], mode = "wb")
  # archive_extract( archive(file_dest[ii]), "temp_dir")
  # # extract with 7z directly. This does extract directly in getwd()
  # system('"C:\\Program Files\\7-Zip\\7z" x "C:\\cloud\\MEGA\\Projects\\sApropos\\analyses\\CHELSA_temp_1979_01.7z"')
  
  # get climate information ----------------------------------------------------------------
  
  # read raster
  raster_file <- grep('.tif$',list.files('C:/temp_dir'), value=T)
  repP        <- raster( paste0('C:/temp_dir/',raster_file) )
  
  # extract info 
  values_clim <- raster::extract(repP, site_coord, method = 'bilinear')
  clim_df     <- data.frame( variable = chelsa_df$variable[ii],
                             year     = chelsa_df$year[ii],
                             month    = chelsa_df$month[ii],
                             value    = values_clim )
  
  file.remove( paste0('C:/temp_dir/',list.files('C:/temp_dir/')) )
  
  print(ii)
  
  return(clim_df)

}

start <- Sys.time()
range_clim  <- 1:60
climate_all <- lapply(range_clim, extract_year_month)
Sys.time() - start

climate_df <- Reduce(function(...) bind_rows(...), climate_all ) %>% 
                mutate( value = (value /10) - 273.15 )

write.csv(climate_df, paste0(sA_dir, '/climate_chelsa_2010_2011.csv'), row.names=F)
# write.csv(climate_df, paste0(sA_dir, '/climate_chelsa_tmin_2010_2011.csv'), row.names=F)
