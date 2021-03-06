# semi-automatic download of PRISM climate data
setwd("C:/")
library(dplyr)
library(tidyr)
devtools::install_github('jimhester/archive')
library(archive)
library(raster)
library(prism)
library(RCurl)

# Point Reynes location
site_coord <- c(-122.959824,38.109108) %>% 
                # in matrix form to appeal raster::extract
                matrix(nrow=1, ncol=2)

# PRISM set up -------------------------------------------------------------------------------

# what do we need from PRISM?
prism_df <- expand.grid( variable = c('ppt','tmean'),
                          year     = c(1987:2018),
                          month    = c(paste0('0',1:9),paste0(10:12)),
                          stringsAsFactors = F) %>% 
                arrange(variable,year,month)

# set up reading path
read_dir  <- 'ftp://prism.oregonstate.edu/monthly/'

# produce file name based on index 'ii'
produce_file_name <- function(ii){
  
  if( prism_df$variable[ii] == 'ppt'){
    file_root  <- paste0(prism_df$variable[ii],'/',prism_df$year[ii],
                         '/PRISM_',prism_df$variable[ii],'_stable_4kmM3_',
                          prism_df$year[ii],prism_df$month[ii],'_bil.zip')
  }
  
  if( prism_df$variable[ii] == 'tmean'){
    file_root  <- paste0(prism_df$variable[ii],'/',prism_df$year[ii],
                         '/PRISM_',prism_df$variable[ii],'_stable_4kmM2_',
                          prism_df$year[ii],prism_df$month[ii],'_bil.zip')
  }
  
  return(file_root)
  
}

# get all file names 
file_names <- lapply(1:nrow(prism_df), produce_file_name) %>% unlist

# get all the file links (from file name)
file_links <- paste0(read_dir,file_names)

# produce file destinations (put it all under C:/)
file_dest  <- gsub("tmean/[0-9]{4}/|ppt/[0-9]{4}/","",file_names) %>% 
                paste0('C:/',.)

file_links_n <- file_links

# update links "by hand"
# the most recent files are "provisional", not "stable".
file_links_n[c(375:378)] <- gsub('_stable_','_provisional_',file_links[c(375:378)])
file_links_n[c(759:762)] <- gsub('_stable_','_provisional_',file_links[c(759:762)])
file_links_good <- file_links_n


# Extract PRISM DATA ------------------------------------------------------------------------

# extract year and monthly data
extract_year_month <- function(ii){
  
  # extac with archive 
  # devtools::install_github('jimhester/archive')
  file_path <- file_links_good[ii]
  
  download.file( file_path, destfile = file_dest[ii], mode = "wb")
  archive_extract( archive(file_dest[ii]), "temp_dir")
  # # extract with 7z directly. This does extract directly in getwd()
  # system('"C:\\Program Files\\7-Zip\\7z" x "C:\\cloud\\MEGA\\Projects\\sApropos\\analyses\\CHELSA_temp_1979_01.7z"')
  
  # get climate information ----------------------------------------------------------------
  
  # read raster
  raster_file <- grep('.bil$',list.files('temp_dir'), value=T)
  rast_stack  <- raster(paste0('temp_dir/',raster_file) )
  
  # extract climatic info 
  values_clim <- raster::extract(rast_stack, site_coord,layer=1) #, method = 'bilinear')
  clim_df     <- data.frame( variable = prism_df$variable[ii],
                             year     = prism_df$year[ii],
                             month    = prism_df$month[ii],
                             value    = values_clim,
                             stringsAsFactors = F)
  
  file.remove( paste0('temp_dir/',list.files('temp_dir/')) )
  file.remove( file_dest[ii] )
  
  print(ii)
  
  return(clim_df)

}

# # extract year and monthly data
# # Service function to check that links work/exist
# check_links <- function(ii){
# 
#   url.exists(file_links_good[ii])
# 
# }
# url.exists(file_links_good[8])
# check that all file links actually exist
# check_links <- sapply(1:length(file_links_good),check_links)

start <- Sys.time()
climate_all <- lapply(1:length(file_links_good), 
                      extract_year_month)
Sys.time() - start

# put it all in one data frame
climate_df <- climate_all %>% bind_rows()

# store the data frame!
# write.csv(climate_df, 
#          'C:/cloud/Dropbox/lupine/data/weather_station/756_762_prism_point_reyes.csv', 
#          row.names=F)
