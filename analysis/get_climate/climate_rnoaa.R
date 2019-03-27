setwd('C:/cloud/Dropbox/lupine/')
library(leaflet)
library(rnoaa)
library(dplyr)
library(tidyr)
library(measurements)

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
conv_plot_coord('38 04', '-122 53', 'deg_dec_min')



# Oceanic Nino index (ONI) --------------------------------------------------------

# oni data
oni      <- read.table('http://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt',
                        header = T, stringsAsFactors = F) 

# pair "SEAS" variable with the actual month they refer to
oni_df    <- data.frame( SEAS  = oni$SEAS[1:12],
                         mon   = month.abb,
                         mon_n = 1:12,
                         stringsAsFactors = F ) %>% 
                # join with original oni data
                left_join( oni ) %>% 
                # remove absolute value (retain only anomalies)
                select( -TOTAL, -SEAS ) %>%
                select( YR, mon, mon_n, ANOM ) %>%
                arrange( YR, mon_n ) %>% 
                setNames( c('year','mon', 'month_num', 'clim_value') ) %>% 
                mutate( clim_var = 'oni' )
                # spread it in a tabular format
                # spread(mon, ANOM) 

# update file in data folder
write.csv(oni_df, 'data/oni_data.csv', row.names=F)


# Southern Oscillation Index --------------------------------------------
soi_df <- read.csv('https://www.ncdc.noaa.gov/teleconnections/enso/indicators/soi/data.csv',
                header=F, skip = 2) %>% 
              setNames( c('date', 'value') ) %>% 
              mutate( year  = substr(date,1,4),
                      mon_n = substr(date,5,6) %>% as.numeric ) %>% 
              select( -date ) %>% 
              left_join( data.frame( mon_n = 1:12,
                                     mon   = month.abb,
                                     stringsAsFactors=F) ) %>% 
              select( year, mon, mon_n, value ) %>% 
              arrange( year, mon_n ) %>% 
              setNames( c('year','mon', 'month_num', 'clim_value') ) %>% 
              mutate( clim_var = 'soi',
                      year = as.numeric(year) )

# update file in data folder
write.csv(soi_df, 'data/soi_data.csv', row.names=F)

# all indexes together
write.csv(bind_rows(oni_df, soi_df),'data/enso_data.csv',row.names=F)

# Trans-Niño Index (TNI)?


# ISD stations ----------------------------------------------------------

# Point Reyes
df <- isd_stations_search(lat = 38.06667, lon = -122.8833, radius = 20)

# plot ISD stations 
leaflet(data = df) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~longitude, ~latitude)

# get data
cia <- isd(df$usaf[1],df$wban[1],2014)
cia <- isd(df$usaf[2],df$wban[2],2011)

# NCDC stations ----------------------------------------------------------

# follow link for token
res <- ncdc_stations(extent = c(35.06667,-130.8833,
                                40.06667,-118.8833),
                     token = "ROPsBcXPMsXQKnXTkBWmTBlToXdRvsqz")

# Plot 
leaflet(data = res$data) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~longitude, ~latitude)

ncdc_datasets(datasetid = "GHCND:SZ000009480",
              token = "ROPsBcXPMsXQKnXTkBWmTBlToXdRvsqz")


# COOP stations ----------------------------------------------------

coops_search(begin_date = 19950202, end_date = 19950220, 
             station_name = 047028,
            'air_temperature', 
            datum = NULL, 
            units = "metric", 
            time_zone = "gmt",
            application = "rnoaa")








# leipzig -------------------------------------------------------------
df <- isd_stations_search(lat = 51.5, lon = 12.5, radius = 30)

leaflet(data = df) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~longitude, ~latitude)

leipzig_df <- isd(df$usaf[3],df$wban[3],1999)

# Sam's code



cat('Downloading data for',
    AllSpecies$plant.species[i],
    'from NOAA API.\n There are',
    length(StartYear:EndYear),
    'years of data for this species.',
    'Adjust your expectations accordingly\n')

# Identify the country in the lookup table and extract it's FIPS ID
CountryID <- CountryCodes$id[which(CountryCodes$name == 
                                     AllSpecies$Country[i])]

# Create a geographical bounding box to search for stations in
LatMax <- 51.33333333 + 2
LatMin <- 51.33333333 - 2
LonMax <- 12.38333333 + 2
LonMin <- 12.38333333 - 2

GeoExtent <- c(LatMin, LonMin, LatMax, LonMax)


# follow link for token
res <- ncdc_stations(extent = GeoExtent,
                     token = "ROPsBcXPMsXQKnXTkBWmTBlToXdRvsqz")
leaflet(data = res$data) %>% 
  addTiles() %>% 
  leaflet::addCircleMarkers(~longitude, ~latitude)


ncdc_datasets(datasetid = "GHCND:GM000003218",
              token = "ROPsBcXPMsXQKnXTkBWmTBlToXdRvsqz")


# Find all stations in the bounding box and filter by
# starting date
Stations <- ncdc_stations(datasetid = 'GHCND',
                          locationid = CountryID,
                          limit = 1000, token = 'ROPsBcXPMsXQKnXTkBWmTBlToXdRvsqz',
                          extent = GeoExtent)$data %>%
  filter(year(mindate) < StartYear &
           year(maxdate) > EndYear)

# Next, check which of the stations has the data we need
# This creates a vector of stations that do have the 
# minimum amount of data we need and then filters
# the object down to those. 
DataTypesIndex <- numeric()
for(x in seq_len(dim(Stations)[1])) {
  Data <- ncdc_datatypes(datasetid = 'GHCND',
                         stationid = Stations$id[x])$data$id
  if(all(c('PRCP', 'TMAX', 'TMIN') %in% Data)) {
    DataTypesIndex <- c(DataTypesIndex, x)
  }
}

Stations <- Stations[DataTypesIndex, ]

# Calculate distances between study site and all stations
# in bounding box, choose the closest one
ClosestStation <- distGeo(p1 = c(AllSpecies$Longitude[i],
                                 AllSpecies$Latitude[i]),
                          p2 = matrix(cbind(Stations$longitude,
                                            Stations$latitude),
                                      ncol = 2)) %>%
  which.min() %>%
  Stations$id[.]

# Also want to save the distance for later
Dist <- distGeo(p1 = c(AllSpecies$Longitude[i],
                       AllSpecies$Latitude[i]),
                p2 = matrix(cbind(Stations$longitude,
                                  Stations$latitude),
                            ncol = 2)) %>%
  .[which.min(.)] 

# Next, download all the data we need!
# If the previous iteration already downloaded data
# for the same site and time frame, then there's no 
# point in doing it again.
if(str_detect(OutData$Source[listInd - 1], 
              ClosestStation) & 
   StartYear == (AllSpecies$StartYear[i - 1] - 30) &
   EndYear == (AllSpecies$EndYear[i - 1])) {
  
  message('Closest station, start and end years are the same. Reusing\n',
          'data from previous iteration')
  WeatherData <- WeatherData
} else {
  # However, if we haven't already downloaded it, then we will
  # again
  for(Year in StartYear:EndYear){
    
    
    if(Year == StartYear) {
      WeatherData <- ncdc(datasetid = 'GHCND',
                          stationid = ClosestStation,
                          startdate = paste0(Year, '-01-01'),
                          enddate = paste0(Year, '-12-31'),
                          limit = 1000)
    } else {
      YearData <- ncdc(datasetid = 'GHCND',
                       stationid = ClosestStation,
                       startdate = paste0(Year, '-01-01'),
                       enddate = paste0(Year, '-12-31'),
                       limit = 1000)
      
      WeatherData <- ncdc_combine(YearData, WeatherData)
      prog <- which(StartYear:EndYear == Year) / 
        length(StartYear:EndYear) * 100
      
      if(prog %% 10 < .1) {
        message(prog, '% of data downloaded')
      }
    } # End year loop
  } # end weather station
}
# Extract data from WeatherData object

# source for transformations :
# https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

# Some of these data sets don't have snow data (equatorial weather 
# stations don't get much ;) ), so here it splits based on
# presence of the variable in the data set.

if(!'SNOW' %in% unique(WeatherData$data$datatype)) {
  VarIndex <- c('PRCP', 'TMAX', 'TMIN')
  ClimData <- WeatherData$data %>%
    filter(datatype %in% VarIndex) %>%
    spread(key = 'datatype', value = 'value') %>%
    mutate(Year = year(date),
           Month = month(date)) %>%
    group_by(Year) %>%
    summarise(AnnPrecip = sum(PRCP / 10, # PRCP is 10ths of mm
                              na.rm = TRUE) / 10, # convert to cm
              MeanTemp = (mean(c(TMAX, TMIN),
                               na.rm = TRUE) / 10) , # T reported in 10ths of degC
              MaxTemp = (max(TMAX,
                             na.rm = TRUE) / 10),
              MinTemp = (min(TMIN,
                             na.rm = TRUE) / 10))
  
  # Extract monthly maximum precipitation and then join
  # it to the other climate data
  ClimData <- WeatherData$data %>%
    filter(datatype %in% VarIndex) %>%
    spread(key = 'datatype', value = 'value') %>%
    mutate(Year = year(date),
           Month = month(date)) %>%
    group_by(Year, Month) %>%
    summarise(MPrecip = sum(PRCP / 10, # PRCP  is in 10ths of mm
                            na.rm = TRUE) / 10) %>% # convert to cm
    ungroup() %>%
    group_by(Year) %>%
    summarise(MaxPrecip = max(MPrecip, na.rm = TRUE)) %>%
    left_join(ClimData, ., by = 'Year')
  
} else {
  
  VarIndex <- c('PRCP', 'SNOW', 'TMAX', 'TMIN')
  
  ClimData <- WeatherData$data %>%
    filter(datatype %in% VarIndex) %>%
    spread(key = 'datatype', value = 'value') %>%
    mutate(Year = year(date),
           Month = month(date)) %>%
    group_by(Year) %>%
    summarise(AnnPrecip = sum(PRCP / 10, # PRCP is 10ths of mm
                              SNOW / 10,      # snow is in mm. 10 mm snow  = 1 mm water
                              na.rm = TRUE) / 10, # convert to cm
              MeanTemp = (mean(c(TMAX, TMIN),
                               na.rm = TRUE) / 10) , # T reported in 10ths of degC
              MaxTemp = (max(TMAX,
                             na.rm = TRUE) / 10),
              MinTemp = (min(TMIN,
                             na.rm = TRUE) / 10))
  
  ClimData <- WeatherData$data %>%
    filter(datatype %in% VarIndex) %>%
    spread(key = 'datatype', value = 'value') %>%
    mutate(Year = year(date),
           Month = month(date)) %>%
    group_by(Year, Month) %>%
    summarise(MPrecip = sum(PRCP / 10,
                            SNOW / 10,
                            na.rm = TRUE) / 10) %>%
    ungroup() %>%
    group_by(Year) %>%
    summarise(MaxPrecip = max(MPrecip, na.rm = TRUE)) %>%
    left_join(ClimData, .)
  
  
}

# css selecotor for ARPA site

# dropdown menu to select data
#grafico_semplice_parametri_tutti_sel3

# Click this to go to station
# .item:nth-child(1) a
# 
