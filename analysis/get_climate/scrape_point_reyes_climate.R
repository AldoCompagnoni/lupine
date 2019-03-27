# Rselenium Scrape of Point Reyes climate data 
# https://www.r-bloggers.com/scraping-with-selenium/
setwd('C:/cloud/Dropbox/lupine')
library(RSelenium) 
library(rvest)
library(dplyr)
library(XML)
library(testthat)

# devtools::install_github('ropensci/RSelenium')

# start Rselenium?
# rD        <- rsDriver( browser = "firefox", port = 4565L)
# remDr     <- remoteDriver( remoteServerAddr = "localhost",
#                            browser = "firefox", port = 4565L )

# start Rselenium?

# which months you need to get
var_df <- expand.grid( months = month.name,
                       year   = 2007:2009,
                       stringsAsFactors = F)



rD       <- rsDriver(port=4444L,browser="firefox")
remDr    <- rD$client
clim_ddd  <- list()

var_df <- expand.grid( months = month.name,
                       year   = 2006:2009,
                       stringsAsFactors = F) %>% subset(year == 2006)


for(ii in 10:12){

  # open up climate page
  open_page <- function(){
    # open connection to page
    remDr$navigate( 'https://wrcc.dri.edu/cgi-bin/wea_monsum.pl?nvprca' )
    Sys.sleep(2)
    
    # select month
    month <- remDr$findElement(using = "xpath", '/html/body/form/select[1]')
    m_in  <- var_df$months[ii]
    month$sendKeysToElement( list(m_in) )
    
    # select year
    year    <- remDr$findElement(using = "xpath", '/html/body/form/select[2]')
    year_in <- var_df$year[ii] %>% as.character
    year$sendKeysToElement( list(year_in) )
    
    
    # get the data table
    button <- remDr$findElement(using = "xpath", "/html/body/form/input[2]")
    button$sendKeysToElement(list(key = "enter"))
    Sys.sleep(15)
    
    # store title name
    tab_title  <- remDr$findElement(using = "xpath", '/html/body/center[2]/h3')
    month_year <- tab_title$getElementAttribute("outerHTML")[[1]] %>% 
                    gsub('<h3>| </h3>','',. ) %>% 
                    strsplit(', ') %>% unlist %>% return()
  }

  # CHECK whether the page opened is correct (RSelenium can malfunction!)
  test_vec <- open_page()
  log_vec  <- c(test_vec[2] != var_df$year[ii],test_vec[1] != var_df$month[ii])
  
  while( any(log_vec) ){
    test_vec <- open_page()
    log_vec  <- c(test_vec[2] != var_df$year[ii],test_vec[1] != var_df$month[ii])
  }
      
  # read the table
  tab       <- remDr$findElement(using = "xpath", "/html/body/center[3]/table[1]")
  # Sys.sleep(1)
  elemtxt   <- tab$getElementAttribute("outerHTML")[[1]] # gets us the HTML
  # Sys.sleep(1)
  elemxml   <- htmlTreeParse(elemtxt, useInternalNodes=T) # parse string into HTML tree to allow for querying with XPath
  # Sys.sleep(1)
  month_raw <- readHTMLTable(elemxml)[[1]]
  r_row     <- which(month_raw == 'MONTHLY STATISTICS'):nrow(month_raw)
  cols_keep <- c(1,2, which(month_raw[2,]=='Mean')[1],
                      which(month_raw[2,]=='Precip.') )
  month_d   <- month_raw[-r_row,cols_keep] %>% 
                 setNames( c('month_day', 'julian_day', 'meant', 'ppt') )
  
  clim_ddd[[ii]] <- month_d %>% 
                      mutate( mon = test_vec[1],
                              yr  = test_vec[2] ) 
}


clim_out <- lapply(clim_ddd, function(x) x[-c(1:3),]) %>% 
              bind_rows

# test there are enough days
expect_true(grepl('365|366',clim_out$julian_day %>% unique %>% length %>% as.character) )
setdiff(1:365, clim_out$julian_day %>% unique)
clim_out$yr %>% unique
yr_out <- clim_out$yr %>% unique
 
# introduce NAs
clim_out <- clim_out %>% 
              mutate( meant = replace(meant, meant == '', NA),
                      ppt   = replace(ppt,   ppt == '', NA) )
                      
# write it out
write.csv(clim_out, 
          paste0('data/weather_station/point_reyes_',
                 substr(yr_out,3,4),
                 '.csv'),
          row.names=F)
