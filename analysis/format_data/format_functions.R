# monthly climate form ----------------------------------------------------------------------------

# detrend based on SEASONS
prism_seas_form <- function(clim_x, clim_var = "ppt", yearz, m_back, m_obs){
  
  # seasons as suggested by Pardini
  seas_df <- data.frame( month_num = 1:12,
                         season    = c(rep('wet_t1',  2),
                                       rep('windy',3),
                                       rep('dry',3),
                                       rep('rainish',3),
                                       'wet_t0') )
  
  # Create 'climate of the year'
  # TRY ALTERNATIVE: split clim_x in two separate objects (t0 and t1)
  clim_s  <- clim_x %>% 
                left_join( seas_df ) %>% 
                select( -clim_var, -month_num ) %>% 
                # wet_t0, dry, and rainish is the following year. 
                # hence, for them, year = year + 1 
                mutate( year = replace(year, 
                                       season %in% c('wet_t0', 'dry', 'rainish'),
                                       year[season %in% c('wet_t0','dry', 'rainish')] + 1) 
                        ) %>% 
                mutate( season = replace(season,
                                         season == 'wet_t0',
                                         'wet_t1') )
  
  # seasonal means/sums
  if( clim_var == "ppt" ){
    clim_m <- clim_s %>% 
                group_by( year, season ) %>% 
                summarise( clim_value = sum(clim_value,na.rm=T) ) %>% 
                ungroup %>% 
                spread( season, clim_value)
  }else{
    clim_m <- clim_s %>% 
                group_by( year, season ) %>% 
                summarise( clim_value = mean(clim_value,na.rm=T) ) %>% 
                ungroup %>% 
                spread( season, clim_value)
  }
  
  # create anomalies from the spreaded months
  d_prec <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>%
              as.data.frame %>%
              bind_cols( clim_m[,"year", drop=F] ) %>%
              dplyr::select( year, dry:windy )
  
  clim_detr <- d_prec
  yr_range  <- range(yearz)
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather(month, clim_value, dry:windy ) %>%
                  setNames(c("year", "season", "clim_value") ) %>% 
                  mutate(season = factor(season, 
                                         levels = c('dry',
                                                    'rainish',
                                                    'wet_t1',
                                                    'windy') )
                         ) %>% 
                  arrange(year, season)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$season == 'windy')
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(yearz, clim_back, m_obs, long_out)
  x_clim    <- year_by_month_mat(clim_x_l, yearz)
  return(x_clim)

}

# detrend population-level climate; put it in "long" form
oni_form <- function(clim_x, yearz, m_back, m_obs){
  
  # "spread" the 12 months. These are already anomalies
  clim_m <- clim_x %>%  
              select( -clim_var, -mon ) %>% 
              spread( month_num, clim_value)
  
  clim_detr <- clim_m
  yr_range  <- range(yearz)
  
  # detrended climate in "long" form
  long_out  <- clim_m %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather_('month', 'clim_value', '1':'12') %>%
                  setNames(c("year", "month_num", "clim_value") ) %>% 
                  mutate(month_num = factor(month_num, levels = paste0(1:12)) ) %>% 
                  arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(yearz, clim_back, m_obs, long_out)
  x_clim    <- year_by_month_mat(clim_x_l, yearz)
  return(x_clim)

}


# detrend population-level climate; put it in "long" form
prism_clim_form <- function(clim_x, clim_var = "ppt", yearz, m_back, m_obs){
  
  # "spread" the 12 months
  clim_m <- clim_x %>%  
              select( -clim_var ) %>% 
              spread( month_num, clim_value)
  
  # create anomalies from the spreaded months
  d_prec <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>%
              as.data.frame() %>%
              bind_cols( clim_m[,"year", drop=F] ) %>%
              dplyr::select( c("year", 1:12) )
  
  clim_detr <- d_prec
  yr_range  <- range(yearz)
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather('month', 'clim_value', '1':'12') %>%
                  setNames(c("year", "month_num", "clim_value") ) %>% 
                  mutate(month_num = factor(month_num, levels = paste0(1:12)) ) %>% 
                  arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(yearz, clim_back, m_obs, long_out)
  x_clim    <- year_by_month_mat(clim_x_l, yearz)
  return(x_clim)

}

# format spei data
spei_clim_form <- function(clim_x, yearz, m_back, m_obs){
  
  # "spread" the 12 months
  clim_m <- clim_x %>%  
              select( -clim_var ) %>% 
              spread( month_num, clim_value)
  
  clim_detr <- clim_m
  yr_range  <- range(yearz)
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather('month', 'clim_value', '1':'12') %>%
                  setNames(c("year", "month_num", "clim_value") ) %>% 
                  mutate(month_num = factor(month_num, levels = paste0(1:12)) ) %>% 
                  arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(yearz, clim_back, m_obs, long_out)
  x_clim    <- year_by_month_mat(clim_x_l, yearz)
  return(x_clim)

}




# detrend population-level climate; put it in "long" form
month_clim_form <- function(clim_x, clim_var = "precip", yearz, m_back, m_obs){
  
  # ONI and SOI data are already detrended: skip directly to output
  if( clim_var == 'oni' | clim_var == 'soi' ){
    
    long_out <- clim_x
    
    # select temporal extent
    clim_back <- function(yrz, m_obs, dat){
      id <- which(dat$year == yrz & dat$month_num == m_obs)
      r  <- c( id:(id - (m_back-1)) )
      return(dat[r,"clim_value"])
    }
    
    # climate data in matrix form 
    year_by_month_mat <- function(dat, years){
      do.call(rbind, dat) %>% 
        as.data.frame %>%
        tibble::add_column(year = years, .before=1)
    }
    
    # calculate monthly precipitation values
    clim_x_l  <- lapply(yearz, clim_back, m_obs, long_out)
    x_clim    <- year_by_month_mat(clim_x_l, yearz)
    return(x_clim)

  }
  
  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y" ) 
  
  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
    as.character %>%
    as.data.frame(stringsAsFactors=F) %>%
    separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
    bind_cols(clim_x) %>% 
    dplyr::select(-year,-day) %>%
    setNames( c("year", "month", "day", "species", "value") )
  
  # # if climate_var airt, then do means, otherwise, do sums! 
  if( clim_var == "airt" | clim_var == "soilmoist" ){
    clim_m  <- clim_d %>%
      group_by(year, month) %>%
      summarise( value = mean(value, na.rm=T) )  %>%
      spread( month, value ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame()
  } else{
    clim_m  <- clim_d %>%
      group_by(year, month) %>%
      summarise( value = sum(value, na.rm=T) )  %>%
      spread( month, value ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame()
  }
  
  # throw error
  if( !any( clim_var %in% c("precip","pet","airt","gdd","soilmoist")) ) {
    stop( paste0(clim_var," is not a supported varible") ) 
  }
  
  # detrend climate - but NOT if you are using GDD 
  if( clim_var != "gdd" ){
    d_prec   <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>%
      as.data.frame() %>%
      bind_cols( clim_m[,"year", drop=F] ) %>%
      dplyr::select( c("year", month.abb) )
  } else{
    d_prec   <- clim_m
  }
  
  # Present 
  nans_n <- d_prec[,-1] %>% as.matrix %>% is.nan %>% sum
  if(nans_n > 0) print( "Warning: NANs present in climate predictor" )
  
  # Make NaNs 0
  for(c_i in 1:ncol(d_prec) ){
    d_prec[,c_i] <- replace(d_prec[,c_i], is.nan(d_prec[,c_i]), 0)
  }
  
  # fecth year range, observation month
  clim_detr <- d_prec
  yr_range  <- range(yearz)
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
    subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
    gather(month, precip, Jan:Dec) %>%
    setNames(c("year", "month", "clim_value") ) %>% 
    mutate(month_num = factor(month, levels = month.abb) ) %>% 
    mutate(month_num = as.numeric(month_num) ) %>% 
    arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrz, m_obs, dat){
    id <- which(dat$year == yrz & dat$month_num == m_obs)
    r  <- c( id:(id - (m_back-1)) )
    return(dat[r,"clim_value"])
  }
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  clim_x_l  <- lapply(yearz, clim_back, m_obs, long_out)
  x_clim    <- year_by_month_mat(clim_x_l, yearz)
  return(x_clim)
  
}



# format climate information based on climate variable, and summary function
clim_form <- function(clim, clim_variable, foo, years, m_obs, m_back){
  
  # detrend data ----------------------------------------------------------
  clim_x    <- subset(clim, clim_var == clim_variable)
  
  # clim_x value rename
  years_vec <- clim$year %>% unique
  
  # prepare indexes
  split_id  <- Filter(function(x) x < 366, which( (0:365 %% 5) == 0 ) ) %>%
                lapply( function(x) x + (0:4) ) %>% 
                setNames( paste0(c(1:length(.))) )
  
  # period summaries
  period_summ <- function(ii,clim_yr, foo){
    
    #fun         <- substitute( foo )
    period_subs <- clim_yr[ii,]
    period_clim <- data.frame( year    = unique(period_subs$year), 
                               value   = foo(period_subs$value, na.rm=T) %>% eval
    )
    return(period_clim)
    
  }
  
  # repeat summaries by year
  years_pop <- function(yrs, clim_x, period_summ, split_id, foo){
    
    clim_yr    <- subset(clim_x, year == yrs)
    # create summaries data frame
    period_l  <- lapply(split_id, period_summ, clim_yr, foo)
    period_df <- Reduce(function(...) rbind(...), period_l) %>% 
                    tibble::add_column( period = 1:length(split_id), .after=2 ) %>%
                    spread( period, value )
    return(period_df)
    
  }
  
  # substitute -Inf with NA
  na_for_inf  <- function(x) replace(x, x==-Inf, NA)
  
  # 5 day interval summaries
  clim_int_l  <- lapply(years_vec, years_pop, clim_x, period_summ, split_id, sum)
  clim_int    <- Reduce(function(...) rbind(...), clim_int_l) %>%
                    apply(2, na_for_inf)
  
  # detrend climate 
  d_prec      <- apply(clim_int[,-1], 2, FUN = scale, center = T, scale = T) %>%
                    as.data.frame() %>%
                    bind_cols( as.data.frame(clim_int[,"year", drop=F]) ) %>%
                    dplyr::select( c("year", 1:(ncol(clim_int)-1) ) )
  
  # Present 
  nans_n <- d_prec[,-1] %>% as.matrix %>% is.na %>% sum
  if(nans_n > 0) print( "Warning: NANs present in climate predictor" )
  
  # Make NAs 0
  for(c_i in 1:ncol(d_prec) ){
    d_prec[,c_i] <- replace(d_prec[,c_i], is.na(d_prec[,c_i]), 0)
  }
  
  
  # data in long form -----------------------------------------------------------
  
  # fecth year range, observation month
  yr_range  <- range(years)
  
  # format day one (just a trick)
  day_one   <- as.Date( paste0("1/1/", first(d_prec$year) ), 
                        format="%d/%m/%Y" ) 
  
  # dates data frame (with julian date)
  dates_df  <- data.frame( date = as.Date(0:364, day_one) ) %>%
                  mutate( date  = as.character(date) ) %>%
                  mutate( jul   = 1:365 ) %>%
                  separate( date, into = c("year","month","day"), sep = "-") %>%
                  mutate( month = as.numeric(month) )
  
  # identify what's the 5-day period where month observation is located 
  five_day_id <- function(ii, x){
    # pick last day of the month
    out <- subset(x, month == ii) %>% .$jul %>% last  
    # number of 5 day period
    return( round(out / 5) ) 
  }
  
  # 5-day period of month (demographic) observation
  five_d_id   <- five_day_id(m_obs, dates_df)
  
  # N. of 5-day periods "back"
  five_d_back <- round( c(m_back * 30.41667) / 5)
  
  # detrended climate in "long" form
  long_out  <- d_prec %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather(five_d, clim_value, paste0(1:73)) %>%
                  mutate(five_d = as.numeric(five_d) ) %>% 
                  arrange(year, five_d)
  
  # select temporal extent
  clim_back <- function(yrz, long_dat, five_d_id, five_d_back){
                  id <- which(long_dat$year == yrz & long_dat$five_d == five_d_id)
                  r  <- c( id:(id - (five_d_back-1)) )
                  return(long_dat[r,"clim_value"])
  }
  clim_x_l  <- lapply(years, clim_back, long_out, five_d_id, five_d_back)
  
  # climate data in matrix form 
  year_by_month_mat <- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # calculate monthly precipitation values
  x_clim    <- year_by_month_mat(clim_x_l, years)
  
  return(x_clim)
  
}

