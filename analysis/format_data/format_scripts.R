# detrend population-level climate; put it in "long" form
clim_form <- function(clim_x, clim_type = "prate", foo){ 
  
  # years
  years_vec <- clim_x$year %>% unique %>% sort
  
  # clim_x value rename
  clim_x    <- clim_x %>% base::subset( clim_var == clim_type )
  
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
  clim_int_l  <- lapply(years_vec, years_pop, clim_x, period_summ, split_id, foo)
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
  
  return(d_prec)
  
}
