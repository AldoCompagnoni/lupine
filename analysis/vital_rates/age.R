# calculate average age at death
# calculate average age at death for non seedlings

# read data
lupine_df   <- read.csv( "data/lupine_all.csv")

# average age at death
trans_d %>% 
  # remove data prior to 2012
  # (prop. of population with age info below 75%)
  subset( year > 2011 ) %>% 
  subset( !is.na(age_t0) ) %>% 
  # subset( age_t0 > 1 ) %>% 
  .$age_t0 %>% mean

# average age at death of adults (past 1 year of age)
trans_d %>% 
  # remove data prior to 2012
  # (prop. of population with age info below 75%)
  subset( year > 2011 ) %>% 
  subset( !is.na(age_t0) ) %>% 
  subset( age_t0 > 1 ) %>%
  .$age_t0 %>% mean
