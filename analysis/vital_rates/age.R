# calculate average age at death
# calculate average age at death for non seedlings
library(dplyr)

# read data
lupine_df   <- read.csv( "data/lupine_all.csv")

# average age at death
lupine_df %>% 
  # remove data prior to 2012
  # (prop. of population with age info below 75%)
  subset( year > 2011 ) %>% 
  subset( !is.na(age_t0) ) %>% 
  # subset( age_t0 > 1 ) %>% 
  .$age_t0 %>% 
  mean

# average age at death of adults (past 1 year of age)
lupine_df %>% 
  # remove data prior to 2012
  # (prop. of population with age info below 75%)
  subset( year > 2007 ) %>% 
  subset( !is.na(age_t0) ) %>% 
  subset( age_t0 > 2 ) %>% 
  .$age_t0 %>% 
  mean


# histogram of dying individuals
age_at_death <- lupine_df %>% 
  # remove data prior to 2012
  # (prop. of population with age info below 75%)
  subset( year > 2007 ) %>% 
  subset( !is.na(age_t0) ) %>% 
  subset( surv_t1 == 0 ) %>% 
  .$age_t0 

# proportion of seeds dying at age...
tab_d <- age_at_death %>% table
(tab_d / sum(tab_d)) %>% round(3)
