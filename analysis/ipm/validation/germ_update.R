rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)


germ   <- read_xlsx('data/seedbaskets.xlsx')

g_test <- germ %>% 
            mutate( g0_calc = seedlings2009 /seeds2008,
                    g1_calc = seedlings2010 /seeds2008,
                    g2_calc = seedlings2011 /seeds2008 ) 

expect_true( all(g_test$g0 == g_test$g0_calc) )
expect_true( all(g_test$g1 == g_test$g1_calc) )
expect_true( all(g_test$g2 == g_test$g2_calc) )

g_updt <- germ %>% 
            mutate( seeds2009 = seeds2008 - seedlings2009,
                    seeds2010 = seeds2008 - (seedlings2009 + seedlings2010) ) %>% 
            mutate( g1_update = seedlings2010 / seeds2009,
                    g2_update = seedlings2011 / seeds2009 )


# compare updated germination rates with observed ones
g_updt %>% 
  select(g0, g1_update, g2_update) %>% 
  colMeans

g_updt %>% 
  select(g0, g1, g2) %>% 
  colMeans
