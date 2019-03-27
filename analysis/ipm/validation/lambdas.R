# Create lambas
# rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"
# "AL (1)" is good
# "TAT (8)" is not


# data
lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                  subset( location %in% c("BS (7)") ) %>% 
                  subset( year > 2008 )



# lambdas for each location ----------------------------------
pop_n     <- lupine_df %>% 
                subset( year > 2008 & !is.na(stage_t0) ) %>% 
                count( year, location )

# matplot of population numbers
pop_n_mat <- spread(pop_n,location,n)
matplot(pop_n_mat$year,pop_n_mat[,-1],type='l',lty=1)


# pop numbers at each time step
n_t1 <- mutate(pop_n, year = year -1 ) %>% 
            rename( n_t1 = n )
n_t0 <- pop_n
n_df <- full_join(n_t0, n_t1)  

mean(log(n_df$n_t1 / n_df$n),na.rm=T) %>% exp



# lambdas for each location ----------------------------------
pop_n     <- lupine_df %>% 
                subset( year > 2007 & !is.na(stage_t0) ) %>% 
                count( year )

# pop numbers at each time step
n_t1 <- mutate(pop_n, year = year -1 ) %>% 
            rename( n_t1 = n )
n_t0 <- pop_n
n_df <- full_join(n_t0, n_t1)  

mean(log(n_df$n_t1 / n_df$n),na.rm=T) %>% exp
