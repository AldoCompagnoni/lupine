# attempt to 
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
source("analysis/format_data/format_functions.R")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)

# data
clim        <- read.csv("data/prism_point_reyes_87_18.csv")

# format climate data ----------------------------------------
years     <- c(1990:2018)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var ){
  
  # set names of climate variables
  clim_names <- paste0( var,c('_t0','_tm1','_t0_tm1','_t0_tm2') )
  
  mutate(x, 
         avgt0     = x %>% select(V1:V12) %>% rowSums,
         avgtm1    = x %>% select(V13:V24) %>% rowSums,
         avgt0_tm1 = x %>% select(V1:V24) %>% rowSums,
         avgt0_tm2 = x %>% select(V1:V36) %>% rowSums ) %>% 
    select(year, avgt0, avgtm1, avgt0_tm1, avgt0_tm2) %>% 
    setNames( c('year',clim_names) )
  
}

# format climate - need to select climate predictor first 
tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              year_anom('tmp') %>% 
              select(year,tmp_tm1)

# calculate empirically observed var-covar  
tmp_df <- cbind(tmp_mat[1:28,2],tmp_mat[2:29,2])
Sig    <- cov(tmp_df)

# # from skratch 
# corm   <- cor(tmp_df)
# sig    <- apply(tmp_df,2,var) %>% sqrt
# Sig    <- (sig * corm * sig) 

library(MASS)
xx     <- mvrnorm(1000, mu = c(0,0), Sigma = Sig)
cov(xx)
cor(xx)

# concatenate observations by row 
get_row <- function(ii) xx[ii,1:2]
vec_xx  <- sapply(1:nrow(xx), get_row) %>% as.numeric    

# correlation is weaker
acf(vec_xx)
acf(tmp_mat$tmp_tm1)



# https://stats.stackexchange.com/questions/29239/creating-auto-correlated-random-values-in-r
# https://stat.ethz.ch/pipermail/r-help/2007-April/128925.html
# 
# # create the initial x variable
# x1 <- rnorm(100, 15, 5)
# 
# # x2, x3, and x4 in a matrix, these will be modified to meet the
# # criteria
# x234 <- scale(matrix( rnorm(300), ncol=3 ))
# 
# # put all into 1 matrix for simplicity
# x1234 <- cbind(scale(x1),x234)
# 
# # find the current correlation matrix
# c1 <- var(x1234)
# 
# # cholesky decomposition to get independence
# chol1 <- solve(chol(c1))
# 
# # what is this?!?
# newx <-  x1234 %*% chol1 
# 
# # check that we have independence and x1 unchanged
# zapsmall(cor(newx))
# all.equal( x1234[,1], newx[,1] )
# 
# # create new correlation structure
# # (zeros can be replaced with other r vals)
# newc <- matrix( 
# c(1  , 0.4, 0.5, 0.6, 
#   0.4, 1  , 0  , 0  ,
#   0.5, 0  , 1  , 0  ,
#   0.6, 0  , 0  , 1  ), ncol=4 )
# 
# # check that it is positive definite
# eigen(newc)
# 
# chol2 <- chol(newc)
# 
# finalx <- newx %*% chol2 * sd(x1) + mean(x1)
# 
# # verify success
# mean(x1)
# colMeans(finalx)
# 
# sd(x1)
# apply(finalx, 2, sd)
# 
# zapsmall(cor(finalx))
# pairs(finalx)
# 
# all.equal(x1, finalx[,1])
