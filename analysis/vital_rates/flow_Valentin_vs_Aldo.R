rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(ggplot2)
library(testthat)
library(bbmle)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")

# data
lupine_df   <- read.csv( "data/lupine_all.csv")
enso        <- read.csv( "data/enso_data.csv")
clim        <- read.csv( "data/prism_point_reyes_87_18.csv")

# data format --------------------------------------------------------------



library(data.table)
library(dplyr)
library(ggplot2)
library(broom)


# Read data ---------------------------------------------------------------
dt_raw  <- read.csv('data/lupine_all.csv')
setDT(dt_raw)

# Add a copy of the table so that the last panel refers to the polled data
# across pupolations for each year.
# `copy` is needed because data.table modifies columns by reference.
dt <- rbindlist(list(dt_raw, copy(dt_raw)[, location := "pooled"]))

# Use an ordered factor forcing "pooled" label at the bottom.
dt[, location := factor(location, 
                        levels = c("AL (1)", 
                                   "ATT (8)",
                                   "BR (6)", 
                                   "BS (7)", 
                                   "DR (3)", 
                                   "NB (2)", 
                                   "POP9 (9)", 
                                   "pooled"),
                        ordered = TRUE)]

# Subset from table with pooled:
flow_v <- dt[stage_t0 != "SL" & !is.na(flow_t0) & area_t0 != 0 & location != "pooled",
            .(newid,location, year, area_t0, flow_t0)]

flow_a <- subset(lupine_df, !is.na(flow_t1) ) %>% 
          subset( area_t1 != 0 ) %>% 
          mutate( log_area_t1  = log(area_t1),
                  log_area_t12 = log(area_t1)^2,
                  year         = year +1  ) %>% 
          select(newid,location,year,area_t0,flow_t1)



props_a  <- flow_a %>% 
                  group_by( year, location ) %>% 
                  summarise( succ  = sum(flow_t1, na.rm=T),
                             fail  = sum(flow_t1 == 0, na.rm=T) ) %>% 
                  ungroup %>% 
                  mutate( tot  = succ + fail,
                          prop_a = round(succ/(succ+fail),3) ) %>% 
                  as.data.frame %>% 
                  select(year, location, prop_a )

props_v  <- flow_v %>% 
                  mutate( location = as.character(location) ) %>% 
                  group_by( year, location ) %>% 
                  summarise( succ  = sum(flow_t0, na.rm=T),
                             fail  = sum(flow_t0 == 0, na.rm=T) ) %>% 
                  ungroup %>% 
                  mutate( tot  = succ + fail,
                          prop_v = round(succ/(succ+fail),3) ) %>% 
                  as.data.frame %>% 
                  select(year, location, prop_v )


prop_compare <- full_join(props_a, props_v) %>% 
                  mutate( test_c = as.numeric(prop_a == prop_v) ) %>% 
                  subset( test_c == 0 | is.na(test_c) )

prop_compare

flow_a %>% 
  subset( stage_t1 == "SL" ) %>% 
  .$flow_t1 %>% 
  sum(na.rm=T)

