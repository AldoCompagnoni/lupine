setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(lme4)
options(stringsAsFactors = F)
source("analysis/format_data/format_scripts.R")
source("analysis/format_data/format_functions.R")


# data
lupine_df   <- read.csv("data/lupine_all.csv")
enso        <- read.csv("data/enso_data.csv")
# clim        <- read.csv("data/lupine_fc_vars.csv")
clim        <- read.csv("data/prism_point_reyes_87_18.csv")


# SURVIVAL ----------------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( stage_t0 != 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 )

# climate format ----------------------------------------------------------------
years     <- unique(surv$year)
m_obs     <- 5
m_back    <- 36

# format climate - need to select climate predictor first 
clim_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs)
# clim_mat <- subset(enso, clim_var == "oni") %>%
#               oni_form(years, m_back, m_obs)


# seedling data
surv_clim <- surv %>% 
                # add 1 to get year_t1
                left_join( clim_mat ) %>%
                subset( !is.na(location) )
                # indices for STAN models

surv_df   <- surv_clim %>% 
                mutate( year_i     = year %>% as.factor %>% as.numeric,
                        site_i     = location %>% as.factor %>% as.numeric,
                        avgt0      = surv_clim %>% select(V1:V12) %>% rowSums,
                        avgtm1     = surv_clim %>% select(V13:V24) %>% rowSums,
                        avgtm2     = surv_clim %>% select(V25:V36) %>% rowSums,
                        avgt0_tm1  = surv_clim %>% select(V1:V24) %>% rowSums,
                        avgt0_tm2  = surv_clim %>% select(V1:V36) %>% rowSums, 
                        avgtm1_tm2 = surv_clim %>% select(V13:V36) %>% rowSums )

# climate data
clim_pred <- dplyr::select( surv_df, year ) %>%
                inner_join( clim_mat ) %>% 
                unique %>%
                arrange( year ) %>% 
                dplyr::select( -year )

# test correspondence between years in clim_pred and years in surv_df.
surv_clim %>% 
  mutate( year_i  = year %>% as.factor %>% as.numeric,
          year_r  = year %>% as.factor ) %>% 
  select( year_i, year_r ) %>% 
  unique

# Random effects model selection -------------------------------------------------
mod             <- glmer(surv_t1 ~ log_area_t0 + (1 | year_i), data=surv_df, family='binomial')
mod2            <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (1 | year_i), data=surv_df, family='binomial')
mod_size        <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 | year_i), data=surv_df, family='binomial')
mod_size2       <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t02 | year_i), data=surv_df, family='binomial')
mod_size_size2  <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + (log_area_t0 + log_area_t02 | year_i), data=surv_df, family='binomial')

# best model is most complex
AIC(mod,
    mod2,
    mod_size,
    mod_size2,
    mod_size_size2)

# random individual effect
mod_size        <- glmer(surv_t1 ~ log_area_t0 + 
                                   log_area_t02 + 
                                    (log_area_t0 | year_i) + (1|newid), data=surv_df, family='binomial')

# Climate effect model selection -------------------------------------------------
mod_t0      <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + avgt0      + (log_area_t0 | year_i) + (1 | site_i), data=surv_df, family='binomial')
mod_tm1     <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + avgtm1     + (log_area_t0 | year_i) + (1 | site_i), data=surv_df, family='binomial')
mod_tm2     <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + avgtm2     + (log_area_t0 | year_i) + (1 | site_i), data=surv_df, family='binomial')
mod_t0_tm1  <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + avgt0_tm1  + (log_area_t0 | year_i) + (1 | site_i), data=surv_df, family='binomial')
mod_tm1_tm2 <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + avgtm1_tm2 + (log_area_t0 | year_i) + (1 | site_i), data=surv_df, family='binomial')
mod_t0_tm2  <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + avgt0_tm2  + (log_area_t0 | year_i) + (1 | site_i), data=surv_df, family='binomial')


relgrad0 <- with(mod_t0@optinfo$derivs,solve(Hessian,gradient))
relgrad1 <- with(mod_tm1@optinfo$derivs,solve(Hessian,gradient))
relgrad2 <- with(mod_tm2@optinfo$derivs,solve(Hessian,gradient))
relgrad0_1 <- with(mod_t0_tm1@optinfo$derivs,solve(Hessian,gradient))
relgrad1_2 <- with(mod_tm1_tm2@optinfo$derivs,solve(Hessian,gradient))
relgrad0_2 <- with(mod_t0_tm2@optinfo$derivs,solve(Hessian,gradient))

max(abs(relgrad0))
max(abs(relgrad1))
max(abs(relgrad2))
max(abs(relgrad0_1))
max(abs(relgrad1_2))
max(abs(relgrad0_2))


# best model is most complex
AIC(mod_t0,
    mod_tm1,
    mod_tm2,
    mod_t0_tm1,
    mod_tm1_tm2,
    mod_t0_tm2)


# GROWTH ---------------------------------------------------------------------------
grow        <- lupine_df %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1 = log(area_t1),
                          log_area_t0 = log(area_t0),
                          year        = year + 1 ) 

# growth data
grow_clim <- left_join(grow, clim_mat) %>%
                subset( !is.na(location) )
                # indices for STAN models
grow_df   <- grow_clim %>% 
                mutate( year_i = year %>% as.factor %>% as.numeric,
                        site_i = location %>% as.factor %>% as.numeric,
                        avgt0   = grow_clim %>% select(V1:V12) %>% rowSums,
                        avgtm1  = grow_clim %>% select(V13:V24) %>% rowSums,
                        avgtm2  = grow_clim %>% select(V25:V36) %>% rowSums,
                        avgt0_tm1  = surv_clim %>% select(V1:V24) %>% rowSums,
                        avgt0_tm2  = surv_clim %>% select(V1:V36) %>% rowSums, 
                        avgtm1_tm2 = surv_clim %>% select(V13:V36) %>% rowSums)

# Random effects model selection -------------------------------------------------
mod         <- lmer(log_area_t1 ~ log_area_t0 + (1 | year_i), data=surv_df)
mod_size    <- lmer(log_area_t1 ~ log_area_t0 + (log_area_t0 | year_i), data=surv_df)
mod_st      <- lmer(log_area_t1 ~ log_area_t0 + (1 | year_i) + (1 | site_i), data=surv_df)
mod_size_st <- lmer(log_area_t1 ~ log_area_t0 + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)


# best model is most complex
AIC(mod,mod_size,mod_st,mod_size_st)


# Climate effect model selection -------------------------------------------------
mod_t0      <- lmer(log_area_t1 ~ log_area_t0 + avgt0      + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)
mod_tm1     <- lmer(log_area_t1 ~ log_area_t0 + avgtm1     + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)
mod_tm2     <- lmer(log_area_t1 ~ log_area_t0 + avgtm2     + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)
mod_t0_tm1  <- lmer(log_area_t1 ~ log_area_t0 + avgt0_tm1  + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)
mod_tm1_tm2 <- lmer(log_area_t1 ~ log_area_t0 + avgtm1_tm2 + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)
mod_t0_tm2  <- lmer(log_area_t1 ~ log_area_t0 + avgt0_tm2  + (log_area_t0 | year_i) + (1 | site_i), data=surv_df)

# best model is most complex
AIC(mod_t0,
    mod_tm1,
    mod_tm2,
    mod_t0_tm1,
    mod_tm1_tm2,
    mod_t0_tm2)
