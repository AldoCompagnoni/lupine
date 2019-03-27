# Simulation and data example to check the bias of random effects on the "overall mean"
# This is idea was pointed out by Crone 2013 (AmNat) and Crone 2016 (JoE)

# simulation --------------------------------------

# simulate a hierarchical Poisson process
set.seed(1779)

# 30 normally distributed random effects
ran_ef <- rnorm(30,2,1)
y      <- list()
for(i in seq_along(ran_ef) ){
  y[[i]] <- rpois(30,exp(ran_ef[i]))
}
y <- unlist(y)

# put everything in a data frame for analysis
df_y <- expand.grid( n  = 1:30,
                     re = 1:30 ) %>% 
          mutate(    y  = y )

# random effect and fixed effects GLMs
re_mod <- glmer( y ~ 1 + (1 | re), data=df_y, family='poisson')
fe_mod <- glm( y ~ 1, data=df_y,family='poisson')

# mean effects differ substantially between two models
fixef(re_mod) 
coef(fe_mod)

# calculate mean of re_mod via back-transformation
coef(re_mod)$re[,1] %>% 
  exp() %>% 
  mean %>% 
  log
coef(fe_mod)

# bias in predictions is enormous
mean(y)
coef(fe_mod) %>% exp
fixef(re_mod) %>% exp

# but not on the mean of the single levels effects
df_y %>% 
  subset( re==1) %>% 
  .$y %>% mean
coef(re_mod)$re[1,1] %>% exp

# data example (lupine fertility) ----------------------

setwd("C:/cloud/Dropbox/lupine")
source("analysis/format_data/format_functions.R")

# lupine data
lupine_df   <- read.csv( "data/lupine_all.csv") #%>% 

# get fertility data
fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) 

# Compare glm with glmm
mod_glm   <- glm(numrac_t0 ~ log_area_t0, 
                  data=fert, family='poisson')
mod_glmm  <- glmer(numrac_t0 ~ log_area_t0 + (1|year), 
                  data=fert, family='poisson')

# get mean size, sd, replication (n)
siz_mean <- subset(fert, !(is.na(log_area_t0) | is.na(numrac_t0)) ) %>% 
              .$log_area_t0 %>% mean
siz_sd   <- subset(fert, !(is.na(log_area_t0) | is.na(numrac_t0)) ) %>% 
              .$log_area_t0 %>% sd
siz_n    <- subset(fert, !(is.na(log_area_t0) | is.na(numrac_t0)) ) %>% 
              .$log_area_t0 %>% length

# simulate mean n. of racemes per indiv 
# we need to simulate this because of Jensen's inequality
set.seed(1776)
exp(coef(mod_fr)[1] + coef(mod_fr)[2]*rnorm(6131,5.934136,0.9072458) ) %>% mean
set.seed(1776)
exp(fixef(mod_fr1)[1] + fixef(mod_fr1)[2]*rnorm(6131,5.934136,0.9072458)) %>% mean

# observed mean n. of racemes per indiv
subset(fert, !(is.na(log_area_t0) | is.na(numrac_t0)) ) %>% 
  .$numrac_t0 %>% mean


# plot differences --------------------------------------

# produce predictions
x_seq     <- seq(min(fert$log_area_t0),
                 max(fert$log_area_t0), length.out=100 )
new_df <- expand.grid( log_area_t0 = x_seq)

par(mfrow=c(1,1),mar=c(3,3,0.1,0.1),mgp=c(1.8,0.7,0))
plot(numrac_t0 ~ log_area_t0, data=fert)
y_glm <- predict(mod_glm, 
                 newdata=new_df,
                 type='response')
lines(x_seq,y_glm, lwd=2, col='red')
y_glmm <- exp(fixef(mod_glmm)[1] + fixef(mod_glmm)[2]*x_seq)
lines(x_seq, y_glmm, lwd=2, col='blue')
