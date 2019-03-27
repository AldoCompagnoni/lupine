# IPM from data
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(readxl)
library(testthat)
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"

all_indiv_sample <- c("BS (7)", 'DR (3)', 'NB (2)')

# data
lupine_df   <- read.csv( "data/lupine_all.csv") #%>% 
                  # subset( !(location %in% c("AL (1)")) ) 
                  # subset( location %in% c("BS (7)") ) %>%
                  # subset( year > 2008 )
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")


# vital rates format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( stage_t0 != 'SL' ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) 

seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2,
                          log_area_t03= log_area_t0^3 )

sl_grow     <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2 ) %>% 
                  # only seedling GROWTH
                  subset( surv_t1 == 1 )

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset( !(stage_t0 %in%  "SL") ) %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM", "SL")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2 )

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 )

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
                  subset( !(numrac_t0 %in% 0) )


# models ---------------------------------------------------------
mod_sl   <- glm(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, data=seedl, family='binomial')
mod_g_sl <- lm( log_area_t1 ~ log_area_t0, data=sl_grow)
mod_s    <- glm(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, data=surv, family='binomial')
mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range( c(grow$log_area_t0, grow$log_area_t1) )
mod_fl   <- glm(flow_t0 ~ log_area_t0, data=flow, family='binomial')
mod_fr   <- glm(numrac_t0 ~ log_area_t0, data=fert, family='poisson')
mod_fr   <- MASS::glm.nb(numrac_t0 ~ log_area_t0, data=fert )
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))
germ_coef<- select(germ, g0:g2) %>% colMeans


# vital rate models 
surv_sl_p <- coef(mod_sl)
grow_sl_p <- coef(mod_g_sl)
grow_sl_p <- c(grow_sl_p, summary(mod_g_sl)$sigma)
surv_p    <- coef(mod_s)
grow_p    <- coef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- coef(mod_fl)
fert_p    <- coef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ_coef



# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }


# list of mean IPM parameters. 
pars_mean   <- list( # seedlings vital rates
                     surv_sl_b0   = surv_sl_p['(Intercept)'],
                     surv_sl_b1   = surv_sl_p['log_area_t0'],
                     surv_sl_b2   = surv_sl_p['log_area_t02'],
                     surv_sl_b3   = surv_sl_p['log_area_t03'],
                     
                     grow_sl_b0   = grow_sl_p["(Intercept)"],
                     grow_sl_b1   = grow_sl_p["log_area_t0"],
                     grow_sl_sig  = grow_sl_p[3], 
                     
                     # adults vital rates           
                     surv_b0      = surv_p['(Intercept)'],
                     surv_b1      = surv_p['log_area_t0'],
                     surv_b2      = surv_p['log_area_t02'],
                     surv_b3      = surv_p['log_area_t03'],
                     
                     grow_b0      = grow_p['(Intercept)'],
                     grow_b1      = grow_p['log_area_t0'],
                     grow_sig     = grow_p[3],

                     flow_b0      = flow_p['(Intercept)'],
                     flow_b1      = flow_p['log_area_t0'],
                     
                     fert_b0      = fert_p['(Intercept)'],
                     fert_b1      = fert_p['log_area_t0'],
                     
                     abort        = 0.22, # hardcoded for now!
                     clip         = 0.57, # hardcoded for now!
                     
                     fruit_rac    = fr_rac_p,
                     seed_fruit   = seed_fr_p,
                     g0           = 0.035,#germ_coef['g0'],#
                     g1           = germ_coef['g1'],
                     g2           = germ_coef['g2'],
                     
                     recr_sz      = size_sl_p$mean_sl_size,
                     recr_sd      = size_sl_p$sd_sl_size,
                     
                     one_sz       = sl_grow$log_area_t1 %>% mean(na.rm=T),
                     one_sd       = sl_grow$log_area_t1 %>% sd(na.rm=T),
                     
                     L_sl         = size_sl_p$min_sl_size,
                     U_sl         = size_sl_p$max_sl_size,
                     
                     
                     L            = g_lim[1],
                     U            = g_lim[2],
                     
                     mat_siz_sl = 100,
                     mat_siz    = 100 )


# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# Transforms all values below/above limits in min/max size
x_range_s <- function(x,pars){
  pmin(pmax(x,pars$L_sl),pars$U_sl)
}

# Survival at size x
sx_s<-function(x,pars){
  xb <- x_range_s(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_sl_b0 + 
                    pars$surv_sl_b1 * xb + 
                    pars$surv_sl_b2 * (xb^2) + 
                    pars$surv_sl_b3 * (xb^3)) )
}

# update kernel functions
grow_sl_sd <- function(x,pars){
  pars$a_sl*(exp(pars$b_sl*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy_s <- function(x,y,pars){
  xb <- x_range_s(x, pars)
  # returns a *probability density distribution* for each x value
  # return( dnorm(y,  mean = pars$grow_sl_b0 + pars$grow_sl_b1*xb, 
  #                   sd   = pars$grow_sl_sig) )
  return( dnorm(y,  mean = pars$one_sz, 
                    sd   = pars$one_sd) )
}

# transition: Survival * growth
pxy_s<-function(x,y,pars){
  xb <- x_range_s(x, pars)
  return( sx_s(xb,pars) * gxy_s(xb,y,pars) )
}

x_range <- function(x,pars){
  pmin(pmax(x,pars$L),pars$U)
}

# Survival at size x
sx<-function(x,pars){
  xb <- x_range(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + 
                    pars$surv_b1 * xb +
                    pars$surv_b2 * xb^2 + 
                    pars$surv_b3 * xb^3 ) )
}

# update kernel functions
grow_sd <- function(x,pars){
  pars$a*(exp(pars$b*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy <- function(x,y,pars){
  xb <- x_range(x, pars)
  # returns a *probability density distribution* for each x value
  return( dnorm(y,  mean = pars$grow_b0 + pars$grow_b1*xb, 
                    sd   = pars$grow_sig) )
}

# transition: Survival * growth
pxy<-function(x,y,pars){
  xb <- x_range(x, pars)
  return( sx(xb,pars) * gxy(xb,y,pars) )
}


# production of seeds from x-sized mothers
fx <-function(x,pars){
  
  xb       <- x_range(x, pars)
  
  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + pars$flow_b1*xb ) * 
              exp(       pars$fert_b0 + pars$fert_b1*xb )
              
  # viable racs
  viab_rac <- tot_rac * ( 1-(pars$abort+pars$clip) )
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars,h){
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd ) * h
}

# seeds to seedling direct transition
fxy_s <- function(x,y,pars,h){
  xb <- x_range_s(x, pars)
  fx(xb,pars) * pars$g0 * recs(y, pars, h)
}


# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(pars){
 
  # set up IPM domains --------------------------------------------------------
  
  # seedlings
  n   <- pars$mat_siz
  L_s <- pars$L_sl 
  U_s <- pars$U_sl
  #these are the upper and lower integration limits
  h_s <- (U_s-L_s)/n                   #Bin size
  b_s <- L_s+c(0:n)*h_s                #Lower boundaries of bins 
  y_s <- 0.5*(b_s[1:n]+b_s[2:(n+1)])   #Bins' midpoints

  # plants
  n   <- pars$mat_siz
  L   <- pars$L 
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # Kernel "pieces" --------------------------------------------------------
  
  # vectors
  plant_s1  <- numeric(n)
  plant_s2  <- numeric(n)
  rec_s1    <- numeric(n)
  rec_s2    <- numeric(n)
  s1_rec    <- numeric(n)
  s2_rec    <- numeric(n)
  s1_plant  <- numeric(n)
  s2_plant  <- numeric(n)
  
  # seeds 
  s_mat     <- matrix(0,2,2)  
   
  # populate kernel ------------------------------------------------------------
  
  # seeds that enter 1 yr-old seed bank 
  plant_s1  <- fx(y,pars) * (1 - pars$g0)
  # no seeds go directly to 2 yr-old seed bank!
  plant_s2  <- numeric(n)

  # Impossible: seedlings to seedling
  sl_sl1     <- matrix(0,n,n)
  
  # seeds that go directly to seedings germinate right away 
  pl_sl2     <- (outer(y,y_s, fxy_s, pars, h_s) ) %>% t 
  # pl_sl2     <- matrix(0,n,n)
  
  # recruits from the 1 yr-old seedbank 
  s1_rec[]   <- recs(y_s, pars, h_s) * pars$g1
  # s1_rec[]   <- rep(0,100)
  
  # seeds that enter 2 yr-old seed bank 
  s_mat[2,1] <- (1 - pars$g1)
  
  # recruits from the 2 yr-old seedbank 
  s2_rec[]   <- recs(y_s, pars, h_s) * pars$g2
  # s2_rec[]   <- rep(0,100)
  
  # seedling survival and growth into plants
  sl_pl3     <- (outer(y_s, y, pxy_s, pars)*h) %>% t
  # sl_pl3     <- matrix(0,n,n) 
  
  # survival and growth of adult plants
  pl_pl4     <- (outer(y,y,pxy,pars)*h) %>% t
  
  # utilities for check
  image_k <- function(x) t(apply(x, 2, rev)) %>% image
  # sl_pl3 %>% image_k
  # pl_pl4 %>% image_k
  
  # Assemble the kernel -------------------------------------------------------------
  
  # "big blocks"
  big_mat    <- rbind( cbind(sl_sl1,pl_sl2),
                       cbind(sl_pl3,pl_pl4) )
  
  # top 2 vectors 
  top_vec    <- rbind( c(rec_s1, plant_s1),
                       c(rec_s2, plant_s2) )

  side_vec   <- rbind( s_mat,
                       cbind(s1_rec,   s2_rec  ),
                       cbind(s1_plant, s2_plant) )
  
  
  k_yx       <- cbind( side_vec, 
                       rbind(top_vec, big_mat) 
                       )
     
  return(k_yx)
   
  # tests "integrating' functions ---------------------------------------------------
  
  # s_sl
  # expect_true( ((outer( rep(1,100), y_s, s_sl, pars, h_s) %>% t %>% colSums) > 0.99) %>% all )

  # gxy 
  # expect_true( ((outer(y,y,gxy,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
  # gxy_s: huge unintentional eviction. Why?
  # expect_true( ((outer(y_s,y,gxy_s,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
}

ker <- kernel(pars_mean)
Re(eigen(ker)$value[1])
