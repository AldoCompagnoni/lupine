# Hyper-simplified IPM for debugging
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
source("analysis/format_data/format_functions.R")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"
# "AL (1)" is good
# "ATT (8)" is not

# data
lupine_df   <- read.csv( "data/lupine_all.csv") #%>% 
                    #subset( (year > 2011 & year < 2016) & location == 'ATT (8)')
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")


# format climate data ----------------------------------------
years     <- c(2005:2018)
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
ppt_mat <- subset(clim, clim_var == "ppt") %>%
              prism_clim_form("precip", years, m_back, m_obs) %>% 
              year_anom('ppt')

tmp_mat <- subset(clim, clim_var == 'tmean') %>% 
              prism_clim_form('tmean', years, m_back, m_obs) %>% 
              year_anom('tmp')
  
enso_mat <- subset(enso, clim_var == 'oni' ) %>%
              month_clim_form('oni', years, m_back, m_obs) %>% 
              year_anom('oni')

# put together all climate
clim_mat <- Reduce( function(...) full_join(...),
                    list(ppt_mat,tmp_mat,enso_mat) )


# vital rates format --------------------------------------------------------------
surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                    subset( area_t0 != 0) %>%
                    mutate( log_area_t0 = log(area_t0),
                            year        = year + 1 ) %>% 
                    mutate( log_area_t02 = log_area_t0^2,
                            log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

seedl       <- subset(lupine_df, !is.na(surv_t1) ) %>% 
                  subset( year != 2005 ) %>% 
                  subset( stage_t0 == 'SL' ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0 = log(area_t0),
                          year        = year + 1 ) %>% 
                  mutate( log_area_t02= log_area_t0^2 ) %>% 
                  left_join( clim_mat ) 

grow        <- lupine_df %>% 
                  subset(!(stage_t0 %in% c("DORM", "NF", "SL")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>% 
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  left_join( clim_mat ) 

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  left_join( clim_mat )

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove 
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year + 1 ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) 


# models ---------------------------------------------------------
mod_s    <- glm(surv_t1 ~ log_area_t0 + tmp_tm1, data=surv, family='binomial')
mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range(c(grow$log_area_t0, grow$log_area_t1))
mod_fl   <- glm(flow_t0 ~ log_area_t0 + tmp_t0, data=flow, family='binomial')
abort    <- glm(flow_t0 ~ log_area_t0, data=fert, family='binomial')
clip     <- glm(flow_t0 ~ log_area_t0, data=fert, family='binomial')
mod_fr   <- glm(numrac_t0 ~ log_area_t0 + tmp_t0, data=fert, family='poisson')
# mod_fr   <- MASS::glm.nb(numrac_t1 ~ log_area_t1, data=fert )
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
surv_p    <- coef(mod_s)
grow_p    <- coef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- coef(mod_fl)
fert_p    <- coef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ_coef

mod_fr   <- glm(numrac_t1 ~ log_area_t1, data=fert, family='poisson')
# mod_fr   <- MASS::glm.nb(numrac_t1 ~ log_area_t1, data=fert )

plot(numrac_t1 ~ log_area_t1, data=fert)
x_seq    <- seq(min(fert$log_area_t1),
                max(fert$log_area_t1),by=0.1)
y_pred   <- exp( coef(mod_fr)[1] + coef(mod_fr)[2]*x_seq )
lines(x_seq, y_pred,lwd=3,col='red')


# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }


# list of mean IPM parameters. 
pars_mean   <- list( # adults vital rates           
                     surv_b0      = surv_p['(Intercept)'],
                     surv_b1      = surv_p['log_area_t0'],
                     surv_clim    = surv_p['tmp_tm1'],
                       
                     grow_b0      = grow_p['(Intercept)'],
                     grow_b1      = grow_p['log_area_t0'],
                     grow_sig     = grow_p[3],

                     flow_b0      = flow_p['(Intercept)'],
                     flow_b1      = flow_p['log_area_t1'],
                     flow_clim    = flow_p['tmp_t0'],
                     
                     fert_b0      = fert_p['(Intercept)'],
                     fert_b1      = fert_p['log_area_t1'],
                     fert_clim    = fert_p['tmp_t0'],
                     
                     abort        = 0.22, # hardcoded for now!
                     clip         = 0.57, # hardcoded for now!
                     
                     fruit_rac    = fr_rac_p,
                     seed_fruit   = seed_fr_p,
                     g0           = germ_coef['g0'],
                     g1           = germ_coef['g1'],
                     g2           = germ_coef['g2'],
                     
                     recr_sz      = size_sl_p$mean_sl_size,
                     recr_sd      = size_sl_p$sd_sl_size,
                     
                     L       = g_lim[1],
                     U       = g_lim[2],
                     
                     mat_siz_sl = 100,
                     mat_siz    = 100 )


# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# x_range <- function(x,pars){
#   pmin(pmax(x,pars$L),pars$U)
# }

# Survival at size x
sx<-function(x,pars,tmp_anom){
  # xb <- x_range(x, pars)
  xb <- x
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + 
                    pars$surv_b1 * xb + 
                    pars$surv_clim * tmp_anom) )
}

# update kernel functions
grow_sd <- function(x,pars){
  pars$a*(exp(pars$b*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy <- function(x,y,pars){
  # xb <- x_range(x, pars)
  xb <- x
  # returns a *probability density distribution* for each x value
  return( dnorm(y,  mean = pars$grow_b0 + pars$grow_b1*xb, 
                    sd   = pars$grow_sig) )
}

# transition: Survival * growth
pxy<-function(x,y,pars,tmp_anom){
  # xb <- x_range(x, pars)
  xb <- x
  return( sx(xb,pars,tmp_anom) * gxy(xb,y,pars) )
}


# production of seeds from x-sized mothers
fx <-function(x,pars,tmp_anom){
  
  # xb        <- x_range(x, pars)
  xb        <- x
  
  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + 
                         pars$flow_b1*xb +
                         pars$flow_clim*tmp_anom ) * 
              exp(       pars$fert_b0 + 
                         pars$fert_b1*xb + 
                         pars$fert_clim*tmp_anom )
              
  # viable racs
  viab_rac <- tot_rac * (1- (pars$abort+pars$clip) )
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars,h){
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd ) * h
}

fxy <- function(x,y,pars,h,tmp_anom){
  # xb        <- x_range(x, pars)
  xb        <- x
  fx(xb,pars,tmp_anom) * recs(y,pars,h) 
}


# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(tmp_anom,pars){
 
  # set up IPM domains --------------------------------------------------------
 
  # plants
  n   <- pars$mat_siz
  L   <- pars$L 
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # populate kernel ------------------------------------------------------------
  
  # seeds mini matrix 
  s_mat     <- matrix(0,2,2)  
   
  # seeds that enter 1 yr-old seed bank 
  plant_s1  <- fx(y,pars,tmp_anom) * (1 - pars$g0)
  # no seeds go directly to 2 yr-old seed bank!
  plant_s2  <- numeric(n)

  # seeds that go directly to seedings germinate right away 
  Fmat       <- (outer(y,y, fxy, pars, h, tmp_anom) * pars$g0) %>% t 

  # recruits from the 1 yr-old seedbank 
  s1_rec     <- recs(y, pars, h) * pars$g1
  
  # seeds that enter 2 yr-old seed bank 
  s_mat[2,1] <- (1 - pars$g1)
  
  # recruits from the 2 yr-old seedbank 
  s2_rec     <- recs(y, pars, h) * pars$g2
  
  # survival and growth of adult plants
  Tmat       <- (outer(y,y,pxy,pars,tmp_anom)*h) %>% t
  
  # rotate <- function(x) t(apply(x, 2, rev))
  # outer(y,y, fxy, pars, h) %>% t %>% rotate %>% image
  
  
  # Assemble the kernel -------------------------------------------------------------
  
  # top 2 vectors 
  from_plant <- rbind( rbind( plant_s1, plant_s2),
                       Tmat )

  from_seed  <- rbind( s_mat,
                       cbind(s1_rec, s2_rec) )
  
  
  k_yx       <- cbind( from_seed, from_plant )
     
  return(k_yx)
   
  # tests "integrating' functions ---------------------------------------------------
  
  # gxy 
  # expect_true( ((outer(y,y,gxy,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
}

ker <- kernel(0,pars_mean)
Re(eigen(ker)$value[1])

clim_x <- seq( min(clim_mat$tmp_t0),
               max(clim_mat$tmp_t0), 
               length.out = 100 )

clim_lambda <- function(xx,pars_mean){
  ker <- kernel(xx,pars_mean)
  Re(eigen(ker)$value[1])
}

lam_vs_clim <- sapply(clim_x, clim_lambda, pars_mean)
plot(clim_x, lam_vs_clim, type='l')
abline(h=1,lty=2,lwd=2)

df <- data.frame( x=clim_x, y=lam_vs_clim)

ggplot(df, aes( x, y) ) +
  geom_line( size = 1 ) + 
  geom_hline( yintercept = 1,
              linetype   = 2 ) + 
  ylab( expression(lambda) ) + 
  xlab( 'Climate anomaly' ) +
  theme( axis.title = element_text( size = 30) ) + 
  ggsave( 'results/lambda_vs_clim.tiff',
          width = 6.3, height = 6.3, compression="lzw" )

  