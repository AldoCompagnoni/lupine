# Compare Dangremond's results with IPM
# 2. Calculate deterministic lambda for IPM/MPMs 
rm(list=ls())
source("analysis/format_data/format_functions.R")
source('analysis/vital_rates/plot_binned_prop.R')
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(readxl)
library(testthat)
library(lme4)
# Sites are: AL (1) (Abbotts Lagoon) and ATT (8) (Radio Tower)
# We have no "NB (2)" data before 2008

# Calculate deterministic lambda for IPM/MPMs -------------------------------------

# Before-2008 site-specific info
al  <- data.frame( year    = c(2005:2007),
                   fr      = c(3.025, 4.599, 1.286),
                   sf       = 3.563,
                   cons     = 1 - c(0.300, 0.257, 0.180),
                   abor     = 0) %>% 
          right_join( data.frame( year    = c(2005:2017),
                                  site_id = 'AL (1)') )
at  <- data.frame( year    = c(2005:2007),
                   fr      = c(2.837, 5.198, 5.872),
                   sf      = 3.563,
                   cons    = 1 - c(0.495, 0.667, 0.586),
                   abor    = 0) %>% 
          right_join( data.frame( year    = c(2005:2017),
                                  site_id = 'ATT (8)') )

# Post 2008 rates
fruit_rac <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr <- read_xlsx('data/seedsperfruit.xlsx')
fr_post08 <- glm(NumFruits ~ 1, 
                data=fruit_rac, 
                family='poisson') %>% coef %>% exp()
sf_post08 <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log")) %>% 
                coef %>% exp

# all site
site_all <- list(al, at) %>% 
              bind_rows %>% 
              mutate( cons = replace(cons,
                                     is.na(cons),
                                     0.57) ) %>%
              mutate( abor = replace(abor,
                                     is.na(abor),
                                     0.22) ) %>%
              mutate( fr   = replace(fr,
                                     is.na(fr),
                                     fr_post08)) %>% 
              mutate( sf   = replace(sf,
                                     is.na(sf),
                                     sf_post08))


# create lambdas for every year
yr_lambdas <-function(ii, germ_est=F, dangrem_sb, fixed_fr, dan_g){
  
  # data
  lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                  subset( year == site_all$year[ii] ) %>% 
                  subset( location == site_all$site_id[ii] )
                  
  # lupine_df   <- subset(lupine_df, log_area_t0 > 1)
  germ        <- read_xlsx('data/seedbaskets.xlsx')
  # sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
  sl_size     <- data.frame( mean_sl_size = 2.725375531,
                             sd_sl_size   = 0.914582829, 
                             max_sl_size	= 6.082794487,
                             min_sl_size  = -0.241564475 )
  
  
  # vital rates format --------------------------------------------------------------
  surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                    subset( area_t0 != 0) %>%
                    mutate( log_area_t0 = log(area_t0) ) %>% 
                    mutate( log_area_t02 = log_area_t0^2,
                            log_area_t03 = log_area_t0^3 )
  
  grow        <- lupine_df %>% 
                    subset(!(stage_t0 %in% c("DORM", "NF") ) & 
                           !(stage_t1 %in% c("D", "NF", "DORM")) ) %>% 
                    # remove zeroes from area_t0 and area_t1
                    subset( area_t0 != 0) %>%
                    subset( area_t1 != 0) %>% 
                    mutate( log_area_t1  = log(area_t1),
                            log_area_t0  = log(area_t0),
                            log_area_t02 = log(area_t0)^2 )
  
  flow <- subset(lupine_df, !is.na(flow_t0) ) %>% 
            subset( area_t0 != 0) %>% 
            mutate( log_area_t0  = log(area_t0),
                    log_area_t02 = log(area_t0)^2 )
  
  fert <- subset(lupine_df, flow_t0 == 1 ) %>% 
                    subset( area_t0 != 0) %>% 
                    subset( !is.na(numrac_t0) ) %>% 
                    # remove non-flowering indiv.
                    subset( !(flow_t0 %in% 0) ) %>% 
                    mutate( log_area_t0  = log(area_t0),
                            log_area_t02 = log(area_t0)^2 ) %>% 
                    # remove zero fertility (becase fertility should not be 0)
                    # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                    subset( !(numrac_t0 %in% 0) )
  
  # models ---------------------------------------------------------
  
  # survival: quadratic predictor or not?
  mod_s1   <- glm(surv_t1 ~ log_area_t0 + log_area_t02, data=surv, family='binomial')
  mod_s2   <- glm(surv_t1 ~ log_area_t0, data=surv, family='binomial')
  mod_l    <- list(mod_s1, mod_s2)
  mod_sel  <- c(AIC(mod_s1),AIC(mod_s2)) 
  mod_s    <- mod_l[[which(mod_sel == min(mod_sel))]]
  
  # other models
  mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow ) 
  g_lim    <- range(lupine_df$log_area_t0,na.rm=T)
  mod_fl   <- glm(flow_t0 ~ log_area_t0, data=flow, family='binomial')
  # abort    <- glm(flow_t0 ~ log_area_t0, data=fert, family='binomial')
  # clip     <- glm(flow_t0 ~ log_area_t0, data=fert, family='binomial')
  mod_fr   <- glm(numrac_t0 ~ log_area_t0, data=fert, family='poisson')
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
  germ_coef['g0'] <- 0.035 
  
  # models parameters -------------------------------------------------
  surv_p    <- coef(mod_s)
  grow_p    <- coef(mod_g) 
  grow_p    <- c(grow_p, summary(mod_g)$sigma)
  flow_p    <- coef(mod_fl)
  fert_p    <- coef(mod_fr)
  size_sl_p <- sl_size
  fr_rac_p  <- coef(fr_rac) %>% exp
  seed_fr_p <- coef(seed_fr) %>% exp
  germ_p    <- germ_coef
  
  # use parameters from Dangremond et al. 2010
  if(dan_g){
    germ_p    <- setNames(c(0.011, 0.034, 0.015),
                          names(germ_coef) )  
  }else{
    germ_p    <- germ_coef
  }

  # if using year-specific dangremond parameters
  if(fixed_fr){
    fixed_p   <- subset(site_all, year == 2008 & site_id == 'AL (1)')
    clip_p    <- fixed_p$cons
    abor_p    <- fixed_p$abor
    fr_rac_p  <- fixed_p$fr
    seed_fr_p <- fixed_p$sf
  }else{
    damgre_p  <- subset(site_all, year    == site_all$year[ii] &
                                  site_id == site_all$site_id[ii] )
    clip_p    <- damgre_p$cons
    abor_p    <- damgre_p$abor
    fr_rac_p  <- damgre_p$fr
    seed_fr_p <- damgre_p$sf
  }
  
  # model validation plots ---------------------------------------------------
  tiff( paste0('results/ipm/validation_mpm/vr/',
               site_all$site_id[ii],'_',
               site_all$year[ii],'.tiff'),
        unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
  
  par( mfrow=c(2,2), mar=c(3,3,0.1,0.1),
       mgp = c(1.8,0.8,0))
  # survival
  plot_binned_prop(surv, 10, log_area_t0, surv_t1)
  coef_s   <- coef(mod_s)
  coef_s[3]<- ifelse(coef_s['log_area_t02'] %>% 
                      is.na,
                     0,
                     coef_s['log_area_t02'])
  x_seq    <- seq(min(surv$log_area_t0),
                  max(surv$log_area_t0),by=0.1)
  y_pred   <- boot::inv.logit( coef_s[1] + 
                               coef_s[2]*x_seq + 
                               coef_s[3]*(x_seq^2) )
  lines(x_seq, y_pred)
  
  # growth
  plot(log_area_t1 ~ log_area_t0, data=grow)
  mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow )
  abline(mod_g)
  
  # flowering
  plot_binned_prop(flow, 10, log_area_t0, flow_t0)
  mod_fl   <- glm(flow_t0 ~ log_area_t0, data=flow, family='binomial')
  x_seq    <- seq(min(flow$log_area_t0),
                  max(flow$log_area_t0),by=0.1)
  y_pred   <- boot::inv.logit( coef(mod_fl)[1] + 
                               coef(mod_fl)[2]*x_seq )
  lines(x_seq, y_pred)
  
  # fertility
  plot(numrac_t0 ~ log_area_t0, data=fert)
  x_seq    <- seq(min(fert$log_area_t0),
                  max(fert$log_area_t0),by=0.1)
  y_pred   <- exp( coef(mod_fr)[1] + coef(mod_fr)[2]*x_seq )
  lines(x_seq, y_pred,lwd=3,col='red')
  
  dev.off()
  
  # IPM parameters -------------------------------------------------------------
  
  # list of mean IPM parameters. 
  pars_mean   <- list( # adults vital rates           
                       surv_b0      = surv_p['(Intercept)'],
                       surv_b1      = surv_p['log_area_t0'],
                       surv_b2      = ifelse(surv_p['log_area_t02'] %>% is.na,
                                             0,
                                             surv_p['log_area_t02']),
                       
                       grow_b0      = grow_p['(Intercept)'],
                       grow_b1      = grow_p['log_area_t0'],
                       grow_sig     = grow_p[3],
  
                       flow_b0      = flow_p['(Intercept)'],
                       flow_b1      = flow_p['log_area_t0'],
                       
                       fert_b0      = fert_p['(Intercept)'],
                       fert_b1      = fert_p['log_area_t0'],
                       
                       abort        = abor_p, 
                       clip         = clip_p, 
                       
                       fruit_rac    = fr_rac_p,
                       seed_fruit   = seed_fr_p,
                       g0           = ifelse(germ_est,0.035,
                                             germ_p['g0']),
                       g1           = germ_p['g1'],
                       g2           = germ_p['g2'],
                       
                       recr_sz      = size_sl_p$mean_sl_size,
                       recr_sd      = size_sl_p$sd_sl_size,
                       
                       L       = g_lim[1], #-0.2415645,
                       U       = g_lim[2], #9.3550582, 
                       
                       mat_siz    = 200 )
  
  # correct germination rates
  pars_mean$g1_c <- pars_mean$g1 / (1 - pars_mean$g0)
  pars_mean$g2_c <- pars_mean$g2 / (1 - (pars_mean$g0 + pars_mean$g1) ) 
  
  # IPM functions ------------------------------------------------------------------------------
  inv_logit <- function(x){ exp(x)/(1+exp(x)) } 
  
  # Survival at size x
  sx<-function(x,pars){
    # survival prob. of each x size class 
    inv_logit( pars$surv_b0 + 
               pars$surv_b1 * x + 
               pars$surv_b2 * (x^2) # + pars$surv_b3 * x^3
              ) 
  }
  
  # growth (transition) from size x to size y
  gxy <- function(y,x,pars){
    # returns a *probability density distribution* for each x value
    dnorm(y,  mean = pars$grow_b0 + pars$grow_b1*x, 
              sd   = pars$grow_sig)
  }
  
  # transition: Survival * growth
  pxy<-function(y,x,pars){
    return( sx(x,pars) * gxy(y,x,pars) )
  }
  
  # production of seeds from x-sized mothers
  fx <-function(x,pars){
    
    # total racemes prod
    tot_rac  <- inv_logit( pars$flow_b0 + pars$flow_b1*x ) * 
                exp(       pars$fert_b0 + pars$fert_b1*x )
                
    # viable racs
    viab_rac <- tot_rac * (1 - (pars$abort + pars$clip) )
    # viable seeds
    viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
    viab_sd
    
  }
  
  # Size distribution of recruits
  recs <-function(y,pars){
    dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd )
  }
  
  fxy <- function(y,x,pars){
    fx(x,pars) * recs(y,pars) 
  }
  
  # IPM kernel/matrix ------------------------------------------------------------
  kernel_sb <- function(pars){
   
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
    plant_s1   <- fx(y,pars) * (1 - pars$g0)
    
    # no seeds go directly to 2 yr-old seed bank!
    plant_s2   <- numeric(n)

    # seeds that go directly to seedlings germinate right away 
    Fmat       <- (outer(y,y, fxy, pars) * pars$g0 * h) 
  
    # recruits from the 1 yr-old seedbank
    s1_rec     <- h * recs(y, pars) * pars$g1_c 

    # seeds that enter 2 yr-old seed bank
    s_mat[2,1] <- (1 - pars$g1_c)

    # recruits from the 2 yr-old seedbank
    s2_rec     <- h * recs(y, pars) * pars$g2_c
    
    # survival and growth of adult plants
    Tmat       <- (outer(y,y,pxy,pars)*h) 
    
    # rotate <- function(x) t(apply(x, 2, rev))
    # outer(y,y, fxy, pars, h) %>% t %>% rotate %>% image
    
    small_K    <- Tmat + Fmat
    
    # Assemble the kernel -------------------------------------------------------------
    
    # top 2 vectors
    from_plant <- rbind( rbind( plant_s1, plant_s2),
                         small_K )

    # leftmost vectors
    from_seed  <- rbind( s_mat,
                         cbind(s1_rec, s2_rec) )

    k_yx       <- cbind( from_seed, from_plant )

    return(k_yx)
     
  }
  
  # kernel dangremond
  kernel_dangre <- function(pars){
   
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

    # seeds that enter 2 yr-old seed bank
    plant_s2   <- fx(y,pars) * pars$g2
    
    # seeds that enter 1 yr-old seed bank
    plant_s1   <- fx(y,pars) * pars$g1

    # seeds that go directly to seedlings germinate right away 
    Fmat       <- (outer(y,y, fxy, pars) * pars$g0 * h) 
  
    # seeds that enter 2 yr-old seed bank
    s_mat[2,1] <- 1

    # recruits from the 1 yr-old seedbank
    s1_rec     <- h * recs(y, pars) 
    
    # recruits from the 2 yr-old seedbank
    s2_rec     <- h * recs(y, pars)
    
    # survival and growth of adult plants
    Tmat       <- (outer(y,y,pxy,pars)*h) 
    
    # rotate <- function(x) t(apply(x, 2, rev))
    # outer(y,y, fxy, pars, h) %>% t %>% rotate %>% image
    
    small_K    <- Tmat + Fmat
    
    # Assemble the kernel -------------------------------------------------------------
    
    # top 2 vectors
    from_plant <- rbind( rbind( plant_s2, plant_s1),
                         small_K )

    # leftmost vectors
    from_seed  <- rbind( s_mat,
                         cbind(s1_rec, s2_rec) )

    k_yx       <- cbind( from_seed, from_plant )

    return(k_yx)
     
  }
  
  if(dangrem_sb){
    ker <- kernel_dangre(pars_mean)
  }else{
    ker <- kernel_sb(pars_mean)
  }

  ker
  
}

# Our IPM (lambdas w/ fixed fertility values, our germ. experiment)
al_ipm_l   <- lapply(1:13,  yr_lambdas, dangrem_sb=T, 
                            fixed_fr=T, dan_g=F)
att_ipm_l  <- lapply(14:26, yr_lambdas, dangrem_sb=T, 
                            fixed_fr=T, dan_g=F)

# IPM lambdas with flexible fertility, our germ. experiment
al_ipm_flex  <- lapply(1:3,   yr_lambdas,  germ_est=F, dangrem_sb=T, 
                              fixed_fr=F,  dan_g=F)
att_ipm_flex <- lapply(14:16, yr_lambdas,  germ_est=F, dangrem_sb=T, 
                              fixed_fr=F , dan_g=F)

# IPM lambdas with fixed fertility but with Dangremond germ params
al_ipm_dang  <- lapply(1:3,   yr_lambdas,  germ_est=F, dangrem_sb=T, 
                              fixed_fr=T,  dan_g=T)
att_ipm_dang <- lapply(14:16, yr_lambdas,  germ_est=F, dangrem_sb=T, 
                              fixed_fr=T , dan_g=T)

# IPM lambdas w/ FIXED fertility values and Dangermond seedbank params
al_ipm_flex_dang <- lapply(1:3,   yr_lambdas,  germ_est=F, dangrem_sb=T, 
                                 fixed_fr=F,  dan_g=T)
att_ipm_flex_dang <- lapply(14:16, yr_lambdas,  germ_est=F, dangrem_sb=T, 
                                 fixed_fr=F,  dan_g=T)

# lambdas from dangremond et al. 
al_mpm_l  <- lapply(1:3, function(ii){ read_xlsx('data/dangremond/abbotts_ambient_mpm.xlsx',
                                                 col_names = F,
                                                 sheet=ii) %>% 
                                                  as.matrix })
att_mpm_l <- lapply(1:3, function(ii){ read_xlsx('data/dangremond/radiotower_ambient_mpm.xlsx',
                                                 col_names = F,
                                                 sheet=ii) %>% 
                                                  as.matrix })

# calculate 
calc_lam <- function(x){
  eig_i <- which( Re(eigen(x)$value) == max(Re(eigen(x)$value)) )
  eigen(x)$value[eig_i] %>% Re
}



tiff('results/ipm/validation_mpm/lams_ipm_vs_mpm.tiff',
     res=300, width=6.3,height=9,unit='in',compression='lzw')

par(mfrow = c(2,1), mar = c(1.8,3,1,0.1),
    mgp   = c(1.5,0.5,0) )
# Abbotts
plot(2005:2017,
     sapply(al_ipm_l, calc_lam),
     ylab = expression('Deterministic '*lambda),
     xlab='',
     type='l',ylim=c(0.4,1.65),
     col='grey',lwd=2,cex.lab=1.5,
     main = 'Abbotts')
lines(2005:2007,
      sapply(al_mpm_l,  calc_lam),
      lty=1,lwd=2)
points(2005:2007,sapply(al_mpm_l,  calc_lam), pch = 2)
lines(2005:2007, 
      sapply(al_ipm_flex,  calc_lam),
      lty=2, lwd=2, col = 'blue')
lines(2005:2007, 
      sapply(al_ipm_dang,  calc_lam),
      lty=4, lwd=2, col = 'red')
lines(2005:2007, 
      sapply(al_ipm_flex_dang,  calc_lam),
      lty=3, lwd=2, col = 'brown')
# radio tower
plot(2005:2017,
     sapply(att_ipm_l, calc_lam),
     ylab = expression('Deterministic '*lambda),
     xlab='',type='l',ylim=c(0.6,1.65),
     col='grey',lwd=2,cex.lab=1.5,
     main = 'Radio Tower')
lines(2005:2007,
      sapply(att_mpm_l,  calc_lam),
      lty=1,lwd=2)
points(2005:2007,
       sapply(att_mpm_l,  calc_lam), pch = 2)
lines(2005:2007, 
      sapply(att_ipm_flex,  calc_lam),
      lty=2, lwd=2, col = 'blue')
lines(2005:2007, 
      sapply(att_ipm_dang,  calc_lam),
      lty=4, lwd=2, col = 'red')
lines(2005:2007, 
      sapply(att_ipm_flex_dang,  calc_lam),
      lty=3, lwd=2, col = 'brown')
legend('bottomright',
       c('MPM by Dangremond',
         'IPM: our model',
         'IPM: Dangremond fert. + germ. params.',
         'IPM: Dangremond fert. params',
         'IPM: Dangremond germ. params.'),
       lty = c(1,1,3,2,4),    lwd=2,
       col = c('black','grey','brown','red','blue'),
       pch = c(2,NA,NA,NA,NA), bty='n')
dev.off()


# mpm stochastic simulation ------------------------------
it      <- 50000
set.seed(1776)
yr_seq  <- sample(1:3,it,replace=T)


# simulate stochastic lambda
lam_stoch_clim <- function(proj_l){

  dimension <- nrow(proj_l[[1]])
  n0        <- rep(1/dimension, dimension)
  r_v       <- rep(0,it)

  # stochastic simulations
  for(ii in 1:it){
    
    #Store kernel
    mat     <- proj_l[[yr_seq[ii]]]
    
    # calculate lambdas
    n0      <- mat %*% n0
    N       <- sum(n0)
    r_v[ii] <- log(N)
    n0      <- n0/N
  
  }

  # r_v[-c(1:1000)] %>% mean %>% exp
  r_v[] %>% mean %>% exp
  
}

lam_stoch_clim(lam_ipm_l)
lam_stoch_clim(abb_mpm_l)

