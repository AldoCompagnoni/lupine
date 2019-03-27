# Hyper-simplified IPM for debugging
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
source('analysis/vital_rates/plot_binned_prop.R')
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"


# create lambdas for every year
det_lambda <-function(germ_est, sb){
  
  # data
  fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
  seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
  germ        <- read_xlsx('data/seedbaskets.xlsx')
  sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
  
  # vital rates format --------------------------------------------------------------
  surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                    subset( area_t0 != 0) %>%
                    mutate( log_area_t0 = log(area_t0),
                            year        = year + 1 ) %>% 
                    mutate( log_area_t02 = log_area_t0^2,
                            log_area_t03 = log_area_t0^3) 
  
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
                    log_area_t02 = log(area_t0)^2,
                    year         = year + 1 )
  
  fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                    subset( area_t0 != 0) %>% 
                    subset( !is.na(numrac_t0) ) %>% 
                    # remove non-flowering indiv.
                    subset( !(flow_t0 %in% 0) ) %>% 
                    mutate( log_area_t0  = log(area_t0),
                            log_area_t02 = log(area_t0)^2,
                            year         = year + 1 ) %>% 
                    # remove zero fertility (becase fertility should not be 0)
                    # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                    subset( !(numrac_t0 %in% 0) )
  
  
  # models ---------------------------------------------------------
  mod_s    <- glm(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, data=surv, family='binomial')
  # mod_s2    <- glm(surv_t1 ~ log_area_t0 + log_area_t02, data=surv, family='binomial')
  mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow) 
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
  
  
  
  # model validation plots ---------------------------------------------------
  # tiff( paste0('results/validation/vr/',
  #              site_all$year[ii],'_',
  #              site_all$site_id[ii],'.tiff'),
  #       unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
  
  par( mfrow=c(2,2) )
  # survival
  plot_binned_prop(surv, 10, log_area_t0, surv_t1)
  x_seq    <- seq(min(surv$log_area_t0),
                  max(surv$log_area_t0),by=0.1)
  y_pred   <- boot::inv.logit( coef(mod_s)[1] + 
                               coef(mod_s)[2]*x_seq + 
                               coef(mod_s)[3]*(x_seq^2) +
                               coef(mod_s)[4]*(x_seq^3) )
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
  
  # dev.off()
  
  
  # IPM parameters -------------------------------------------------------------
  
  # function to extract values
  extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }
  
  
  # list of mean IPM parameters. 
  pars_mean   <- list( # adults vital rates           
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
                       g0           = ifelse(germ_est,0.035,
                                             germ_coef['g0']),
                       g1           = germ_coef['g1'],
                       g2           = germ_coef['g2'],
                       
                       recr_sz      = size_sl_p$mean_sl_size,
                       recr_sd      = size_sl_p$sd_sl_size,
                       
                       L       = g_lim[1], #-0.2415645,
                       U       = g_lim[2], #9.3550582, 
                       
                       mat_siz    = 200 )
  
  
  # IPM functions ------------------------------------------------------------------------------
  
  inv_logit <- function(x){ exp(x)/(1+exp(x)) } 
  
  # Survival at size x
  sx<-function(x,pars){
    # survival prob. of each x size class 
    inv_logit( pars$surv_b0 + 
               pars$surv_b1 * x + 
               pars$surv_b2 * (x^2) + 
               pars$surv_b3 * (x^3)
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
    plant_s1  <- fx(y,pars) * (1 - pars$g0)

    # no seeds go directly to 2 yr-old seed bank!
    plant_s2  <- numeric(n)

    # seeds that go directly to seedlings germinate right away 
    Fmat       <- (outer(y,y, fxy, pars) * pars$g0 * h) 
  
    # recruits from the 1 yr-old seedbank
    s1_rec     <- h * recs(y, pars) * pars$g1

    # seeds that enter 2 yr-old seed bank
    s_mat[2,1] <- (1 - pars$g1)

    # recruits from the 2 yr-old seedbank
    s2_rec     <- h * recs(y, pars) * pars$g2
    
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
     
    # return(small_K)
    
  }
  
  # simple IPM kernel/matrix  ------------------------------------------------------------
  kernel_s <- function(pars){
   
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
    
    # seeds that go directly to seedlings germinate right away 
    Fmat       <- (outer(y,y, fxy, pars) * pars$g0 * h) 
  
    # survival and growth of adult plants
    Tmat       <- (outer(y,y,pxy,pars)*h) 
    
    # rotate <- function(x) t(apply(x, 2, rev))
    # outer(y,y, fxy, pars, h) %>% t %>% rotate %>% image
    
    small_K    <- Tmat + Fmat
    
    return(small_K)
    
  }
  
  if(sb){
    ker <- kernel_sb(pars_mean)
  }else{
    ker <- kernel_s( pars_mean)
  }

  Re(eigen(ker)$value[1])
  
}

lupine_df   <- read.csv( "data/lupine_all.csv") #%>% 
                  # subset( !(location == 'AL (1)') )

det_lambda(germ_est = F, sb = T)
det_lambda(germ_est = F, sb = F)
det_lambda(germ_est = T, sb = T)
det_lambda(germ_est = T, sb = F)
