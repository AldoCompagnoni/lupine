# Hyper-simplified IPM for validation
# I compare observed lambda to deterministic lambda
# only for sites BS, NB, and DR
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)
library(bbmle)
source('analysis/vital_rates/plot_binned_prop.R')
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"

bs <- data.frame( site_id = 'BS (7)',
                  year    = c(2009:2017) )
dr <- data.frame( site_id = 'DR (3)',
                  year    = c(2009:2014,2016:2017) )
nb <- data.frame( site_id = 'NB (2)',
                  year    = c(2010,2012:2016) )

site_all <- list(bs, dr, nb) %>% bind_rows 


# create lambdas for every year
yr_lambdas <-function(ii, germ_est, sb){
# for(ii in 17:nrow(site_all)){
  
  # data
  lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                  subset( year == site_all$year[ii] ) %>% 
                  # subset( location == 'DR (3)')
                  subset( location == site_all$site_id[ii] )
                  # subset( location == 'NB (2)')
                  
  hist(lupine_df$log_area_t0, freq=F)
  abline(v=1)
           
  # lupine_df   <- subset(lupine_df, log_area_t0 > 1)
  
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
  
  # survival: quadratic predictor or not?
  mod_s1   <- glm(surv_t1 ~ log_area_t0 + log_area_t02, data=surv, family='binomial')
  mod_s2   <- glm(surv_t1 ~ log_area_t0, data=surv, family='binomial')
  mod_l    <- list(mod_s1, mod_s2)
  mod_sel  <- c(AIC(mod_s1),AIC(mod_s2)) 
  mod_s    <- mod_l[[which(mod_sel == min(mod_sel))]]
  
  # other models
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
  tiff( paste0('results/validation/vr/',
               site_all$year[ii],'_',
               site_all$site_id[ii],'.tiff'),
        unit="in", width=6.3, height=6.3, res=600,compression="lzw" )
  
  par( mfrow=c(2,2) )
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
  
  # function to extract values
  extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }
  
  
  # list of mean IPM parameters. 
  pars_mean   <- list( # adults vital rates           
                       surv_b0      = surv_p['(Intercept)'],
                       surv_b1      = surv_p['log_area_t0'],
                       surv_b2      = ifelse(surv_p['log_area_t02'] %>% 
                                               is.na,
                                             0,
                                             surv_p['log_area_t02']),
                       
                       #surv_b3      = surv_p['log_area_t03'],
                       
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
  
  eK=eigen(ker); lambda=Re(eK$values[1]);
  w=Re(eK$vectors[,1]); w=w/sum(w);
  v=Re(eigen(t(ker))$vectors[,1]); v=v/v[1];
  par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5))
  plot(w)
  
  # calculate h
  n   <- pars_mean$mat_siz
  L   <- pars_mean$L 
  U   <- pars_mean$U
  h   <- (U-L)/n                  
  
  # set up observed frequencies of 
  freq_v <- lupine_df$log_area_t0 %>% 
              cut( breaks=seq(g_lim[1]-0.00001, 
                              g_lim[2]+0.00001, 
                              length.out=(pars_mean$mat_siz+1)) ) %>% 
              table %>% 
              as.vector 
  
  # realized lambda: simulated
  # lam_r_s <- sum(ker %*% (freq_v))/sum(freq_v)
  lam_r_s <- NA
  
  pop_n   <- read.csv( "data/lupine_all.csv") %>% 
                  subset( location == site_all$site_id[ii] ) %>% 
                  subset( !is.na(stage_t0) ) %>% 
                  subset( log_area_t0 > 1 ) %>% 
                  count( year, location ) 
                  
  # realized lambda: observed
  lam_r   <- subset(pop_n, year == site_all$year[ii]+1)$n /
             subset(pop_n, year == site_all$year[ii])$n
  
  legend('topright',
         as.character(round(Re(eigen(ker)$value[1]),3)),
          bty='n')
  legend('bottomleft',
         as.character( lam_r ),
         bty='n')
  
  print(paste0('end of ', ii))
  
  list( lam_det = Re(eigen(ker)$value[1]),
        lam_r   = lam_r,
        lam_r_s = lam_r_s,
        germ_est= germ_est,
        sb      = sb )

}

# F,F
# F,T
# T,T
# T,F

# store lambdas
lam_df   <- lapply(1:nrow(site_all), yr_lambdas, 
                   germ_est=T, sb=T) %>% bind_rows

# potential output labes
label_df <- data.frame( title = c('Seedbank model, experimental data',
                                  'No seedbank model, experimental data',
                                  'Seedbank model, estimated data',
                                  'No seedbank model, estimated data'),
                        file  = c( 'lambda_seedbank_experimental_g',
                                   'lambda_Noseedbank_experimental_g',
                                   'lambda_seedbank_estimated_g',
                                   'lambda_Noseedbank_estimated_g'),
                        germ_est= c(F,F,T,T), 
                        sb    =   c(T,F,T,F),
                        stringsAsFactors = F
                      )
  
# store labels
plot_lab <- label_df %>% 
              subset( germ_est == first(lam_df$germ_est) &
                      sb       == first(lam_df$sb) )

# # store R squared  
# mod <- lm(lam_df$lam_det~lam_df$lam_r, offset=)
# R2  <- round(summary(mod)$r.squared, 3)

# plot results
ggplot(lam_df, aes(x=lam_r, y=lam_det)) +
  geom_point() +
  geom_abline( size=1.2 ) +
  xlab(expression(lambda['observed'])) +
  ylab(expression(lambda)) +
  # annotate( 'text', label = paste0('R2= ',R2), x=0.5, y=2) +
  theme( axis.title = element_text( size=25) ) + 
  ggtitle( plot_lab$title ) +
  ggsave( paste0('results/validation/',plot_lab$file,'.tiff'),
          width = 6.3, height = 6.3, compression="lzw" )


# Compare a validation in the literature
tenhumberg <- read.csv('C:/cloud/Dropbox/lupine/data/tenhumberg_validation.csv')

ggplot(tenhumberg, aes(lam_r, lam_s) ) +
  geom_point( ) + 
  xlab( expression(lambda['observed']) ) +
  ylab( expression(lambda['s']) ) +
  theme( axis.title = element_text( size=25) ) + 
  geom_abline( size = 1.2 ) +
  ggsave( 'results/validation/tenhumberg_validation.csv.tiff',
          width = 6.3, height = 6.3, compression="lzw" )

# nPar = length(pars_mean);
# sPar = numeric(nPar); # vector to hold parameter sensitivities
# dp = 0.01;            # perturbation for calculating sensitivities
# 
# for(j in 1:nPar){
#     m.par = pars_mean;
#     m.par[[j]]=m.par[[j]] - dp;
#     IPM.down = kernel(m.par);
#     lambda.down = Re(eigen(IPM.down)$values[1]);
#     m.par[[j]]=m.par[[j]] + 2*dp;
#     IPM.up = kernel(m.par);
#     lambda.up = Re(eigen(IPM.up)$values[1]);
#     sj = (lambda.up-lambda.down)/(2*dp);
#     sPar[j]=sj;
#     cat(j,names(pars_mean)[j],sj,"\n");
# }
# 
# graphics.off(); dev.new(width=11,height=6);
# par(mfrow=c(2,1),mar=c(4,2,2,1),mgp=c(2,1,0));
# barplot(sPar,names.arg=names(pars_mean),main="Parameter sensitivities");
# barplot(sPar*abs(pars_mean)/lambda,names.arg=names(pars_mean),main="Parameter elasticities");
