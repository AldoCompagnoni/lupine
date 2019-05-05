# Hyper-simplified IPM for debugging
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(readxl)
library(testthat)
source('analysis/vital_rates/plot_binned_prop.R')
# "AL (1)"   "ATT (8)"  "BR (6)"   "BS (7)"   "DR (3)"   "NB (2)"   "POP9 (9)"
# "AL (1)" is good
# "ATT (8)" is not


year_trans  <- c(2009:2017)
all_indiv_sample <- c("BS (7)", 'DR (3)', 'NB (2)')
site_id     <- 'BS (7)'

# data
lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                # analyze a specific site only
                subset( location == site_id)
hist(lupine_df$log_area_t0, freq=F)
abline(v=1)
         
lupine_df   <- subset(lupine_df, log_area_t0 > 1)
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
germ        <- read_xlsx('data/seedbaskets.xlsx')
# sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
sl_size     <- data.frame( mean_sl_size = 2.725375531,
                           sd_sl_size   = 0.914582829, 
                           max_sl_size	 = 6.082794487,
                           min_sl_size  = -0.241564475 )
			

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
                  log_area_t02 = log(area_t0)^2 )


fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
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
mod_s    <- glmer(surv_t1 ~ log_area_t0 + 
                            log_area_t02 + 
                            (1|year), data=surv, family='binomial')
mod_g    <- lmer( log_area_t1 ~ log_area_t0 + 
                            (1|year), data=grow) 
g_lim    <- range(lupine_df$log_area_t0,na.rm=T)
mod_fl   <- glmer(flow_t0 ~ log_area_t0 + 
                            (1|year), data=flow, family='binomial')
# abort    <- glm(flow_t0 ~ log_area_t0, data=fert, family='binomial')
# clip     <- glm(flow_t0 ~ log_area_t0, data=fert, family='binomial')
mod_fr   <- glmer(numrac_t0 ~ log_area_t0 +
                              (1|year), data=fert, family='poisson')
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
surv_p    <- coef(mod_s)$year
grow_p    <- coef(mod_g)$year
# grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- coef(mod_fl)$year
fert_p    <- coef(mod_fr)$year
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ_coef



# model validation plots ---------------------------------------------------

par( mfrow=c(2,2), mar=c(3,3,0.5,0.5), 
     mgp=c(1.8,0.7,0) )

# survival
plot_binned_prop(surv, 10, log_area_t0, surv_t1)
x_seq    <- seq(min(surv$log_area_t0),
                max(surv$log_area_t0),by=0.1)
y_pred   <- boot::inv.logit( fixef(mod_s)[1] + 
                             fixef(mod_s)[2]*x_seq + 
                             fixef(mod_s)[3]*(x_seq^2) )
lines(x_seq, y_pred)

# growth
plot(log_area_t1 ~ log_area_t0, data=grow)
x_seq    <- seq(min(grow$log_area_t0),
                max(grow$log_area_t0),by=0.1)
y_pred   <- fixef(mod_g)[1] + fixef(mod_g)[2]*x_seq
lines(x_seq, y_pred, lwd=2, col='red')
abline(0,1,lty=2)

# flowering
plot_binned_prop(flow, 10, log_area_t0, flow_t0)
x_seq    <- seq(min(flow$log_area_t0),
                max(flow$log_area_t0),by=0.1)
y_pred   <- boot::inv.logit( fixef(mod_fl)[1] + 
                             fixef(mod_fl)[2]*x_seq )
lines(x_seq, y_pred)

# fertility
plot(numrac_t0 ~ log_area_t0, data=fert)
x_seq    <- seq(min(fert$log_area_t0),
                max(fert$log_area_t0),by=0.1)
y_pred   <- exp( fixef(mod_fr)[1] + fixef(mod_fr)[2]*x_seq )
lines(x_seq, y_pred,lwd=3,col='red')


# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# x_range <- function(x,pars){
#   pmin(pmax(x,pars$L),pars$U)
# }

# Survival at size x
sx<-function(x,pars){
  # survival prob. of each x size class 
  inv_logit( pars$surv_b0 + 
             pars$surv_b1 * x + 
             pars$surv_b2 * (x^2) # + pars$surv_b3 * x^3
            ) 
}

# # update kernel functions
# grow_sd <- function(x,pars){
#   pars$a*(exp(pars$b*x)) %>% sqrt 
# }

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
  viab_rac <- tot_rac * (1-pars$abort) * (1-pars$clip)
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

# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }

pars_mean_ii <- function(ii){

  # list of mean IPM parameters. 
  list( # adults vital rates           
        surv_b0      = surv_p[ii,'(Intercept)'],
        surv_b1      = surv_p[ii,'log_area_t0'],
        surv_b2      = surv_p[ii,'log_area_t02'],
        #surv_b3      = surv_p['log_area_t03'],
       
       grow_b0      = grow_p[ii,'(Intercept)'],
       grow_b1      = grow_p[ii,'log_area_t0'],
       grow_sig     = summary(mod_g)$sigma,

       flow_b0      = flow_p[ii,'(Intercept)'],
       flow_b1      = flow_p[ii,'log_area_t0'],
       
       fert_b0      = fert_p[ii,'(Intercept)'],
       fert_b1      = fert_p[ii,'log_area_t0'],
       
       abort        = 0.22, # hardcoded for now!
       clip         = 0.57, # hardcoded for now!
       
       fruit_rac    = fr_rac_p,
       seed_fruit   = seed_fr_p,
       g0           = 0.05,#germ_coef['g0'],
       g1           = germ_coef['g1'],
       g2           = germ_coef['g2'],
       
       recr_sz      = size_sl_p$mean_sl_size,
       recr_sd      = size_sl_p$sd_sl_size,
       
       L       = g_lim[1], #-0.2415645,
       U       = g_lim[2], #9.3550582, 
       
       mat_siz    = 200 )

}


# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(pars){
 
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
  
  # # seeds mini matrix 
  # s_mat     <- matrix(0,2,2)  
  #  
  # # seeds that enter 1 yr-old seed bank 
  # plant_s1  <- fx(y,pars) * (1 - pars$g0)
  # 
  # # no seeds go directly to 2 yr-old seed bank!
  # plant_s2  <- numeric(n)

  # seeds that go directly to seedlings germinate right away 
  Fmat       <- (outer(y,y, fxy, pars) * pars$g0 * h) 

  # # recruits from the 1 yr-old seedbank 
  # s1_rec     <- h * recs(y, pars) * pars$g1
  
  # # seeds that enter 2 yr-old seed bank 
  # s_mat[2,1] <- (1 - pars$g1)
  
  # # recruits from the 2 yr-old seedbank 
  # s2_rec     <- recs(y, pars) * pars$g2 * h
  
  # survival and growth of adult plants
  Tmat       <- (outer(y,y,pxy,pars)*h) 
  
  # rotate <- function(x) t(apply(x, 2, rev))
  # outer(y,y, fxy, pars, h) %>% t %>% rotate %>% image
  
  small_K    <- Tmat + Fmat
  
  # Assemble the kernel -------------------------------------------------------------
  
  # # top 2 vectors 
  # from_plant <- rbind( rbind( plant_s1, plant_s2),
  #                      small_K )
  # 
  # from_seed  <- rbind( s_mat,
  #                      cbind(s1_rec, s2_rec) )
  # 
  
  # k_yx       <- cbind( from_seed, from_plant )
     
  # return(k_yx)
   
  return(small_K)
  
  # tests "integrating' functions ---------------------------------------------------
  
  # gxy 
  # expect_true( ((outer(y,y,gxy,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
}



for(ii in 1:9){
  
  pars_mean <- pars_mean_ii(ii)
    
  ker <- kernel(pars_mean)
  Re(eigen(ker)$value[1])
  
  eK=eigen(ker); lambda=Re(eK$values[1]);
  w=Re(eK$vectors[,1]); w=w/sum(w);
  v=Re(eigen(t(ker))$vectors[,1]); v=v/v[1];
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
  
  # realized lambda
  lam_r_s[ii] <- sum(ker %*% (freq_v))/sum(freq_v)
  
  pop_n   <- read.csv( "data/lupine_all.csv") %>% 
                  subset( location == site_id ) %>% 
                  subset( !is.na(stage_t0) ) %>% 
                  subset( log_area_t0 > 1 ) %>% 
                  count( year, location ) 
                  
    
  lam_r   <- subset(pop_n, year == year_trans[ii]+1)$n /
             subset(pop_n, year == year_trans[ii])$n
  
  par(mar=c(3.5,3.5,0.5,0.5))
  plot(1:10,type='n')
  legend('topright',
         as.character(round(Re(eigen(ker)$value[1]),3)),
          bty='n')
  legend('bottomleft',
         as.character( lam_r ),
         bty='n')
  
  lam_real[ii] <- lam_r
  lam_det[ii]  <- Re(eigen(ker)$value[1])
  
}


par( mfrow=c(2,2), 
     mar=c(2.5,2.5,0.2,0.2),
     mgp=c(1.8,0.7,0) )
plot(lam_real, lam_det)
abline(0,1)
plot(lam_real, lam_r_s)
abline(0,1)
plot(lam_det, lam_r_s)
abline(0,1)


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
