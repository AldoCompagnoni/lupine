# Compare deterministic lambda bw IPMs created:
# 1. Subsetting site- and year-specific data
# 2. Fitting a GLMM, and getting site- and year-specific params.
# 3. Plot lambdas
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

# Read and format data ------------------------------------------------------

# demographic data
lupine_df <- read.csv( "data/lupine_all.csv")

# Post 2008 rates
fruit_rac <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr <- read_xlsx('data/seedsperfruit.xlsx')
cons      <- read_xlsx('data/consumption.xlsx') %>% 
                mutate( Mean_consumption = Mean_consumption %>% as.numeric) %>% 
                select( Year, Site, Mean_consumption) %>% 
                # expand potential "cases"
                complete( Site, Year) %>% 
                # update name
                mutate( Site = toupper(Site) ) %>% 
                mutate( Mean_consumption = replace(Mean_consumption,
                                                   is.na(Mean_consumption),
                                                   mean(Mean_consumption,na.rm=T) 
                                                   ) )
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
pred_g    <- read_xlsx('data/post predation_lupinus tidestromii.xlsx')
abor_raw  <- subset(lupine_df, !is.na(flow_t0) & flow_t0 == 1 ) %>% 
                subset( !is.na(numrac_t0) ) %>% 
                # remove non-flowering individuals
                subset( !(flow_t0 %in% 0) ) %>% 
                # remove zero fertility (becase fertility should not be 0)
                subset( !(numrac_t0 %in% 0) ) %>% 
                # only years indicated by Tiffany
                subset( year %in% c(2010, 2011, 2013:2017) ) %>% 
                # calculate abortion rates
                mutate( ab_r      = numab_t0 / numrac_t0 ) %>% 
                group_by( location, year ) %>% 
                summarise( ab_r_m = mean(ab_r) ) %>% 
                right_join(
                  expand.grid( year     = 2005:2017,
                               location = unique(lupine_df$location),
                               stringsAsFactors = F ) 
                          ) 
# absurd workaround to get correct data in.
# UNDERSTAND WHAT IS GOING ON HERE.
abor    <- select(abor_raw, -ab_r_m)
ab_repl <- replace(abor_raw$ab_r_m, 
                   is.na(abor_raw$ab_r_m),
                   mean(abor_raw$ab_r_m, na.rm=T)) 
abor$ab_r_m <- ab_repl

# store all site- and year- combinations
site_all  <- expand.grid( year    = c(2008:2017),
                          site_id = unique(lupine_df$location),
                          stringsAsFactors = F ) %>% 
                subset( !(year == 2008 & site_id == 'BS (7)') ) %>%
                subset( !(year == 2008 & site_id == 'DR (3)') ) %>% 
                subset( !(year == 2008 & site_id == 'NB (2)') ) %>% 
                subset( !(year == 2009 & site_id == 'NB (2)') ) %>% 
                subset( !(year == 2015 & site_id == 'DR (3)') ) 


# test_co_ab <- function(ii){
#   cons_df     <- subset(cons, Year == site_all$year[ii] &
#                               Site == gsub(' \\([0-9]\\)','',site_all$site_id[ii]) )
#   abor_df     <- subset(abor, year     == site_all$year[ii] &
#                               location == site_all$site_id[ii])
#   cons_df$Mean_consumption + abor_df$ab_r_m
# }
# site_all_repr <- mutate( site_all, 
#                          tot_cons = sapply(1:nrow(site_all), 
#                                            test_co_ab) )

           
# create lambdas for every year
yr_lambdas <-function(ii){
  
  # data
  lupine_df   <- read.csv( "data/lupine_all.csv") %>% 
                   subset( year     == site_all$year[ii] ) %>% 
                   subset( location == site_all$site_id[ii] )
                  
  germ        <- read_xlsx('data/seedbaskets.xlsx')
  sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
  cons_df     <- subset(cons, Year == site_all$year[ii] &
                              Site == gsub(' \\([0-9]\\)','',site_all$site_id[ii]) )
  abor_df     <- subset(abor, year     == site_all$year[ii] &
                              location == site_all$site_id[ii])
  
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
  mod_s1   <- glm(surv_t1 ~ log_area_t0, data=surv, family='binomial')
  mod_s2   <- glm(surv_t1 ~ log_area_t0 + log_area_t02, data=surv, family='binomial')
  mod_s3   <- glm(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03, 
                  data=surv, family='binomial')
  mod_l    <- list(mod_s1, mod_s2, mod_s3)
  mod_sel  <- c(AIC(mod_s1),AIC(mod_s2),AIC(mod_s3)) 
  mod_s    <- mod_l[[which(mod_sel == min(mod_sel))]]
  if( ii == 51 ){
    mod_s  <- mod_s1
  }
  
  # other models
  mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow ) 
  g_lim    <- range(lupine_df$log_area_t0,na.rm=T)
  mod_fl   <- glm(flow_t0 ~ log_area_t0, data=flow, family='binomial')
  mod_fr   <- glm(numrac_t0 ~ log_area_t0, data=fert, family='poisson')
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
  germ_p    <- germ_coef * (1 - 0.43) # post-dispersal predation
  clip_p    <- cons_df$Mean_consumption
  abor_p    <- abor_df$ab_r_m
  
  # if clip + abor are > 1, make them both 0.5
  if( (clip_p + abor_p) > 1 ){
    clip_p <- 
    abor_p <- 0.5
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
  coef_s[3]<- ifelse(coef_s['log_area_t02'] %>% is.na,
                     0,
                     coef_s['log_area_t02'])
  coef_s[4]<- ifelse(coef_s['log_area_t03'] %>% is.na,
                     0,
                     coef_s['log_area_t03'])
  x_seq    <- seq(min(surv$log_area_t0),
                  max(surv$log_area_t0),by=0.1)
  y_pred   <- boot::inv.logit( coef_s[1] + 
                               coef_s[2]*x_seq + 
                               coef_s[3]*(x_seq^2) + 
                               coef_s[4]*(x_seq^3) )
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
                       surv_b3      = ifelse(surv_p['log_area_t03'] %>% is.na,
                                             0,
                                             surv_p['log_area_t03']),
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
                       g0           = germ_p['g0'],
                       g1           = germ_p['g1'],
                       g2           = germ_p['g2'],
                       
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
    boot::inv.logit( pars$surv_b0 + 
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
    tot_rac  <- boot::inv.logit( pars$flow_b0 + pars$flow_b1*x ) * 
                exp( pars$fert_b0 + pars$fert_b1*x )
                
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
  
  # kernel
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
  
  kernel(pars_mean)
  
}

# Our IPM (lambdas w/ fixed fertility values, our germ. experiment)
ker_sbst_l   <- lapply(1:65,  yr_lambdas)

# calculate deterministic lambda
calc_lam <- function(x){
  eig_i <- which( Re(eigen(x)$value) == max(Re(eigen(x)$value)) )
  eigen(x)$value[eig_i] %>% Re %>% unique
}

sapply(ker_sbst_l, function(x) is.na(x) %>% sum)
sapply(ker_sbst_l, calc_lam)

lam_df <- site_all %>% 
            mutate( det_lambda_subset = sapply(ker_sbst_l, 
                                               calc_lam) )

# 2. glmm lambdas ---------------------------------------------

# demographic data
lupine_df <- read.csv( "data/lupine_all.csv")

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

# models 
mod_s    <- glmer(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=surv, family='binomial')
mod_g    <- lmer( log_area_t1 ~ log_area_t0 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location),
                  data=grow)
mod_g2   <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range( c(grow$log_area_t0, grow$log_area_t1) )
mod_fl   <- glmer(flow_t0 ~ log_area_t0 + 
                  (1 | year) + (0 + log_area_t0 | year) +
                  (1 | location) + (0 + log_area_t0 | location),
                  data=flow, family='binomial')
mod_fr   <- glmer(numrac_t0 ~ log_area_t0 + 
                  (1 | year) + (0 + log_area_t0 | year) + 
                  (1 | location) + (0 + log_area_t0 | location), 
                  data=fert, family='poisson')
# mod_fr   <- MASS::glm.nb(numrac_t0 ~ log_area_t0 , data=fert )
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
surv_sl_p <- fixef(mod_sl)
grow_sl_p <- coef(mod_g_sl)
grow_sl_p <- c(grow_sl_p, summary(mod_g_sl)$sigma)
surv_p    <- fixef(mod_s)
grow_p    <- fixef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- fixef(mod_fl)
fert_p    <- fixef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ_coef * (1 - 0.43)


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
                     g0           = germ_p['g0'],
                     g1           = germ_p['g1'],
                     g2           = germ_p['g2'],
                     
                     recr_sz      = size_sl_p$mean_sl_size,
                     recr_sd      = size_sl_p$sd_sl_size,
                     
                     L_sl         = size_sl_p$min_sl_size,
                     U_sl         = size_sl_p$max_sl_size,
                     
                     L            = g_lim[1],
                     U            = g_lim[2],
                     
                     mat_siz_sl = 100,
                     mat_siz    = 100 )

# test that no NAs present
expect_equal(pars_mean %>% 
               unlist %>% 
               is.na() %>% 
               sum, 0)


# update with yearly parameters
update_par_spec <- function(year_n, loc_n, vr){
  
  pars_yr <- pars_mean
  
  # get year-specific parameters
  get_yr <- function(mod_obj,year_n){
  
    ranef(mod_obj)$year %>% 
      tibble::add_column(.,
                         .before=1,
                         year=row.names(.)) %>% 
      subset( year == year_n) %>% 
      select(-year)
  }
  
  # get location-specific parameters
  get_loc <- function(mod_obj,loc_n){
    ranef(mod_obj)$location %>% 
      tibble::add_column(.,
                         .before=1,
                         location=row.names(.)) %>% 
      subset( location == loc_n) %>% 
      select(-location)
  }
  
  # get info on raceme "loss"
  get_rac <- function(cons, abor, year_n, loc_n){
    
    # raceme consumption
    locc     <- gsub(' \\([0-9]\\)','',loc_n)
    out_cons <- cons %>% 
                  subset(Site == locc & Year == year_n ) %>% 
                  .$Mean_consumption %>% 
                  as.numeric
    # raceme abortion
    out_abor <- abor %>% 
                  subset( location == loc_n & year == year_n ) %>% 
                  .$ab_r_m
    
    if( (out_abor + out_cons) > 1 ){
      out_abor <-
      out_cons <- 0.5
    }
    
    data.frame( clip = out_cons,
                abor = out_abor)
    
  }
  
  if( vr == 'surv'){
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_yr(mod_s,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_yr(mod_s,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_loc(mod_s,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_loc(mod_s,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'grow'){
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_yr(mod_g,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_yr(mod_g,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_loc(mod_g,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_loc(mod_g,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'flow'){
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_yr(mod_fl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_yr(mod_fl,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_loc(mod_fl,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_loc(mod_fl,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'fert'){
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_yr(mod_fr,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_yr(mod_fr,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_loc(mod_fr,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_loc(mod_fr,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  if( vr == 'all' ){
    # years
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_yr(mod_s,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_yr(mod_s,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_yr(mod_g,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_yr(mod_g,year_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_yr(mod_fl,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_yr(mod_fl,year_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_yr(mod_fr,year_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_yr(mod_fr,year_n)['log_area_t0']) %>% as.numeric
    
    # locations
    pars_yr$surv_b0    <- (pars_yr$surv_b0 + get_loc(mod_s,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$surv_b1    <- (pars_yr$surv_b1 + get_loc(mod_s,loc_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$grow_b0    <- (pars_yr$grow_b0 + get_loc(mod_g,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$grow_b1    <- (pars_yr$grow_b1 + get_loc(mod_g,loc_n)['log_area_t0']) %>% as.numeric
    
    pars_yr$flow_b0    <- (pars_yr$flow_b0 + get_loc(mod_fl,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$flow_b1    <- (pars_yr$flow_b1 + get_loc(mod_fl,loc_n)['log_area_t0']) %>% as.numeric
  
    pars_yr$fert_b0    <- (pars_yr$fert_b0 + get_loc(mod_fr,loc_n)['(Intercept)']) %>% as.numeric
    pars_yr$fert_b1    <- (pars_yr$fert_b1 + get_loc(mod_fr,loc_n)['log_area_t0']) %>% as.numeric
  }
  
  # update data on clipped racemes
  pars_yr$clip       <- get_rac(cons, abor, year_n, loc_n)$clip
  pars_yr$abort      <- get_rac(cons, abor, year_n, loc_n)$abor
  
  pars_yr
  
}

# test function
update_par_spec(2006,"BS (7)",'all')
update_par_spec(2010,"POP9 (9)",'all')$abor
update_par_spec(2010,"POP9 (9)",'all')$clip


# IPM functions ------------------------------------------------------------------------------

# Survival at size x
sx<-function(x,pars){
  # survival prob. of each x size class 
  boot::inv.logit(pars$surv_b0 + 
                  pars$surv_b1 * x + 
                  pars$surv_b2 * x^2 +
                  pars$surv_b3 * x^3) 
}

# update kernel functions
grow_sd <- function(x,pars){
  pars$a*(exp(pars$b*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy <- function(y,x,pars){
  # returns a *probability density distribution* for each x value
  dnorm(y,  mean = pars$grow_b0 + pars$grow_b1*x, 
            sd   = pars$grow_sig)
}

# transition: Survival * growth
pxy<-function(y,x,pars){
  sx(x,pars) * gxy(y,x,pars) 
}

# production of seeds from x-sized mothers
fx <-function(x,pars){

  # total racemes prod
  tot_rac  <- boot::inv.logit( pars$flow_b0 + 
                               pars$flow_b1*x ) * 
              exp(       pars$fert_b0 + 
                         pars$fert_b1*x )
              
  # viable racs
  viab_rac <- tot_rac * (1- (pars$abort+pars$clip) )
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars){
  dnorm(y, mean = pars$recr_sz, 
           sd   = pars$recr_sd )
}

fxy <- function(y,x,pars){
  fx(x,pars) * recs(y,pars) 
}


# IPM kernel/matrix ------------------------------------------------------------
kernel_glmm <- function(pars){
 
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
  Tmat       <- outer(y,y,pxy,pars) * h
  
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

  # tests "integrating' functions ---------------------------------------------------
  
  # gxy 
  # expect_true( ( (outer(y,y,gxy,pars)*  h) %>% colSums > 0.97) %>% all)
  
  # only return if tests are true
  return(k_yx)
  
}


# calculate stochastic lambda for each climate anomaly
ker_glmm <- function(ii){
  
  # use 'site_all'
  # yr- and site-specific glmm parametrs
  par_yr_loc <- update_par_spec(site_all$year[ii],
                                site_all$site_id[ii],
                                'all')
  
  kernel_glmm(par_yr_loc)
  
}

ker_glmm_l <- lapply(1:nrow(site_all), ker_glmm)

# update data frame
lam_df <- lam_df %>% 
            mutate( det_lambda_glmm   = sapply(ker_glmm_l, 
                                               calc_lam) ) %>% 
            complete( year, site_id ) 
  


# 3. Plots lambdas -------------------------------------------------------

# format for plotting
lam_df <- lam_df %>% 
            # gather( type, lambda, 
            #         det_lambda_subset, 
            #         det_lambda_glmm ) %>% 
            mutate( year = as.numeric(as.character(year)) )
            


# expand

ggplot(lam_df) +
  geom_line( aes(x=year, 
                 y=lambda,
                 color= site_id),
             size = 1,
             alpha= 0.8) +
  scale_x_continuous(breaks = 0:2100) +
  ylab( expression(lambda) ) +
  scale_color_viridis_d() +
  theme( axis.title.y = element_text( size=30),
         axis.text.x  = element_text( angle=60) ) +
  facet_grid( 1 ~ type) +
  ggsave('results/ipm/validation/subset_vs_glmm/glmm_vs_subset.tiff',
         width=6.3, height=3, compression='lzw')

# 
ggplot(lam_df) +
  geom_line( aes(x=year, 
                 y=lambda,
                 color    = site_id,
                 linetype = type),
             size  = 1.5) +
  scale_x_continuous(breaks = 0:2100) +
  ylab( expression(lambda) ) +
  scale_color_viridis_d() +
  theme( axis.title.y = element_text( size=30),
         axis.text.x  = element_text( angle=60) ) +
  facet_grid( site_id ~ 1 ) +
  ggsave('results/ipm/validation/subset_vs_glmm/glmm_vs_subset_bysite.tiff',
         width=6.3, height=9, compression='lzw')
