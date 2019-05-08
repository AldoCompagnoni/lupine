# Calculate sensitivities of the mean model
rm(list=ls())
source("analysis/format_data/format_functions.R")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(testthat)

# data
lupine_df   <- read.csv( "data/lupine_all.csv") 
fruit_rac   <- read_xlsx('data/fruits_per_raceme.xlsx')
seed_x_fr   <- read_xlsx('data/seedsperfruit.xlsx')
pred_g      <- read_xlsx('data/post predation_lupinus tidestromii.xlsx')
sl_size     <- read.csv('results/ml_mod_sel/size_sl/seedl_size.csv')
clim        <- read.csv("data/prism_point_reyes_87_18.csv")
enso        <- read.csv("data/enso_data.csv")
germ        <- read_xlsx('data/seedbaskets.xlsx') %>% 
                 select(g0:g2) %>% 
                 colMeans
germ_adj    <- read.csv('results/ml_mod_sel/germ/germ_adj.csv')


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


# format climate data ----------------------------------------
years     <- c(2005:2018)
m_obs     <- 5
m_back    <- 36

# calculate yearly anomalies
year_anom <- function(x, var){
  
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

# first, format site/year combinations
site_df     <- select(lupine_df, year, location) %>% 
                 unique %>% 
                 # create 'Site' column to merge with consumption dat
                 mutate( Site = gsub(' \\([0-9]\\)','',location) ) %>% 
                 subset( year > 2004 ) %>% 
                 arrange( location, year ) %>% 
                 complete(location,year)

surv        <- subset(lupine_df, !is.na(surv_t1) ) %>%
                  subset( area_t0 != 0) %>%
                  mutate( log_area_t0 = log(area_t0),
                          year        = year ) %>% 
                  mutate( log_area_t02 = log_area_t0^2,
                          log_area_t03 = log_area_t0^3) %>% 
                  left_join( clim_mat ) 

grow        <- lupine_df %>% 
                  # remove sleedings at stage_t0
                  subset(!(stage_t0 %in% c("DORM", "NF")) & 
                         !(stage_t1 %in% c("D", "NF", "DORM")) ) %>%
                  # remove zeroes from area_t0 and area_t1
                  subset( area_t0 != 0) %>%
                  subset( area_t1 != 0) %>% 
                  mutate( log_area_t1  = log(area_t1),
                          log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year ) %>% 
                  left_join( clim_mat ) 

flow        <- subset(lupine_df, !is.na(flow_t0) ) %>% 
                  subset( area_t0 != 0) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year ) %>% 
                  left_join( clim_mat ) 

fert        <- subset(lupine_df, flow_t0 == 1 ) %>% 
                  subset( area_t0 != 0) %>% 
                  subset( !is.na(numrac_t0) ) %>% 
                  # remove non-flowering individuals
                  subset( !(flow_t0 %in% 0) ) %>% 
                  mutate( log_area_t0  = log(area_t0),
                          log_area_t02 = log(area_t0)^2,
                          year         = year ) %>% 
                  # remove zero fertility (becase fertility should not be 0)
                  # NOTE: in many cases, notab_t1 == 0, because numab_t1 == 0 also
                  subset( !(numrac_t0 %in% 0) ) %>% 
                  left_join( clim_mat ) 

abor_df  <- subset(lupine_df, !is.na(flow_t0) & flow_t0 == 1 ) %>% 
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
                summarise( ab_r_m = mean(ab_r, na.rm=T) ) %>% 
                ungroup %>% 
                right_join( select(site_df,-Site) ) %>% 
                mutate( ab_r_m = replace(ab_r_m,
                                         is.na(ab_r_m),
                                         mean(ab_r_m, 
                                              na.rm=T)) )

cons_df    <- read_xlsx('data/consumption.xlsx') %>% 
                mutate( Mean_consumption = Mean_consumption %>% as.numeric) %>% 
                select( Year, Site, Mean_consumption) %>% 
                # update contents/names of variables for merging
                mutate( Site = toupper(Site) ) %>% 
                rename( year = Year ) %>% 
                # expand to all site/year combinations
                right_join( site_df ) %>% 
                mutate( Mean_consumption = replace(Mean_consumption,
                                                   is.na(Mean_consumption),
                                                   mean(Mean_consumption,na.rm=T) 
                                                   ) ) %>% 
                # remove NA locations
                subset( !is.na(location) ) %>% 
                # remove annoying code
                select( -Site ) %>% 
                rename( cons = Mean_consumption ) %>% 
                arrange(location,year)

germ_df     <- site_df %>% 
                 select( location ) %>%
                 unique %>% 
                 arrange( location ) %>% 
                 left_join(germ_adj) %>% 
                 # create germ_obs for AL (1)
                 mutate( germ_obs = replace( germ_obs, 
                                             location == 'AL (1)',
                                             mean(germ_obs[location %in% c('ATT (8)',
                                                                           'POP9 (9)')])) ) %>% 
                 # create germ_obs for BR (6)
                 mutate( germ_obs = replace( germ_obs, 
                                             location == 'BR (6)',
                                             mean(germ_obs[location %in% c('BS (7)',
                                                                           'DR (3)')])) ) %>%   
                 # post-dispersal predation
                 mutate( post_d_p = (germ['g0'] - germ_obs) / germ['g0'] ) %>% 
                 # add germination rats
                 mutate( g0 = germ['g0'],
                         g1 = germ['g1'],
                         g2 = germ['g2'] )


# models ---------------------------------------------------------
mod_s    <- glm(surv_t1 ~ log_area_t0 + log_area_t02 + log_area_t03 + 
                  tmp_tm1, data=surv, family='binomial')
mod_g    <- lm( log_area_t1 ~ log_area_t0, data=grow) 
g_lim    <- range(c(grow$log_area_t0, grow$log_area_t1))
mod_fl   <- glm(flow_t0 ~ log_area_t0 + tmp_tm1, data=flow, family='binomial')
mod_fr   <- glm(numrac_t0 ~ log_area_t0 + tmp_tm1, data=fert, family='poisson')
fr_rac   <- glm(NumFruits ~ 1, data=fruit_rac, family='poisson')        
seed_fr  <- glm(SEEDSPERFRUIT ~ 1, 
                data=mutate(seed_x_fr,  
                            # substitute 0 value with really low value (0.01)
                            SEEDSPERFRUIT = replace(SEEDSPERFRUIT, 
                                                    SEEDSPERFRUIT == 0, 
                                                    0.01) ),
                family=Gamma(link = "log"))


# vital rate models 
surv_p    <- coef(mod_s)
grow_p    <- coef(mod_g) 
grow_p    <- c(grow_p, summary(mod_g)$sigma)
flow_p    <- coef(mod_fl)
fert_p    <- coef(mod_fr)
size_sl_p <- sl_size
fr_rac_p  <- coef(fr_rac) %>% exp
seed_fr_p <- coef(seed_fr) %>% exp
germ_p    <- germ * (1 - 0.9)
cons_p    <- subset(cons_df, location != 'AL (1)')$cons %>% mean(na.rm=T)
abor_p    <- abor_df$ab_r_m %>% mean

# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }


# list of mean IPM parameters. 
pars_mean   <- list( # adults vital rates           
                     surv_b0      = surv_p['(Intercept)'],
                     surv_b1      = surv_p['log_area_t0'],
                     surv_b2      = surv_p['log_area_t02'],
                     surv_b3      = surv_p['log_area_t03'],
                     surv_clim    = surv_p['tmp_tm1'],
                       
                     grow_b0      = grow_p['(Intercept)'],
                     grow_b1      = grow_p['log_area_t0'],
                     grow_sig     = grow_p[3],

                     flow_b0      = flow_p['(Intercept)'],
                     flow_b1      = flow_p['log_area_t0'],
                     flow_clim    = flow_p['tmp_tm1'],
                     
                     fert_b0      = fert_p['(Intercept)'],
                     fert_b1      = fert_p['log_area_t0'],
                     fert_clim    = fert_p['tmp_tm1'],
                     
                     abort        = abor_p,
                     clip         = cons_p,
                     
                     fruit_rac    = fr_rac_p,
                     seed_fruit   = seed_fr_p,
                     g0           = germ_p['g0'],
                     g1           = germ_p['g1'],
                     g2           = germ_p['g2'],
                     
                     recr_sz      = size_sl_p$mean_sl_size,
                     recr_sd      = size_sl_p$sd_sl_size,
                     
                     L       = g_lim[1],
                     U       = g_lim[2],
                     
                     mat_siz_sl = 100,
                     mat_siz    = 100 )

# test that no NAs present
expect_equal(pars_mean %>% 
               unlist %>% 
               is.na() %>% 
               sum, 0)


# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# Survival at size x
sx<-function(x,pars,tmp_anom){
  # survival prob. of each x size class 
  inv_logit(pars$surv_b0 + 
            pars$surv_b1 * x + 
            pars$surv_b2 * x^2 +
            pars$surv_b3 * x^3 +
            pars$surv_clim * tmp_anom) 
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
pxy<-function(y,x,pars,tmp_anom){
  sx(x,pars,tmp_anom) * gxy(y,x,pars) 
}


# production of seeds from x-sized mothers
fx <-function(x,pars,tmp_anom){

  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + 
                         pars$flow_b1*x +
                         pars$flow_clim*tmp_anom ) * 
              exp(       pars$fert_b0 + 
                         pars$fert_b1*x + 
                         pars$fert_clim*tmp_anom )
              
  # viable racs
  viab_rac <- tot_rac * (1- pars$abort) * (1-pars$clip)
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars){
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd )
}

fxy <- function(y,x,pars,tmp_anom){
  fx(x,pars,tmp_anom) * recs(y,pars) 
}


# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(tmp_anom, pars){
 
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
  plant_s2   <- fx(y,pars,tmp_anom) * pars$g2
  
  # seeds that enter 1 yr-old seed bank
  plant_s1   <- fx(y,pars,tmp_anom) * pars$g1

  # seeds that go directly to seedlings germinate right away 
  Fmat       <- (outer(y,y, fxy, pars, tmp_anom) * pars$g0 * h) 

  # seeds that enter 2 yr-old seed bank
  s_mat[2,1] <- 1

  # recruits from the 1 yr-old seedbank
  s1_rec     <- h * recs(y, pars) 
  
  # recruits from the 2 yr-old seedbank
  s2_rec     <- h * recs(y, pars)
  
  # survival and growth of adult plants
  Tmat       <- outer(y,y,pxy,pars,tmp_anom) * h
  
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
   
  # tests "integrating' functions ---------------------------------------------------
  
  # s_sl
  # expect_true( ((outer( rep(1,100), y_s, s_sl, pars, h_s) %>% t %>% colSums) > 0.99) %>% all )

  # gxy 
  # expect_true( ((outer(y,y,gxy,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
  # gxy_s: huge unintentional eviction. Why?
  # expect_true( ((outer(y_s,y,gxy_s,pars)*h) %>% t %>% colSums > 0.97) %>% all)
  
}

ker <- kernel(0,pars_mean)
lambda <- Re(eigen(ker)$value[1])

nPar <- length(pars_mean) - 4
sPar <- numeric(nPar)  # vector to hold parameter sensitivities
dp   <- 0.01           # perturbation for calculating sensitivities

for(j in 1:nPar){
    m.par       <- pars_mean
    m.par[[j]]  <- m.par[[j]] - dp
    IPM.down    <- kernel(0.5,m.par)
    lambda.down <- Re(eigen(IPM.down)$values[1])
    m.par[[j]]  <- m.par[[j]] + 2*dp
    IPM.up      <- kernel(0.5,m.par)
    lambda.up   <- Re(eigen(IPM.up)$values[1])
    sj          <- (lambda.up-lambda.down)/(2*dp)
    sPar[j]     <- sj
    cat(j,names(pars_mean)[j],sj,"\n");
}

graphics.off(); dev.new(width=11,height=6);
par(mfrow=c(2,1),mar=c(4,2,2,1),mgp=c(2,1,0));

# graph the sensitivity
sens_df <- data.frame( parameter   = names(pars_mean)[1:nPar],
                       sensitivity = sPar,
                       elasticity  = sPar*abs(as.numeric(pars_mean[1:nPar]))/lambda ) %>% 
            gather( measure, value, sensitivity:elasticity)

ggplot(sens_df) +
  geom_bar( aes( x = parameter,
                 y = value ),
                 stat = 'identity') +
  theme( axis.text.x = element_text( angle = 80,
                                     vjust = 0.5) ) +
  facet_grid( measure ~ 1 ) +
  ylab( 'Sensitivity/Elasticity') + 
  xlab( 'Parameter' ) +
  ggsave( 'results/ipm/ltre/sens_elast.tiff',
          width=6.3,height=6.3,compression='lzw')




