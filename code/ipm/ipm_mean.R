# IPM from data
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(testthat)

# vital rates
vr_list <- list.files('results/ml_mod_sel/') %>% 
              setdiff( grep('.tiff',list.files('results/ml_mod_sel/'), value=T)  ) %>% 
              setdiff( grep('.png',list.files('results/ml_mod_sel/'), value=T)  )

vr_file <- lapply(vr_list, function(x) 
                           paste0('results/ml_mod_sel/',x,
                                  '/',x,'_best_mod.csv') ) %>% 
              lapply( read.csv ) %>% 
              setNames( vr_list )

# vital rate models 
surv_sl_p <- vr_file$surv_sl
grow_sl_p <- vr_file$grow_sl
surv_p    <- vr_file$surv
grow_p    <- vr_file$grow
flow_p    <- vr_file$flow
fert_p    <- vr_file$fert
size_sl_p <- vr_file$size_sl


# IPM parameters -------------------------------------------------------------

# function to extract values
extr_value <- function(x, field){ subset(x, type_coef == 'fixef' & ranef == field )$V1 }


# list of mean IPM parameters. 
pars_mean   <- list( # seedlings vital rates
                     surv_sl_b0   = extr_value(surv_sl_p, "(Intercept)"),
                     surv_sl_b1   = extr_value(surv_sl_p, "log_area_t0"),
                     surv_sl_b1_2 = extr_value(surv_sl_p, "log_area_t02"),
                     surv_sl_clim = extr_value(surv_sl_p, "tmp_tm1"),
                    
                     grow_sl_b0   = extr_value(grow_sl_p, "(Intercept)"),
                     grow_sl_b1   = extr_value(grow_sl_p, "log_area_t0"),
                     a_sl         = extr_value(grow_sl_p, "a"),
                     b_sl         = extr_value(grow_sl_p, "b"),
                      
                     # adults vital rates           
                     surv_b0      = extr_value(surv_p, "(Intercept)"),
                     surv_b1      = extr_value(surv_p, "log_area_t0"),
                     surv_clim    = extr_value(surv_p, "tmp_tm1"),
                     
                     grow_b0      = extr_value(grow_p, "(Intercept)"),
                     grow_b1      = extr_value(grow_p, "log_area_t0"),
                     a            = extr_value(grow_p, "a"),
                     b            = extr_value(grow_p, "b"),
                     
                     flow_b0      = extr_value(flow_p, "(Intercept)"),
                     flow_b1      = extr_value(flow_p, "log_area_t1"),
                     flow_b1_2    = extr_value(flow_p, "log_area_t12"),
                     flow_clim    = extr_value(flow_p, "tmp_t0"),
                     
                     fert_b0      = extr_value(fert_p, "(Intercept)"),
                     fert_b1      = extr_value(fert_p, "log_area_t1"),
                     fert_b1_2    = extr_value(fert_p, "log_area_t12"),
                     fert_clim    = extr_value(fert_p, "tmp_t0"),
                     
                     abort        = 0.22, # hardcoded for now!
                     clip         = 0.57, # hardcoded for now!
                     
                     fruit_rac    = extr_value(fert_p, "fruit_per_raceme"),
                     seed_fruit   = extr_value(fert_p, "seed_per_fruit"),
                     g0           = extr_value(fert_p, "g0"),
                     g1           = extr_value(fert_p, "g1"),
                     g2           = extr_value(fert_p, "g2"),
                     
                     recr_sz      = size_sl_p$mean_sl_size,
                     recr_sd      = size_sl_p$sd_sl_size,
                     
                     L_sl         = size_sl_p$min_sl_size,
                     U_sl         = size_sl_p$max_sl_size,
                     
                     L       = extr_value(grow_p, "min_size"),
                     U       = extr_value(grow_p, "max_size"),
                     
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
                    pars$surv_sl_b1_2 * (xb^2)) )
}

# update kernel functions
grow_sl_sd <- function(x,pars){
  pars$a_sl*(exp(pars$b_sl*x)) %>% sqrt 
}

# growth (transition) from size x to size y
gxy_s <- function(x,y,pars){
  xb <- x_range_s(x, pars)
  # returns a *probability density distribution* for each x value
  return( dnorm(y,  mean = pars$grow_sl_b0 + pars$grow_sl_b1*xb, 
                    sd   = grow_sl_sd(xb,pars)) )
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
  return( inv_logit(pars$surv_b0 + pars$surv_b1 * xb) )
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
                    sd   = grow_sd(xb,pars)) )
}

# transition: Survival * growth
pxy<-function(x,y,pars){
  xb <- x_range(x, pars)
  return( sx(xb,pars) * gxy(xb,y,pars) )
}


# production of seeds from x-sized mothers
fx <-function(x,pars){
  
  xb        <- x_range(x, pars)
  
  # total racemes prod
  tot_rac  <- inv_logit( pars$flow_b0 + pars$flow_b1*xb + pars$flow_b1_2*(xb^2) ) * 
              exp(       pars$fert_b0 + pars$fert_b1*xb + pars$fert_b1_2*(xb^2) )
              
  # viable racs
  viab_rac <- tot_rac * (1-pars$abort) * (1-pars$clip)
  # viable seeds
  viab_sd  <- viab_rac * pars$fruit_rac * pars$seed_fruit
  return(viab_sd)
  
}

# Size distribution of recruits
recs <-function(y,pars,h){
  dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd ) * h
}

# seeds to seedling direct transition
s_sl <- function(seed_x,y,pars,h){
  seed_x * recs(y, pars, h)
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
  
  # blocks: 
  block_1   <- matrix(0,n,n)
  block_2   <- matrix(0,n,n)
  block_3   <- matrix(0,n,n)
  block_4   <- matrix(0,n,n)
  
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

  # seeds that go directly to seedings germinate right away 
  block_2    <- outer( (fx(y,pars) * pars$g0), y_s, s_sl, pars, h_s) %>% t 

  # recruits from the 1 yr-old seedbank 
  s1_rec[]   <- recs(y_s, pars, h_s) * pars$g1
  
  # seeds that enter 2 yr-old seed bank 
  s_mat[2,1] <- (1 - pars$g1)
  
  # recruits from the 2 yr-old seedbank 
  s2_rec[]   <- recs(y_s, pars, h_s) * pars$g2
  
  # seedling survival and growth into plants
  block_3[]  <- (outer(y_s,y,pxy_s,pars)*h) %>% t
  
  # survival and growth of adult plants
  block_4[]  <- (outer(y,y,pxy,pars)*h) %>% t
  
  # rotate <- function(x) t(apply(x, 2, rev))
  # block_3 %>% rotate %>% image
  # block_4 %>% rotate %>% image
  
  # Assemble the kernel -------------------------------------------------------------
  
  # "big blocks"
  big_mat    <- rbind( cbind(block_1,block_2),
                       cbind(block_3,block_4) )
  
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

