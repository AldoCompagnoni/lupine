# IPM from data
rm(list=ls())
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(testthat)

# vital rates
vr_list <- list.files('results/ml_mod_sel/') %>% 
              setdiff( grep('.tiff',list.files('results/ml_mod_sel/'), value=T)  )

vr_file <- lapply(vr_list, function(x) 
                           paste0('results/ml_mod_sel/',x,'/',x,'_best_mod.csv')) %>% 
            lapply( read.csv )
paste [1]

# vital rate models 
abor_p    <- read.csv( paste0(do_spp,"_grow.csv") )                  
clip_p    <- read.csv( paste0(do_spp,"_surv.csv") )                  
fert_p    <- read.csv( paste0(do_spp,"_recr.csv") )                  
flow_p    <- read.csv( paste0(do_spp,"_other.csv") )               
grow_p    <- read.csv( paste0(do_spp,"_other.csv") )               
surv_p    <- read.csv( paste0(do_spp,"_other.csv") )               
surv_sl_p <- read.csv( paste0(do_spp,"_other.csv") )               

                  
# columns for COMPADRE's metadata
spp_dscr  <- read_xlsx('C:/cloud/Dropbox/sApropos/data/Ulex_gallii_3.xlsx')
mat_stck  <- read_xlsx('C:/cloud/Dropbox/sApropos/data/Ulex_gallii_3.xlsx', sheet=2)
meta_col  <- grep('^((?!A[0-9]{1,3}).)*$',names(mat_stck), value=T, perl=T)


# IPM parameters -------------------------------------------------------------

extr_value <- function(x, field){ subset(x, coefficient == field)$value }

# list of mean IPM parameters. 
pars_mean   <- list( surv_b0 = extr_value(surv_p, "b0"),
                     surv_b1 = extr_value(surv_p, "logarea"),
                     grow_b0 = extr_value(grow_p, "b0"),
                     grow_b1 = extr_value(grow_p, "logarea.t0"),
                     a       = extr_value(grow_p, "a"),
                     b       = extr_value(grow_p, "b"),
                     fecu_b0 = extr_value(recr_p, "b0"),
                     recr_sz = extr_value(othr_p, "rec_siz"),
                     recr_sd = extr_value(othr_p, "rec_sd"),
                     L       = (extr_value(othr_p, "min_siz") * 1 ) %>% log,
                     U       = (extr_value(othr_p, "max_siz") * 1 ) %>% log,
                     mat_siz = 10 )


# IPM functions ------------------------------------------------------------------------------

inv_logit <- function(x){ exp(x)/(1+exp(x)) } 

# Transforms all values below/above limits in min/max size
x_range <- function(x,pars){
  pmin(pmax(x,pars$L),pars$U)
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

# Survival at size x
sx<-function(x,pars){
  xb <- x_range(x, pars)
  # survival prob. of each x size class 
  return( inv_logit(pars$surv_b0 + pars$surv_b1 * xb) )
}

# transition: Survival * growth
pxy<-function(x,y,pars){
  xb <- x_range(x, pars)
  return( sx(xb,pars) * gxy(xb,y,pars) )
}

# # Production of recruits from x-sized moms
# fx<-function(x,pars){
#   xb      <- x_range(x, pars)
#   # NOTE: exp(xb), because xb is log size, but model fit using UNTRANSFORMED cover
#   n_recr  <- (exp(xb)/100) * exp( pars$fecu_b0 )
#   return( n_recr )
# }

# Size distribution of recruits
# This is required in case the majority/all of the size distribution is concentrated between
# the midpoint of a bin, and its upper/lower bounds. In that case, recruits go to 0!
recruits<-function(y,pars,h){
  
  # calculate proportion of seedlings within each bin
  below_upper <- pnorm(y+(h*0.5), pars$recr_sz, pars$recr_sd) 
  below_lower <- pnorm(y-(h*0.5), pars$recr_sz, pars$recr_sd)
  below_upper - below_lower
  
}

# recruitment function: production of y-sized recruits from x-sized mothers
fxy<-function(x,y,pars,h){
  
  xb      <- x_range(x, pars)
  # 1. exp(xb), because xb is log size, but model fit using UNTRANSFORMED cover. 
  # 2. 100 because "per 100cm2"
  n_recr  <- (exp(xb)/100) * exp( pars$fecu_b0 )
  # distribution of recruitment size
  recr_y  <- recruits(y,pars,h)
  # recr_y  <- dnorm(y, mean = pars$recr_sz, sd = pars$recr_sd)
  # recr_y  <- dnorm(y,mean = pars$recr_sz, sd = pars$recr_sd)/(1-pnorm(-1.61,pars$recr_sz,pars$recr_sd))
  # fac2    <- ((y-pars$recr_sz)^2)/(2*pars$recr_sd^2)
  # fac1    <- sqrt(2*pi)*pars$recr_sd #using kidsize mean & SD from above
  #f       <- n_recr*exp(-fac2)/fac1;
  # recruits and size distribution 
  # return(recr_y*n_recr);
  f       <- n_recr*recr_y
  return(f)
  
}



# IPM kernel/matrix ------------------------------------------------------------
kernel <- function(pars){
  
  n   <- pars$mat_siz
  L   <- pars$L 
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins 
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # Fertility matrix
  Fmat        <- matrix(0, n, n)
  
  # Populate fecundity matrix
  # Fmat[1,1:n] <- fx(y,pars)
  Fmat[]      <- t( outer(y,y,fxy,pars,h) )
  
  # Growth/survival transition matrix
  Tmat        <- matrix(0,n,n)
  
  # Growth/survival transitions among cts sizes
  Tmat[]      <- t( outer(y,y,pxy,pars) )*h
  
  # recruits  
  # Tmat[1:n,1] <- Tmat[1:n,1] + recruits(y,pars)*h   
  
  # Full Kernel is simply a summation of fertility and transition matrix
  k_yx        <- Fmat + Tmat
  
  return(list(k_yx    = k_yx,
              Fmat    = Fmat,
              Tmat    = Tmat,
              meshpts = y))
  
}

# check lambda of mean matrix
lam_mean <- Re(eigen(kernel(pars_mean)$k_yx)$value[1])


# year- and group-specific matrices -----------------------------------------------------------


# change group name for recruitment models

if( do_spp == 'ARTR'){
# start from recruitment data
recD        <- read.csv( paste(cur_dir,"/speciesData/",doSpp,"/recArea.csv",sep="") ) 
gr_ordrd    <- recD$Group %>% as.factor %>% unique %>% sort %>% as.character
new_groups  <- paste0('group_',gr_ordrd)

# ugly, but it works
recr_p      <- recr_p %>% 
            mutate( coefficient = replace(coefficient, coefficient == 'intcpt.gr1', new_groups[1]),
                    coefficient = replace(coefficient, coefficient == 'intcpt.gr2', new_groups[2]),
                    coefficient = replace(coefficient, coefficient == 'intcpt.gr3', new_groups[3]),
                    coefficient = replace(coefficient, coefficient == 'intcpt.gr4', new_groups[4]),
                    coefficient = replace(coefficient, coefficient == 'intcpt.gr5', new_groups[5]),
                    coefficient = replace(coefficient, coefficient == 'intcpt.gr6', new_groups[6]) )
}

# extract years fitted in models 
extract_yr <- function(param_df){
  subset(param_df, grepl("year_",coefficient) ) %>%
  mutate( coefficient = gsub("year_","",coefficient) ) %>%
  .$coefficient
}

# extract group names fitted for survival and growth model 
extract_gr <- function(param_df){
  subset(param_df, grepl("group_|intcpt.gr",coefficient) ) %>%
  mutate( coefficient = gsub("group_|intcpt.gr","",coefficient) ) %>%
  .$coefficient
}

# create grid of yr and groups
expand_yr_gr <- function(coef_df){
  expand.grid( yr = extract_yr(coef_df), 
               gr = extract_gr(coef_df), 
               stringsAsFactors = F ) %>% 
    arrange(yr,gr)
}

# list of year by group combinations
yr_by_gr_l    <- list( expand_yr_gr(grow_p), 
                       expand_yr_gr(surv_p), 
                       expand_yr_gr(recr_p) )

# test: coefficients associated with 3 vital rates all have same year/group combinations 
expect_true( sapply(yr_by_gr_l, # nested s apply
                    function(x) sapply(yr_by_gr_l, 
                                       function(y) all.equal(x,y))) %>% all )
 

# year- and group-specific parameters
yr_par <- function(pars_mean, yr_gr){
  
  # find year and group
  yr_grab           <- yr_gr$yr
  gr_grab           <- yr_gr$gr
  
  # extract year- and group- specific intercepts
  surv_yr_b0        <- subset(surv_p, coefficient == paste0("year_",   yr_grab) )$value
  surv_gr_b0        <- subset(surv_p, coefficient == paste0("group_",  gr_grab) )$value
  surv_yr_b1        <- subset(surv_p, coefficient == paste0("logarea_",yr_grab) )$value
  
  grow_yr_b0        <- subset(grow_p, coefficient == paste0("year_",   yr_grab) )$value
  grow_gr_b0        <- subset(grow_p, coefficient == paste0("group_",  gr_grab) )$value
  grow_yr_b1        <- subset(grow_p, coefficient == paste0("logarea.t0_",yr_grab) )$value
  
  recr_yr_b0        <- subset(recr_p, coefficient == paste0("year_",   yr_grab) )$value
  recr_gr_b0        <- subset(recr_p, coefficient == paste0("group_",  gr_grab) )$value
    
  # assign year- and group- specific parameters
  pars_mean$surv_b0 <- extr_value(surv_p, "b0") + surv_yr_b0 + surv_gr_b0
  pars_mean$surv_b1 <- extr_value(surv_p, "logarea") + surv_yr_b1
  pars_mean$grow_b0 <- extr_value(grow_p, "b0") + grow_yr_b0 + grow_gr_b0
  pars_mean$grow_b1 <- extr_value(grow_p, "logarea.t0") + grow_yr_b1
  pars_mean$fecu_b0 <- extr_value(recr_p, "b0") + recr_yr_b0 + recr_gr_b0
  pars_out          <- pars_mean
  
  return(pars_out)
  
}

# store matrices on your computer
store_mats <- function(ii, pars_mean, do_spp, yr_gr_combin){
  
  # identify year and group 
  yr_gr   <- yr_gr_combin[ii,]
  
  # get year- and group-specific parameters, and produce matrices
  pars_in <- yr_par(pars_mean, yr_gr)
  mat_A   <- kernel(pars_in)$k_yx
  mat_T   <- kernel(pars_in)$Tmat
  mat_F   <- kernel(pars_in)$Fmat
  
  # past year and group
  yr_out  <- paste(yr_gr, collapse='_')
  
  # create species folder
  write.csv(mat_A, paste0(do_spp,"/matA_",yr_out,".csv"), row.names=F)
  write.csv(mat_T, paste0(do_spp,"/matT_",yr_out,".csv"), row.names=F)
  write.csv(mat_F, paste0(do_spp,"/matF_",yr_out,".csv"), row.names=F)
  
}

# all year by group combinations
yr_gr_combin <- yr_by_gr_l[[1]] %>% left_join( treat )

# store matrices!
# lapply(1:nrow(yr_gr_combin), store_mats, pars, do_spp, yr_gr_combin)


# tests -----------------------------------------------------

# calculate year-specific lambdas on the go 
det_lam <- function(ii, pars_mean, yr_gr_combin){
  yr_gr     <- yr_gr_combin[ii,]
  pars_feed <- yr_par(pars_mean, yr_gr)
  Re(eigen(kernel(pars_feed)$k_yx)$value[1])
}
lam_vec <- sapply(1:nrow(yr_gr_combin), det_lam, pars_mean, yr_gr_combin)

# # calculate stored lambdas
# file_l <- Filter(function(x) grepl("matA_",x), list.files(do_spp) )
# mat_l  <- lapply(file_l, function(x) read.csv(paste0(do_spp,"/",x)) )
# lam_st <- sapply(mat_l, function(x) Re(eigen(x)$value[1]) )
# 
# # test
# expect_true( all(round(lam_vec,8) == round(lam_st,8) ) )

# # calculate MEAN MATRIX and associated lambda
# mat_arr   <- array(unlist(mat_l), dim = c(10, 10, length(mat_l)) )
# mean_mat  <- apply(mat_arr, c(1,2), mean)
# 
# Re(eigen(mean_mat)$value[1])

# metadata names of stacked matrices
stck_nams <- mat_stck %>% names %>% .[1:32]

# store matrices on your computer
compadre_form <- function(ii, pars_mean, do_spp, yr_gr_combin){
  
  # identify year and group 
  yr_gr   <- yr_gr_combin[ii,]
  tr      <- subset(treat, gr == yr_gr$gr)$MatrixTreatment
  
  # metadata of stacked matrices 
  meta_df <- data.frame ( EnteredBy	  = 'Aldo Compagnoni',
                          EnteredDate = '8.3.2018',
                          Source      = 'Other source',
                          SpeciesAuthor = subset(spp_df, spp_code == do_spp)$SpeciesAuthor,
                          StudiedSex  = "A",
                          
                          MatrixComposite = "Individual",
                          MatrixTreatment = tr,
                          MatrixCaptivity = "W",
                          #MatrixStartYear = 1900 + as.numeric(yr_gr$yr) - 1,
                          MatrixStartSeason = 2, 
                          MatrixStartMonth = 6,
                          #MatrixEndYear = 1900 + as.numeric(yr_gr$yr),
                          MatrixEndSeason = 2,
                          MatrixEndMonth = 6,
                          MatrixPopulation = yr_gr$gr,
                          
                          LatDeg = 44,
                          LatMin = 1,
                          LatSec = 59.9988,
                          LatNS = NA,
                          LonDeg = 112,
                          LonMin = 1,
                          LonSec = 0.0012,
                          LonWE  = NA,
                          Altitude = NA, 
                          Country = "USA", 
                          Continent = "N America",
                          
                          MatrixSplit = "Divided",
                          MatrixFec = "Yes",
                          Observation = NA )

  # introduce NAs
  meta_df[2:10,] <- NA
  
  # get year- and group-specific parameters, and produce matrices
  pars_in <- yr_par(pars_mean, yr_gr)
  
  # matrices
  form_m  <- function(mat) mat %>% as.data.frame %>% tibble::add_column( A11 = NA )
  mat_A   <- kernel(pars_in)$k_yx %>% form_m
  mat_T   <- kernel(pars_in)$Tmat %>% form_m
  mat_F   <- kernel(pars_in)$Fmat %>% form_m
  
  # lower/upper mesh points
  mesh_p  <- kernel(pars_in)$meshpts
  h2      <- (mesh_p[2] - mesh_p[1])/2
  mesh_low<- mesh_p - h2
  mesh_up <- mesh_p + h2
  
  # matrix data frames
  mats_df <- Reduce(function(...) bind_cols(...), list(mat_A,mat_T,mat_F) ) %>% 
                setNames( paste0('A',1:33) ) %>% 
                mutate( MatrixClassNumber = c(1:10),
                        MatrixClassOrganized = 'active',
                        MatrixStartYear   = paste0('19', yr_gr$yr),
                        MatrixEndYear     = paste0('19', as.numeric(yr_gr$yr) + 1),
                        MatrixClassAuthor = paste0('log(area): ',
                                                   round(mesh_low,4), ' - ', round(mesh_up,4) )
                        )
  
  # introduce missing values (if needed)
  out     <- meta_df %>% 
                bind_cols( mats_df ) %>% 
                select( c(stck_nams, paste0('A',1:33)) )
  
  return(out)
  
}

# store matrices!
compadre_l  <- lapply(1:nrow(yr_gr_combin), compadre_form, pars_mean, do_spp, yr_gr_combin)
compadre_df <- do.call(rbind, compadre_l)

setdiff(meta_col, names(compadre_df) )

# final matrix stacked file
write.csv(compadre_df, paste0(do_spp,'/',do_spp,'_COMPADRE.csv'), row.names=F)


# species descriptors data frame
spp_d_df <- data.frame( Student = 'Aldo Compagnoni',
                        Notes   = 'Density independent matrices created from published data',         
                        VersionRelease = 'X',     
                        Transferred    = 'Yes',
                        MissingData    = 'S',
                        DateDigitization = '11.21.2018', 
                        Embargo          = NA,
                        StudyChecked     = 0,
                        StudyDifficulty  = 2, 
                        OrganismType     = 'Herbaceous perennial',
                        SpeciesAuthor    = subset(spp_df, spp_code == do_spp)$SpeciesAuthor,
                        CommonName       = NA,
                        Kingdom          = NA,
                        Phylum           = NA,
                        Class            = NA, 
                        Order            = NA,
                        Family           = NA,
                        Genus            = NA,
                        SpeciesAccepted  = NA,
                        DateOfCheck      = NA,
                        AngioGymno       = 'Angiosperm',
                        DicotMonoc       = 'Monocot',
                        Authors          = 'Chu; Kleinhesselink; Havstad; McClaran; Peters; Vermeire; Wei; Adler',
                        Journal          = 'Nat Commun',
                        YearPublication  = 2016, 
                        CorrespondingAuthor = 'Chengjin Chu',
                        Email            = 'chuchjin@mail.sysu.edu.cn',
                        Contacted        = NA,
                        Contacted_again  = NA,
                        ContentEmail     = NA,
                        Reply            = NA,
                        DOI.ISBN         = '10.1038/ncomms11766',     
                        AdditionalSource = NA,
                        stringsAsFactors = F)
  
# final matrix stacked file
write.csv(spp_d_df, paste0(do_spp,'/',do_spp,'_spp_descriptors.csv'), row.names=F)
