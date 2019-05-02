# 1. Reproduce Dangremond's MPM ABBOTTS parameters
# 2. Reproduce Dangremond's MPM NORTH TOWER parameters
# 3. Can I recreate average number of racemes using "lupine_all.csv", 
# AND Pop1_2005-2008.csv?
# 4. Reproduce MPMs: ABBOTTS and RADIO TOWER from 2005-2008.
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

# 1. Reproduce Dangremond's MPM parameters ------------------------------

# ABBOTTS (AL (1)) 2005-2008 data
pop1    <- read.csv("data/pre_2008/Pop1_2005-2008.csv", stringsAsFactors = F) %>% 
              # fix clear mistake
              mutate( width = replace(width, width == 660, 66) )

# length with brances 
a_vs_b  <- read.csv('data/lupine_area_vs_branches.csv') %>% 
              mutate( area = (LENGTH/2) * (WIDTH/2) * pi ) %>% 
              mutate( NUMBRANCHES2 = NUMBRANCHES^2 )

# models: area versus branches
m_lin  <- lm(area ~ NUMBRANCHES, data = a_vs_b)
m_quad <- lm(area ~ NUMBRANCHES + NUMBRANCHES2, data = a_vs_b)

# quadratic is better
AIC(m_lin, m_quad)

# model on size
mod    <- m_quad

# size (branches)
size  <- pop1 %>%
            # format width/length
            mutate( length = replace(length, length == "",  NA) ) %>% 
            mutate( length = replace(length, length == ".", NA) ) %>%
            mutate( length = replace(length, length == "barely alive", NA) ) %>%
            mutate( width  = replace(width,  width  == "",  NA) ) %>%
            mutate( width  = replace(width,  width  == ".", NA) ) %>%
            mutate( length = as.numeric(length),
                    width  = as.numeric(width) ) %>%
            # grab size-relevant vars
            dplyr::select(UID, population, location, GPS,
                          X..branches.2005,branches.2006,branches,length,width) %>%
                      rename( brances.2007 = branches ) %>%
                      mutate( area.2008    = (length/2) * (width/2) * pi )

# stage
stage <- pop1 %>% dplyr::select(UID, population, location, GPS,
                                grep("stage", names(.)) )

# flowers
flow  <- pop1 %>% dplyr::select(UID, population, location, GPS,
                                X.racemes.2005, X.racemes.2006, racemes, racemes.1) %>%
            rename( notab.2005 = X.racemes.2005, 
                    notab.2006 = X.racemes.2006, 
                    notab.2007 = racemes, 
                    notab.2008 = racemes.1 )

# DO NUMBER CLIPPED (numcl)
clip  <- pop1 %>% dplyr::select(UID, population, location, GPS,
                                X.rac_eaten, 
                                number.racemes.eaten.by.mice, 
                                number.racemes.eaten.by.mice.1, 
                                number.racemes.eaten.by.mice.2) %>% 
            rename( numcl.2005 = X.rac_eaten, 
                    numcl.2006 = number.racemes.eaten.by.mice, 
                    numcl.2007 = number.racemes.eaten.by.mice.1, 
                    numcl.2008 = number.racemes.eaten.by.mice.2 )               
     

# get data in long form
long_data <- function(x, value, vars){
  x %>% gather_( key = "year",
                 value = value,
                 vars )
}

# remove "." and "" observations
dot_blank_na <- function(x){
  replace(x, x == "." | x == "", NA) 
}

# format year: make it numeric by removing punctuation and letters
year_form <- function(x){ 
  gsub("[A-z]|[[:punct:]]","",x) %>% as.numeric
}


# stage: long format
stage_long  <- long_data(stage, "stage", paste0('stage.',2005:2008) ) %>%
                  mutate( stage = dot_blank_na(stage) ) %>%
                  mutate( year  = year_form(year) )

# size long format; removed "regular; large; XL" 
size_long   <- size %>% 
                  dplyr::select(-width, -length) %>%
                  long_data("size",    c("X..branches.2005", "branches.2006", 
                                         "brances.2007",     "area.2008") ) %>%
                  mutate( size = dot_blank_na(size) ) %>%
                  mutate( size = replace(size, size == "regular", NA) ) %>%
                  mutate( size = replace(size, size == "large",   NA) ) %>%
                  mutate( size = replace(size, size == "XL",      NA) ) %>%
                  mutate( size = as.numeric(size) ) %>%
                  # transform n. of branches in area
                  mutate( area = coef(mod)[1] + 
                                 coef(mod)[2] * size + 
                                 coef(mod)[3] * size^2 ) %>% 
                  mutate( year = year_form(year) ) %>% 
                  # size == area in 2008
                  mutate( area = replace(area, 
                                         year == 2008, 
                                         size[year == 2008]) )

# fecundity data in long form  
fecu_long   <- long_data( flow,  "notab", paste0("notab.", 2005:2008) ) %>%
                  mutate( notab = dot_blank_na(notab) ) %>%
                  mutate( year  = year_form(year),
                          notab = as.numeric(notab) )

# clipped data in long form  
clip_long   <- long_data( clip,  "numcl", paste0("numcl.", 2005:2008) ) %>%
                  mutate( numcl = dot_blank_na(numcl) ) %>%
                  mutate( year  = year_form(year),
                          numcl = as.numeric(numcl) )

# combine all demographic data in one file
demo_pop1   <- Reduce(function(...) full_join(...), list( stage_long, 
                                                          size_long,
                                                          fecu_long,
                                                          clip_long) )
# demography at time t0
demo_t0    <- demo_pop1 %>%
                setNames( c(names(.)[1:5],
                            paste0(names(.)[-c(1:5)],"_t0")) )
# demography at time t1
demo_t1    <- demo_pop1 %>%
                setNames( c(names(.)[1:5],
                            paste0(names(.)[-c(1:5)],"_t1")) ) %>%
                mutate( year = year - 1)

# transition data frame 
trans_pop1 <- full_join( demo_t0, demo_t1 ) %>%
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "non", "Non") ) %>%  
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "dead", "Dead") ) %>%  
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "rep", "Rep") ) %>%
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "non", "Non") ) %>%  
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "dead", "Dead") ) %>%  
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "rep", "Rep") ) %>%
                # introduce survival at time t1
                mutate( surv_t1 = 1 ) %>%
                mutate( surv_t1 = replace(surv_t1, stage_t1 == "Dead", 0)  ) %>%
                mutate( surv_t1 = replace(surv_t1, is.na(stage_t0),    NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, stage_t1 == "SL",   NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, is.na(stage_t1),    NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, stage_t1 == "not found", NA) ) %>%
                # introduce flowering at time t0 and t1
                mutate( flow_t0 = 0,
                        flow_t1 = 0) %>% 
                mutate( flow_t0 = replace(flow_t0, stage_t0 == "Rep", 1),
                        flow_t1 = replace(flow_t1, stage_t1 == "Rep", 1) ) %>%
                # introduce NAs: everything that is NOT "Rep" or "Non",
                mutate( flow_t0 = replace(flow_t0, 
                                          !(stage_t0 == "Rep" | stage_t0 == "Non"), 
                                          NA),
                        flow_t1 = replace(flow_t1, 
                                          !(stage_t1 == "Rep" | stage_t1 == "Non"), 
                                          NA) ) %>%
                # NAs whenever is.na(stage_t) == TRUE
                mutate( flow_t0 = replace(flow_t0, is.na(stage_t0), NA) ) %>%
                mutate( flow_t1 = replace(flow_t1, is.na(stage_t1), NA) )
    

# transition rates in 2005-2006 
transition <- function(x, start, finish, tr_year){
  
  # starting abundances
  init_n <- x %>% 
              subset( year == tr_year ) %>% 
              subset( stage_t0 == start) %>% 
              nrow
  
  # ending abundances
  fin_n  <- x %>% 
              subset( year == tr_year  ) %>% 
              subset( stage_t0 == start ) %>% 
              subset( stage_t1 == finish ) %>% 
              nrow
  
  fin_n / init_n
  
}

# 2006 transitions are dead on
transition(trans_pop1, "SL","Non", 2006)
transition(trans_pop1, "SL","Rep", 2006)
transition(trans_pop1, "Non","Non", 2006)
transition(trans_pop1, "Non","Rep", 2006)
transition(trans_pop1, "Rep","Non", 2006)
transition(trans_pop1, "Rep","Rep", 2006)

trans_pop1$stage_t0 %>% unique %>% sort
trans_pop1$stage_t1 %>% unique %>% sort

rr <- trans_df %>% 
        subset( year == year_in ) %>% 
        subset( stage_t0 == 'Rep') %>% 
        .$notab_t0 %>% 
        mean(na.rm=T)
    

g1  <- 0.011
g2  <- 0.034
g3  <- 0.015
vv  <- 3.563
ff  <- 4.599
cm1 <- 0.257
            

rr * ff * vv * cm1 * g1
rr * ff * vv * cm1 * g2
rr * ff * vv * cm1 * g3
  
c(0.3708 / (7.161*4.599*3.563*0.257),
  1.1250 / (7.161*4.599*3.563*0.257),
  0.5057 / (7.161*4.599*3.563*0.257) ) %>% round(4)



# 2. Reproduce Dangremond's MPM NORTH TOWER parameters -----------------------------
pop8    <- read.csv("data/pre_2008/Pop8_2005-2008.csv",  
                    stringsAsFactors = F)[1:168,]

# stage
stage8_long <- pop8 %>% 
                dplyr::select(UID, new.tag.number, population, plant.., GPS, 
                              c("stage.2005", "stage.2006", 
                                "X2007", "stage.2008") ) %>%
  
                rename( stage.2007 = X2007 ) %>%
                long_data( "stage", paste0("stage.", 2005:2008) ) %>%
                mutate( stage = dot_blank_na(stage) ) %>%
                mutate( year = year_form(year) )

# size (branches)
size8_long  <- pop8 %>% # grab 
                  dplyr::select(UID, new.tag.number, population, plant.., GPS, 
                                X..racemes.2005, racemes.2006, racemes, length,width) %>%
                  rename( racemes.2005 = X..racemes.2005,
                          racemes.2007 = racemes ) %>%
                  mutate( area.2008    = (length/2) * (width/2) * pi ) %>%
                  long_data( "size", c( paste0("racemes.",2005:2007), "area.2008") ) %>%
                  mutate( year = year_form(year) ) %>% 
                  # transform n. of branches in area
                  mutate( area = coef(mod)[1] + 
                                 coef(mod)[2] * size + 
                                 coef(mod)[3] * size^2 ) %>% 
                  # size == area in 2008
                  mutate( area = replace(area, 
                                         year == 2008, 
                                         size[year == 2008]) )
            
# flowers
fecu8_long  <- pop8 %>% 
                  dplyr::select(UID, new.tag.number, population, plant.., GPS, 
                                X.fl.stalks.2005,  X.fl.stalks.2006, flower.stalks, racemes.1) %>%
                  rename( notab.2005 = X.fl.stalks.2005,  
                          notab.2006 = X.fl.stalks.2006, 
                          notab.2007 = flower.stalks, 
                          notab.2008 = racemes.1 ) %>% 
                  long_data("notab", paste0("notab.", 2005:2008) ) %>%
                  mutate( notab   = dot_blank_na(notab) ) %>%
                  mutate( year    = year_form(year) )


# clipped data in long form  
clip_long8  <- pop8 %>% 
                  dplyr::select(UID, new.tag.number, population, plant.., GPS, 
                                number.stalks.eaten.by.mice,
                                number.stalks.eaten.by.mice.1,
                                number.stalks.eaten.by.mice.2) %>% 
                  rename( numcl.2005 = number.stalks.eaten.by.mice,
                          numcl.2006 = number.stalks.eaten.by.mice.1,
                          numcl.2007 = number.stalks.eaten.by.mice.2 ) %>%  
                  long_data( "numcl", paste0("numcl.", 2005:2007) ) %>%
                  mutate( numcl = replace(numcl, numcl == '?', NA) ) %>% 
                  mutate( numcl = dot_blank_na(numcl) ) %>%
                  mutate( year  = year_form(year),
                          numcl = as.numeric(numcl) )

# combine all demographic data in one file
demo_pop8   <- Reduce( function(...) full_join(...), list( stage8_long, 
                                                           size8_long,
                                                           fecu8_long,
                                                           clip_long8) ) 
# demography at time t0
demo8_t0    <- demo_pop8 %>%
                setNames( c(names(.)[1:6],
                            paste0(names(.)[-c(1:6)],"_t0")) ) %>%
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "non", "Non") ) %>%  
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "rep", "Rep") )

# demography at time t1
demo8_t1    <- demo_pop8 %>%
                setNames( c(names(.)[1:6],
                            paste0(names(.)[-c(1:6)],"_t1")) ) %>%
                mutate( year = year - 1 )
                
# transition data frame 
trans_pop8 <- full_join(dplyr::select(demo8_t0, -width_t0, -length_t0),
                        dplyr::select(demo8_t1, -width_t1, -length_t1) ) %>% 
                rename( plant_n = plant..) %>%
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "seedling", "SL") ) %>%
                mutate( stage_t0 = replace(stage_t0, stage_t0 == "dead", "Dead") ) %>%
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "seedling", "SL") ) %>%
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "non", "Non") ) %>%  
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "dead", "Dead") ) %>%  
                mutate( stage_t1 = replace(stage_t1, stage_t1 == "rep", "Rep") ) %>%
                # introduce survival at time t1
                mutate( surv_t1 = 1 ) %>%
                mutate( surv_t1 = replace(surv_t1, stage_t1 == "Dead", 0) ) %>%
                mutate( surv_t1 = replace(surv_t1, is.na(stage_t1), NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, stage_t1 == "not found", NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, stage_t1 == "SL", NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, 
                                          stage_t0 == "Dead" & stage_t1 == "Dead", 
                                          NA) ) %>%
                mutate( surv_t1 = replace(surv_t1, 
                                          is.na(stage_t0) & (stage_t1 == "Non" | stage_t1 == "Rep"), 
                                          NA) ) %>%
                # introduce flowering at time t0 and t1
                mutate( flow_t0 = 0,
                        flow_t1 = 0) %>% 
                mutate( flow_t0 = replace(flow_t0, stage_t0 == "Rep", 1),
                        flow_t1 = replace(flow_t1, stage_t1 == "Rep", 1) ) %>%
                # introduce NAs: everything that is NOT "Rep" or "Non",
                mutate( flow_t0 = replace(flow_t0, 
                                          !(stage_t0 == "Rep" | stage_t0 == "Non"), 
                                          NA),
                        flow_t1 = replace(flow_t1, 
                                          !(stage_t1 == "Rep" | stage_t1 == "Non"), 
                                          NA) ) %>%
                # NAs whenever is.na(stage_t) == TRUE
                mutate( flow_t0 = replace(flow_t0, is.na(stage_t0), NA) ) %>%
                mutate( flow_t1 = replace(flow_t1, is.na(stage_t1), NA) )


trans_pop8$stage_t0 %>% unique
trans_pop8$stage_t1 %>% unique

# 2005 transitions are dead on
# 2007 transitions are very off.
transition(trans_pop8, "SL","Non", 2007)
transition(trans_pop8, "SL","Rep", 2007)
transition(trans_pop8, "Non","Non", 2007)
transition(trans_pop8, "Non","Rep", 2007)
transition(trans_pop8, "Rep","Non", 2007)
transition(trans_pop8, "Rep","Rep", 2007)


# check whether there is something wrong with original "stages" data
debug_8 <- stage8_long %>% 
              mutate( stage = replace(stage,
                                      stage == 'rep',
                                      'Rep') ) %>% 
              mutate( stage = replace(stage,
                                      stage == 'non',
                                      'Non') ) %>%
              mutate( stage = replace(stage,
                                      stage == 'seedling',
                                      'SL') ) %>% 
              mutate( stage = replace(stage,
                                      stage == 'dead',
                                      'Dead') ) %>% 
              mutate( year = paste0('yr_',year) ) %>% 
              spread(year,stage)


# SL to Non
init_n <- debug_8 %>% 
            subset( yr_2007 == 'SL' ) %>% 
            nrow
# ending abundances
fin_n  <- debug_8 %>% 
            subset( yr_2007 == 'SL' & yr_2008 == 'Non' ) %>% 
            nrow
fin_n / init_n
transition(trans_pop8, "SL","Non", 2007)

# SL to Rep
fin_n  <- debug_8 %>% 
            subset( yr_2007 == 'SL' & yr_2008 == 'Rep' ) %>% 
            nrow
fin_n / init_n
transition(trans_pop8, "SL","Rep", 2007)

# Non to Non
init_n <- debug_8 %>% 
            subset( yr_2007 == 'Non' ) %>% 
            nrow
# ending abundances
fin_n  <- debug_8 %>% 
            subset( yr_2007 == 'Non' & yr_2008 == 'Non' ) %>% 
            nrow
fin_n / init_n
transition(trans_pop8, "Non","Non", 2007)

# Non to Rep
fin_n  <- debug_8 %>% 
            subset( yr_2007 == 'Non' & yr_2008 == 'Rep' ) %>% 
            nrow
fin_n / init_n
transition(trans_pop8, "Non","Rep", 2007)


# 3. recreate average n. of racemes -----------------------------
year_in <- 2007

# lupine_all.csv
read.csv( "data/lupine_all.csv") %>% 
  subset( year == year_in ) %>% 
  subset( location == 'AL (1)' ) %>% 
  subset( !(numrac_t0 %in% 0) ) %>% 
  # subset( area_t0 != 1 ) %>% 
  .$numrac_t0 %>% 
  mean(na.rm=T)

replace(pop1$X.racemes.2005, 
        pop1$X.racemes.2005 == '.' | 
        pop1$X.racemes.2005 == '', 
        NA) %>% 
  as.numeric %>% 
  mean(na.rm=T)

replace(pop1$X.racemes.2006, 
        pop1$X.racemes.2006 == '.' | 
        pop1$X.racemes.2006 == '', 
        NA) %>% 
  as.numeric %>% 
  mean(na.rm=T)

replace(pop1$racemes, 
        pop1$racemes %in% c('.',''), 
        NA) %>% 
  as.numeric %>% 
  mean(na.rm=T)

# 4. ------------------

# combine two data sets
al_att_df <- mutate(trans_pop1, 
                    population = as.character(population) ) %>% 
              bind_rows( mutate(trans_pop8, 
                                population = 'ATT (8)') ) %>% 
              select( -location ) %>% 
              rename( location = population ) %>% 
              # pop 1 and 8 are "AL (1)" and "ATT (8)"
              mutate( location = replace(location, 
                                         location == '1',
                                         'AL (1)') )

  
# site-year specific fertility parameters
site_fec_df <- data.frame( site = c( rep('AL (1)', 3),
                                     rep('ATT (8)',3) ),
                           year = c(2005, 2006, 2007,
                                    2005, 2006, 2007 ),
                           ff   = c(3.025, 4.599, 1.286, 
                                    2.837, 5.198, 5.872),
                           cm1  = c(0.300, 0.257, 0.180,
                                    0.495, 0.667, 0.586) ) 


# produce matrix
mpm <- function(site_in, year_in){
  
  # select site and year info
  trans_df <- subset(al_att_df, location == site_in )
  fert_df  <- subset(site_fec_df, site == site_in & year == year_in )
    
  # common parameters
  g1  <- 0.011
  g2  <- 0.034
  g3  <- 0.015
  vv  <- 3.563
    
  # site- and year-specific fertility parameters
  rr  <- trans_df %>% 
           subset( year == year_in ) %>% 
           subset( stage_t0 == 'Rep') %>% 
           .$notab_t0 %>% 
           mean(na.rm=T)
  ff  <- fert_df$ff
  cm1 <- fert_df$cm1
  
  # set up the matrix
  mat       <- matrix(0, 5, 5)
  mat[2,1]  <- 1
  mat[3,2]  <- 1
  
  # transition
  mat[4,3]  <- transition(trans_df, "SL","Non",  year_in)
  mat[5,3]  <- transition(trans_df, "SL","Rep",  year_in)
  mat[4,4]  <- transition(trans_df, "Non","Non", year_in)
  mat[5,4]  <- transition(trans_df, "Non","Rep", year_in)
  mat[4,5]  <- transition(trans_df, "Rep","Non", year_in)
  mat[5,5]  <- transition(trans_df, "Rep","Rep", year_in)

  # fertility parameters
  mat[1,5]  <- rr * ff * vv * cm1 * g3
  mat[2,5]  <- rr * ff * vv * cm1 * g2
  mat[3,5]  <- rr * ff * vv * cm1 * g1
  
  mat
}


# compare with Dangermond 2010

# download published Dangremond et al. 2010 MPMs
al_mpm_l  <- lapply(1:3, function(ii){ read_xlsx('data/dangremond/abbotts_ambient_mpm.xlsx',
                                                 col_names = F,
                                                 sheet=ii) %>% 
                                                  as.matrix })
att_mpm_l <- lapply(1:3, function(ii){ read_xlsx('data/dangremond/radiotower_ambient_mpm.xlsx',
                                                 col_names = F,
                                                 sheet=ii) %>% 
                                                  as.matrix })

# calculate lambda from list of matrices
calc_lam <- function(x){
  eig_i <- which( Re(eigen(x)$value) == max(Re(eigen(x)$value)) )
  eigen(x)$value[eig_i] %>% Re
}

# Our test matrices  
al_test_l  <- list( mpm('AL (1)', 2005), 
                    mpm('AL (1)', 2006), 
                    mpm('AL (1)', 2007) )

att_test_l <- list( mpm('ATT (8)', 2005),
                    mpm('ATT (8)', 2006),
                    mpm('ATT (8)', 2007) )


# make data frame to "bind" for plotting   
make_df <- function(site_in, source_in, lambda_v){
  
  data.frame( site   = site_in,
              source = source_in,
              year   = c(2005:2007),
              lambda = lambda_v )
  
}
  

lam_df <- list( make_df('AL',  'published',  sapply(al_mpm_l,   calc_lam) ),
                make_df('ATT', 'published',  sapply(att_mpm_l,  calc_lam) ),
                make_df('AL',  'reproduced', sapply(al_test_l,  calc_lam) ),
                make_df('ATT', 'reproduced', sapply(att_test_l, calc_lam) ) ) %>% 
            bind_rows 

# store results   
ggplot(lam_df) +
  geom_line( aes(x        = year,
                 y        = lambda,
                 color    = source,
                 linetype = site ),
             size    = 1.5,
             lineend = 'round' ) +
  scale_color_viridis_d() + 
  ylab( expression(lambda) ) + 
  theme( axis.text.x = element_text( angle = 70) ) +
  ggsave('results/ipm/validation_mpm/mpm_published_reconstructed.tiff',
         width = 6.3, height = 6.3, compression = 'lzw')

# check what's wrong with 2005 in Abbots
pop1 %>% 
  select(stage.2005, stage.2006)

al_mpm_l[[1]]
mpm('AL (1)', 2005)

sum(pop1$stage.2005  == 'Rep')
sum(pop1$stage.2005  == 'Rep' & 
    (pop1$stage.2006 == 'Rep' | pop1$stage.2006 == '') )

