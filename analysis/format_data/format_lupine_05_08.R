library(dplyr)
library(tidyr)

# 2005-2008 data
pop1    <- read.csv("data/pre_2008/Pop1_2005-2008.csv", stringsAsFactors = F) %>% 
              # fix clear mistake
              mutate( width = replace(width, width == 660, 66) )
pop8    <- read.csv("data/pre_2008/Pop8_2005-2008.csv",  stringsAsFactors = F)[1:168,]

# post 2008 data, for comparison
post_d  <- read.csv("data/from_2008_to_2017/LUTI_0817_AllPops_July212017.csv")

# length with brances 
a_vs_b  <- read.csv('data/lupine_area_vs_branches.csv') %>% 
              mutate( area = (LENGTH/2) * (WIDTH/2) * pi ) %>% 
              mutate( NUMBRANCHES2 = NUMBRANCHES^2 )

# Mismatches between plant.. and plant...1
sum(pop8$plant.. == pop8$plant...1) - nrow(pop8)
id <- which(pop8$plant.. != pop8$plant...1)
head(pop8[id,])
pop8$new.tag.number %>% unique %>% length


# area versus branches ---------------------------------------------------------

# models
m_lin  <- lm(area ~ NUMBRANCHES, data = a_vs_b)
m_quad <- lm(area ~ NUMBRANCHES + NUMBRANCHES2, data = a_vs_b)

# quadratic is better
AIC(m_lin, m_quad)

# model on size
mod    <- m_quad


# Population 1 AL (1) ----------------------------------------------------------------------

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
    
# Questions:
# Why stage_t0 == NA, stage_t1 == "Non"?
# Why stage_t0 == NA, stage_t1 == "Rep"?
subset(trans_pop1, is.na(stage_t0) & stage_t1 == "Non") 
subset(trans_pop1, is.na(stage_t0) & stage_t1 == "Rep")
subset(trans_pop1, is.na(stage_t0) & stage_t1 == "SL")
subset(trans_pop1, is.na(stage_t0) & is.na(stage_t1) )
subset(trans_pop1, !is.na(notab_t0) ) 


# Chek survival 
trans_pop1 %>% 
  dplyr::select(stage_t0, stage_t1, surv_t1) %>% 
  unique %>% 
  arrange(stage_t0) 

# Chek flow_t1 
trans_pop1 %>% 
  dplyr::select(stage_t0, stage_t1, flow_t1) %>% 
  unique %>% 
  arrange(stage_t1) 

# Chek flow_t0
trans_pop1 %>% 
  dplyr::select(stage_t0, stage_t1, flow_t0) %>% 
  unique %>% 
  arrange(stage_t0) 

# Chek stages
trans_pop1 %>% 
  dplyr::select(stage_t0, stage_t1) %>% 
  unique %>% 
  arrange(stage_t0) 




# plot size 
plot(area_t1 ~ area_t0, data=subset(trans_pop1))
# plot racemes
plot(notab_t1 ~ area_t1, data = trans_pop1)
plot(notab_t0 ~ area_t0, data = trans_pop1)
# plot survival 
plot( jitter(surv_t1) ~ area_t0, data=subset(trans_pop1, stage_t0 != "SL"))
plot( jitter(surv_t1) ~ area_t0, data=subset(trans_pop1, stage_t0 == "SL"))
# plot flowering
plot( jitter(flow_t0) ~ area_t0, data=subset(trans_pop1) )
plot( jitter(flow_t1) ~ area_t0, data=subset(trans_pop1) )
plot( jitter(flow_t1) ~ area_t1, data=subset(trans_pop1) )

par(mfrow=c(2,2))
# plot size 
plot(area_t1 ~ area_t0, data=subset(trans_pop1))
# plot racemes
plot(notab_t0 ~ area_t0, data = subset(trans_pop1) )
# plot survival 
plot( jitter(surv_t1) ~ area_t0, data=subset(trans_pop1,stage_t0 != "SL"))
# plot flowering
plot( jitter(flow_t0) ~ area_t0, data=subset(trans_pop1))



plot(area_t1 ~ log(area_t0), data = trans_pop1)

# size based on counts
# plot(area_t1 ~ log(area_t0), data = subset(trans_pop1) )
# mod_po    <- lm(    area_t1 ~ area_t0, data = subset(trans_pop1), family="poisson" )
# mod_po_l  <- lm(    area_t1 ~ log(area_t0), data = subset(trans_pop1, year != 2007), family="poisson" )
# mod_nb    <- glm.nb( area_t1 ~ area_t0, data = subset(trans_pop1, year != 2007) )
# mod_nb_l  <- glm.nb( area_t1 ~ log(area_t0), data = subset(trans_pop1, year != 2007) )
# AIC(mod_po,mod_po_l,mod_nb,mod_nb_l)
# subset(trans_pop1, year != 2007)$size_t0 %>% hist


# Population 8 -----------------------------------------------------------------------------

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


# SOLVE ISSUES --------------------------------------------------------------


# ISSUES?!?
subset(trans_pop8, is.na(stage_t0) & stage_t1 == "Non")
subset(trans_pop8, is.na(stage_t0) & stage_t1 == "SL" )
subset(trans_pop8, is.na(stage_t0) & is.na(stage_t1)  )


# Chek survival 
trans_pop8 %>% 
  dplyr::select(stage_t0, stage_t1, surv_t1) %>% 
  unique %>% 
  arrange(stage_t0) 

# Chek flow_t1 
trans_pop8 %>% 
  dplyr::select(stage_t0, stage_t1, flow_t1) %>% 
  unique %>% 
  arrange(stage_t1) 

# Chek flow_t0
trans_pop8 %>% 
  dplyr::select(stage_t0, stage_t1, flow_t0) %>% 
  unique %>% 
  arrange(stage_t0) 

# Chek stages
trans_pop8 %>% 
  dplyr::select(stage_t0, stage_t1) %>% 
  unique %>% 
  arrange(stage_t0) 


# plot size 
plot(area_t1 ~ area_t0, data=subset(trans_pop8, year != 2007))
# plot racemes
plot(racemes_t1 ~ area_t1, data = subset(trans_pop8, year != 2007) )
plot(racemes_t0 ~ area_t0, data = trans_pop8)
# plot survival 
plot( jitter(surv_t1) ~ area_t0, data=subset(trans_pop8, year != 2007))
# plot flowering
plot( jitter(flow_t0) ~ area_t0, data=subset(trans_pop8, year != 2007))
plot( jitter(flow_t1) ~ area_t0, data=subset(trans_pop8, year != 2007))
plot( jitter(flow_t1) ~ area_t1, data=subset(trans_pop8, year != 2007))


# look at possible problematic cases
trans_pop8$stage_t1 %>% unique

subset(trans_pop8, stage_t1 == "Dead") %>% 
  dplyr::select(stage_t0, stage_t1) %>%
  unique
subset(trans_pop8, stage_t1 == "Rep") %>% 
  dplyr::select(stage_t0, stage_t1) %>%
  unique
subset(trans_pop8, stage_t1 == "Non") %>% 
  dplyr::select(stage_t0, stage_t1) %>%
  unique
subset(trans_pop8, stage_t1 == "SL") %>% 
  dplyr::select(stage_t0, stage_t1) %>%
  unique


# population 1 and 8 together -----------------------------------------------  
  
names(trans_pop1)
names(trans_pop8)

# create plant id, convert population to character
trans_pop1 <- trans_pop1 %>%
                mutate( plant_id   = paste( UID, location, sep="_"),
                        population = as.character(population) )
trans_pop8 <- trans_pop8 %>%
                mutate( plant_id = paste( UID, plant_n, sep="_") )

# put it all together and send store
trans_0508 <- bind_rows(trans_pop1, trans_pop8) 

par(mfrow=c(2,2), mar = c(2.5,2.5,.1,0.1))
plot(jitter(surv_t1) ~ area_t0, data = trans_0508)
plot(area_t1    ~ area_t0, data = trans_0508)
plot(jitter(flow_t1) ~ area_t0, data = trans_0508)
plot(notab_t1 ~ log(area_t1), data = trans_0508)
plot(jitter(numcl_t1) ~ jitter(notab_t1), data = trans_0508)
abline(0,1,lwd=2)

trans_0508$notab_t0 %>% hist
trans_0508$notab_t1 %>% hist

# write out -------------------------------------------------

# homogenize 
trans_0508_out <- trans_0508 %>% 
                    # "location" here is a subsample of two populations
                    rename( sub_sample = location ) %>% 
                    # rename location accordingly
                    rename( location = population ) %>% 
                    # pop 1 and 8 are "AL (1)" and "ATT (8)"
                    mutate( location = replace(location, 
                                               location == '1',
                                               'AL (1)') ) %>% 
                    mutate( location = replace(location, 
                                               location == '8a',
                                               'ATT (8)') ) %>% 
                    # is.na(location) is still 'ATT (8)'
                    mutate( location = replace(location, 
                                               is.na(location),
                                               'ATT (8)') ) %>% 
                    # label year as "transition year"
                    mutate( transition = paste(year, year+1, sep = '-') ) %>% 
                    # create NEW ID!
                    rename( newid = new.tag.number ) %>% 
                    mutate( newid = replace(newid, is.na(newid), 
                                            paste0('PRE_08_',1:nrow(subset(., is.na(newid))))) )

# store it
write.csv(trans_0508_out, "data/lupine_05_08.csv", row.names=F)

  


# Older code to identify categories ---------------------------------------------------------
# 
# # pattern is "stage." 
# Filter(function(x) grepl("stage.",x), names(pop1) )
# # X2007, plus pattern "stage."
# Filter(function(x) grepl("X2007",x), names(pop8) )
# Filter(function(x) grepl("stage.",x), names(pop8) )
# 
# # size
# 
# # pattern "branches 
# Filter(function(x) grepl("branches",x), names(pop1) )
# Filter(function(x) grepl("branches",x), names(pop8) )
# Filter(function(x) grepl("racemes.",x), names(pop8) )
