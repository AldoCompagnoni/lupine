setwd("C:/cloud/Dropbox/lupine")
library(dplyr)
library(tidyr)
library(testthat)
library(readxl)

# 2005-2008 data
pop1    <- read.csv("data/pre_2009/Pop1_2005-2008.csv", stringsAsFactors = F) %>% 
              # fix clear mistake
              mutate( width = replace(width, width == 660, 66) )
pop8    <- read.csv("data/pre_2009/Pop8_2005-2008.csv",  stringsAsFactors = F)[1:168,]
pop_nb  <- read_xlsx("data/pre_2009/RawData_Pop2_NorthBeach_2005-08.xlsx") %>% 
              mutate(UID = paste0('NB_', 1:nrow(.)) )

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


# Population North Beach -------------------------------------------------------------------

# put number at the end of the string
ids   <- grep('[0-9]', names(pop_nb))
years <- grep('[0-9]', names(pop_nb), value=T) %>% 
            strsplit('[[:upper:]]') %>% 
            lapply( function(x) grep('[0-9]', x, value=T) ) %>% 
            unlist
vars  <- grep('[0-9]', names(pop_nb), value=T) %>% 
            strsplit('[0-9]') %>% 
            lapply( function(x) grep('[[:upper:]]', x, value=T) ) %>% 
            unlist
new_n <- paste0(vars,years)

# check the dataframe
check_df <- data.frame( new_n = strsplit(names(pop_nb)[ids],'[0-9]') %>% 
                                  unlist %>% 
                                  grep('.',.,value=T),
                        old_n = strsplit(new_n,'[0-9]') %>% 
                                  unlist %>% 
                                  grep('.',.,value=T)
                       )

# do test
expect_true( all(check_df$new_n == check_df$old_n) )

# change names
names(pop_nb)[ids] <- new_n


# size (branches)
size  <- pop_nb %>%
            # replace . with NAs, for consistency
            lapply(function(x) x = replace(x, x == '.', NA) ) %>% 
            bind_cols %>% 
            # grab size-relevant variables
            dplyr::select(UID, NEWID,BRANCHES05,BRANCHES06,
                          BRANCHES07,BRANCHES08,LENGTH08,WIDTH08 ) %>%
              mutate( LENGTH08  = as.numeric(LENGTH08),
                      WIDTH08   = as.numeric(WIDTH08) ) %>% 
              mutate( area08 = (LENGTH08/2) * (WIDTH08/2) * pi )

# stage
stage <- pop_nb %>% dplyr::select(UID,NEWID, grep("STAGE", names(.)) )

# flowers
flow  <- pop_nb %>% 
            lapply(function(x) x = replace(x, x == '.', NA) ) %>%
            bind_cols %>% 
            dplyr::select(UID,NEWID,grep("NUMRAC", names(.)) )
      
      
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
stage_long  <- stage %>% 
                  gather( "year", "stage", grep("STAGE",names(stage), value=T) ) %>%
                  # mutate( stage = dot_blank_na(stage) ) %>%
                  mutate( year  = paste0('20',gsub('[A-z]','',year)) ) %>% 
                  as.data.frame

# size long format
size_long   <- size %>% 
                  dplyr::select(-LENGTH08, -WIDTH08,-BRANCHES08) %>%
                  gather("year", "size", c('BRANCHES05', 'BRANCHES06', 'BRANCHES07',"area08") ) %>%
                  mutate( size = replace(size, size == "NA", NA) ) %>% 
                  mutate( year  = paste0('20',gsub('[A-z]','',year)) ) %>% 
                  mutate( size = as.numeric(size) ) %>% 
                  # transform n. of branches in area
                  mutate( area = coef(mod)[1] + 
                                 coef(mod)[2] * size + 
                                 coef(mod)[3] * size^2 ) %>% 
                  # size == area in 2008
                  mutate( area = replace(area, 
                                         year == 2008, 
                                         size[year == 2008]) )

# fecundity data in long form  
fecu_long   <- gather( flow, "year",  "racemes", grep('NUMRAC',names(flow),value=T) ) %>%
                  mutate( racemes = replace(racemes,
                                            racemes %in% c('NA','near EBG','111'),
                                            NA)
                         ) %>%
                  mutate( racemes = as.numeric(racemes) ) %>% 
                  mutate( year    = paste0('20',gsub('[A-Z]','',year)) )
                  

# combine all demographic data in one file
demo_pop1   <- Reduce(function(...) full_join(...), list( stage_long, 
                                                          size_long,
                                                          fecu_long ) ) %>% 
                  mutate( year  = as.numeric(year) ) %>% 
                  # homogenize stage "codes"
                  mutate( stage = replace(stage, stage == 'REP','R'),
                          stage = replace(stage, stage %in% c("Can't Connect","NA"), NA) )

# demography at time t0
demo_t0    <- demo_pop1 %>%
                setNames( c(names(.)[1:3],
                            paste0(names(.)[-c(1:3)],"_t0")) )
# demography at time t1
demo_t1    <- demo_pop1 %>%
                setNames( c(names(.)[1:3],
                            paste0(names(.)[-c(1:3)],"_t1")) ) %>%
                mutate( year = year - 1)

# transition data frame 
trans_pop_nb <- inner_join(demo_t0, demo_t1) %>%
                  rename( newid = NEWID ) %>% 
                  # introduce survival at time t1
                  mutate( surv_t1 = 1 ) %>%
                  mutate( surv_t1 = replace(surv_t1, stage_t1 == "D", 0)  ) %>%
                  mutate( surv_t1 = replace(surv_t1, is.na(stage_t0),  NA) ) %>%
                  mutate( surv_t1 = replace(surv_t1, stage_t1 == "SL", NA) ) %>%
                  mutate( surv_t1 = replace(surv_t1, is.na(stage_t1),  NA) ) %>%
                  # introduce flowering at time t0 and t1
                  mutate( flow_t0 = 0,
                          flow_t1 = 0) %>% 
                  mutate( flow_t0 = replace(flow_t0, stage_t0 == "R", 1),
                          flow_t1 = replace(flow_t1, stage_t1 == "R", 1) ) %>%
                  # introduce NAs: everything that is NOT "Rep" or "Non",
                  mutate( flow_t0 = replace(flow_t0, 
                                            stage_t0 == "SL" | stage_t0 == "D" | is.na(stage_t0), 
                                            NA),
                          flow_t1 = replace(flow_t1, 
                                            stage_t1 == "SL" | stage_t1 == "D" | is.na(stage_t1),
                                            NA) ) %>%
                  # NAs whenever is.na(stage_t) == TRUE
                  mutate( flow_t0 = replace(flow_t0, is.na(stage_t0), NA) ) %>%
                  mutate( flow_t1 = replace(flow_t1, is.na(stage_t1), NA) )
                  

# Population 1 -----------------------------------------------------------------------------

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
            rename( X.racemes.2007 = racemes,
                    X.racemes.2008 = racemes.1 )


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
stage_long  <- long_data(stage, "stage", Filter(function(x) grepl("stage.",x), names(stage)) ) %>%
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
fecu_long   <- long_data( flow,  "racemes", paste0("X.racemes.", 2005:2008) ) %>%
                  mutate( racemes = dot_blank_na(racemes) ) %>%
                  mutate( year    = year_form(year),
                          racemes = as.numeric(racemes) )

# combine all demographic data in one file
demo_pop1   <- Reduce(function(...) full_join(...), list( stage_long, 
                                                          size_long,
                                                          fecu_long ) )
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
trans_pop1 <- inner_join(demo_t0, demo_t1) %>%
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
                                          stage_t0 == "Dead" | is.na(stage_t0) |
                                          stage_t0 == "SL",
                                          NA),
                        flow_t1 = replace(flow_t1, 
                                          stage_t1 == "Dead" | is.na(stage_t1) |
                                          stage_t1 == "SL" | stage_t1 == "not found", 
                                          NA) ) %>%
                # NAs whenever is.na(stage_t) == TRUE
                mutate( flow_t0 = replace(flow_t0, is.na(stage_t0), NA) ) %>%
                mutate( flow_t1 = replace(flow_t1, is.na(stage_t1), NA) ) #%>% 
                # remove all cases in which stage_t0 and t1 are NA (NO DATA provided in this cases)
                subset( !(is.na(stage_t0) & is.na(stage_t1)) )


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
plot(area_t1 ~ area_t0, data=subset(trans_pop1)
# plot racemes
plot(racemes_t1 ~ area_t1, data = trans_pop1)
plot(racemes_t0 ~ area_t0, data = trans_pop1)
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
plot(racemes_t0 ~ area_t0, data = subset(trans_pop1) )
# plot survival 
plot( jitter(surv_t1) ~ area_t0, data=subset(trans_pop1,stage_t0 != "SL"))
# plot flowering
plot( jitter(flow_t0) ~ area_t0, data=subset(trans_pop1))

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
                  rename( X.fl.stalks.2007 = flower.stalks,
                          X.fl.stalks.2008 = racemes.1 ) %>%
                  long_data("racemes", paste0("X.fl.stalks.", 2005:2008) ) %>%
                  mutate( racemes = dot_blank_na(racemes) ) %>%
                  mutate( year    = year_form(year) )


# combine all demographic data in one file
demo_pop8   <- Reduce( function(...) full_join(...), list( stage8_long, 
                                                           size8_long,
                                                           fecu8_long ) ) 
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
trans_pop8 <- inner_join(demo8_t0, demo8_t1) %>%
                rename( plant_n = plant..) %>%
                rename( newid = new.tag.number ) %>% 
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
                                          is.na(stage_t0) & stage_t1 == "Non", 
                                          NA) ) %>%
                # introduce flowering at time t0 and t1
                mutate( flow_t0 = 0,
                        flow_t1 = 0) %>% 
                mutate( flow_t0 = replace(flow_t0, stage_t0 == "Rep", 1),
                        flow_t1 = replace(flow_t1, stage_t1 == "Rep", 1) ) %>%
                # introduce NAs: everything that is NOT "Rep" or "Non",
                mutate( flow_t0 = replace(flow_t0, 
                                          stage_t0 == "Dead" | is.na(stage_t0) |
                                          stage_t0 == "SL",
                                          NA),
                        flow_t1 = replace(flow_t1, 
                                          stage_t1 == "Dead" | is.na(stage_t1) |
                                          stage_t1 == "SL" | stage_t1 == "not found", 
                                          NA) ) %>%
                # NAs whenever is.na(stage_t) == TRUE
                mutate( flow_t0 = replace(flow_t0, is.na(stage_t0), NA) ) %>%
                mutate( flow_t1 = replace(flow_t1, is.na(stage_t1), NA) ) %>% 
                # remove all cases in which stage_t0 and t1 are NA (NO DATA provided in this cases)
                subset( !(is.na(stage_t0) & is.na(stage_t1)) ) %>% 
                # remove all cases in which stage_t0 and t1 are 'Dead' (NO DATA provided in this cases)
                subset( !(stage_t0 == 'Dead' & stage_t1 == 'Dead' & 
                          !is.na(stage_t0) & !is.na(stage_t1) ) )


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


# population 1, 8, and NB together -----------------------------------------------  
  
names(trans_pop1)
names(trans_pop8)

# create plant id, convet population to character
trans_pop1   <- trans_pop1 %>%
                  mutate( plant_id   = paste( UID, location, sep="_"),
                          population = as.character(population),
                          UID        = as.character(UID) )
trans_pop8   <- trans_pop8 %>%
                  mutate( plant_id = paste( UID, plant_n, sep="_"),
                          population = 'ATT (8)',
                          UID      = as.character(UID) )
trans_pop_nb <- trans_pop_nb %>%
                  mutate( plant_id = paste( UID, newid, sep="_"),
                          population = 'NB (2)',
                          newid    = as.numeric(newid) )

# put it all together and send store
trans_0508 <- bind_rows(trans_pop1, trans_pop8, trans_pop_nb) 

par(mfrow=c(2,2))
plot(jitter(surv_t1) ~ area_t0, data = trans_0508)
plot(area_t1    ~ area_t0, data = trans_0508)
plot(jitter(flow_t1) ~ area_t0, data = trans_0508)
plot(racemes_t1 ~ area_t0, data = trans_0508)


# write out -------------------------------------------------

# homogenize 
trans_0508_out <- trans_0508 %>% 
                    rename( numrac_t1  = racemes_t1, 
                            numrac_t0  = racemes_t0,
                            # "location" here is a subsample of two populations
                            sub_sample = location ) %>% 
                    # remove rows where stage_t0 and stage_t1 are both NA!
                    subset( !(is.na(stage_t0) & is.na(stage_t1)) ) %>% 
                    # rename location accordingly
                    rename( location = population ) %>% 
                    # pop 1 is "AL (1)"
                    mutate( location = replace(location, 
                                               location == '1',
                                               'AL (1)') ) %>% 
                    # label year as "transition year"
                    mutate( transition = paste(year, year+1, sep = '-') ) %>% 
                    # homogenize stage names
                    mutate( stage_t0 = replace( stage_t0, stage_t0 == 'Dead', 'D'),
                            stage_t0 = replace( stage_t0, stage_t0 == 'Rep', 'R'),
                            stage_t0 = replace( stage_t0, stage_t0 == 'Non', 'NR'),
                            stage_t1 = replace( stage_t1, stage_t1 == 'Dead', 'D'),
                            stage_t1 = replace( stage_t1, stage_t1 == 'Rep', 'R'),
                            stage_t1 = replace( stage_t1, stage_t1 == 'Non', 'NR'),
                            stage_t1 = replace( stage_t1, stage_t1 == 'not found', NA) ) %>% 
                    # create NEW ID!
                    mutate( newid = replace(newid, is.na(newid), 
                                            paste0('PRE_08_',1:nrow(subset(., is.na(newid))))) ) %>% 
                    # useless columns used in "merges"
                    select(-UID, -plant_id, -plant_n)

# store it
write.csv(trans_0508_out, "data/lupine_05_08.csv", row.names=F)
