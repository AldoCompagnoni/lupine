setwd("C:/cloud/Dropbox/lupine")
library(tidyverse)
library(testthat)

dat_17  <- read.csv("data/from_2008_to_2017/LUTI_0817_AllPops_July212017.csv", stringsAsFactors = F)
dat_mar <- read.csv('data/from_2008_to_2017/LUTI_0817_AllPopsMergedMarch2018.csv', stringsAsFactors = F)
meta    <- read.csv("data/from_2008_to_2017/LUTI_0817_AllPops_METADATA.csv", stringsAsFactors = F)

# try to update 
dat     <- dat_mar %>% 
              dplyr::select(Y12LEN,	Y12WID, FID, NEWID) %>% 
              subset( FID != "TC\"") %>% 
              subset( FID != "dead in 2014 but TC in 2015\"" ) %>% 
              subset( NEWID != "NO" ) %>% 
              subset( NEWID != "D" ) %>% 
              mutate( FID   = as.numeric(FID),
                      NEWID = as.numeric(NEWID) ) %>% 
              right_join(dat_17) %>% 
              mutate( Y12LEN = replace(Y12LEN, is.na(Y12LEN), 0),
                      Y12WID = replace(Y12WID, is.na(Y12WID), 0) )

# not needed columns
d_need <- dat %>%
            dplyr::select( -c(TAGYEAR:FLWRLEN) ) %>% # not useful
            dplyr::select(-grep("[0-9][0-9]HT", names(.)) ) %>%  # remove HT info - only for LUCH
            dplyr::select(-grep("SQUOOP|SQPCOMM|CHCK|CHK", names(.)) ) %>%
            dplyr::select(-grep("SUCCSTG", names(.)) ) %>%
            dplyr::select(-grep("CHCK", names(.)) ) %>%
            subset( SPECIES == "LUTI" ) %>% # retai nonly LUTI individuals
            dplyr::select( -c(GNTPYR, X11numrac) )

# Separate data frame per each variable 
by_variable <- function(var){
  
  d_need %>%
    dplyr::select( c(FID:LOCATION, grep(paste0("Y[0-9][0-9]",var), names(.))) ) %>%
    gather( year, stage, -c(FID, NEWID, LOCATION) ) %>%
    mutate( year = gsub("Y", "", year) ) %>%
    mutate( year = gsub(var, "", year) ) %>%
    mutate( year = as.numeric(paste0("20", year)) ) %>%
    mutate( id = paste(FID, NEWID, sep = "_") ) %>%
    rename_( .dots = setNames("stage", var)) 
    
}

# individual level measures
indiv_l <- lapply( c("STAGE","LEN","WID","NUMAB","NUMCL","NUMINT"), by_variable) # c("WID","NUMAB","NUMCL","NUMINT")

# too many lines?
indiv_d <- Reduce( function(...) full_join(...), indiv_l)


# introduce NAs ---------------------------------------------------------------------------

# specific variable where to introduce NAs
na_introduce <- function(var, vec){
  
  indiv_d[,var] <- gsub(paste(vec, collapse="|"), NA, indiv_d[,var] )
  return(indiv_d)
  
}

# STAGE (all those stages are supposed to be NA)
indiv_d <- na_introduce("STAGE", c("^999$","^999/Elim$","^BLNOVR$",
                                   "^NF/BLNOVR$","^NF/E$", "R/Elim") )

# LEN and WID: Introduce NA in place of 9999 and 99
indiv_d <- na_introduce("LEN","9999|999")
indiv_d <- na_introduce("WID","9999|999")

# NUMCL, NUMAB, NUMINT: introuce NA where 9999 or 999
indiv_d <- na_introduce("NUMCL","9999|999")
indiv_d <- na_introduce("NUMAB","9999|999")
indiv_d <- na_introduce("NUMINT","9999|999")

# make raceme numbers "numeric"
indiv_d <- indiv_d %>%
              mutate( NUMAB  = as.numeric(NUMAB),
                      NUMCL  = as.numeric(NUMCL),
                      NUMINT = as.numeric(NUMINT) ) %>%
              subset( !(is.na(STAGE) & is.na(LEN) & is.na(WID) & is.na(NUMAB) & is.na(NUMCL) & is.na(NUMINT)) )

# transitions ----------------------------------------------------------------

# vital rates labels
vr      <- c("STAGE","LEN", "WID", "NUMAB", "NUMCL", "NUMINT")
vr_t0   <- paste0(vr, "_t0")
vr_t1   <- paste0(vr, "_t1")

# order and rename variables
indiv_t0  <- indiv_d %>%
                dplyr::select(id, FID, NEWID, LOCATION, 
                              year, STAGE, LEN, WID, 
                              NUMAB, NUMCL, NUMINT) %>%
                setNames( c("id", "FID", "NEWID", "LOCATION", "year",
                            vr_t0) )

# year t+1
indiv_t1  <- indiv_t0 %>%
                mutate( year = year - 1 ) %>%
                setNames( c("id", "FID", "NEWID", "LOCATION", "year",
                            vr_t1) )

# transition data
trans_d   <- left_join(indiv_t0, indiv_t1)


# len_t0 here is ==2, not 7836 ---------------------------------------------------------
subset(trans_d, id  == "4759_7836")

# len_t0 here is ==2, not 7836
trans_d <- trans_d %>%
              mutate( LEN_t1 = replace(LEN_t1, which(LEN_t1 == "7836"), 2) ) %>%
              mutate( LEN_t0 = replace(LEN_t0, which(LEN_t0 == "7836"), 2) )


# remove "99" in LEN,WID, NUMCL, NUMAB, NUMINT ------------------------------------------
# when "99" corresponds to SL or NR.
# only one R individual that actually has 99 racemes!

# function to remove instances of "99"
remove_99 <- function(var){
  
  var_t0  <- paste0(var,"_t0")
  var_t1  <- paste0(var,"_t1")
  st_t0   <- "STAGE_t0"
  st_t1   <- "STAGE_t1"
  
  # Seedling instances
  trans_d[,var_t0]  <- replace(trans_d[,var_t0],
                                which(trans_d[,st_t0] == "SL" & trans_d[,var_t0] == "99"),
                                NA)
  trans_d[,var_t1]  <- replace(trans_d[,var_t1],
                               which(trans_d[,st_t1] == "SL" & trans_d[,var_t1] == "99"),
                               NA)
  
  # Non-reproductive individuals instances
  trans_d[,var_t0]  <- replace(trans_d[,var_t0],
                               which(trans_d[,st_t0] == "NR" & trans_d[,var_t0] == "99"),
                               NA)
  trans_d[,var_t1]  <- replace(trans_d[,var_t1],
                               which(trans_d[,st_t1] == "NR" & trans_d[,var_t1] == "99"),
                               NA)
  return(trans_d)
  
}

#LEN,WID, NUMCL, NUMAB, NUMINT
trans_d <- remove_99("LEN")
trans_d <- remove_99("WID")
trans_d <- remove_99("NUMCL")
trans_d <- remove_99("NUMAB")
trans_d <- remove_99("NUMINT")

# a few tests
expect_equal( nrow(subset(trans_d, STAGE_t0 == "SL" & LEN_t0 =="99")), 0)
expect_equal( nrow(subset(trans_d, STAGE_t0 == "NR" & WID_t0 =="99")), 0)
expect_equal( nrow(subset(trans_d, STAGE_t1 == "NR" & NUMAB_t1 =="99")), 0)
expect_equal( nrow(subset(trans_d, STAGE_t0 == "NR" & NUMINT_t0 =="99")), 0)


# check weird instances ------------------------------------------------------------
trans_d %>% 
    dplyr::select(STAGE_t0, STAGE_t1) %>% 
    unique %>% 
    arrange(STAGE_t0, STAGE_t1)

trans_d %>% 
    dplyr::select(STAGE_t0, STAGE_t1) %>% 
    unique %>%
    subset( STAGE_t0 != "Elim" & STAGE_t0 != "NF" & STAGE_t0 != "TAGLATER") %>%
    arrange( STAGE_t0, STAGE_t1)


# ISSUES ---------------------------------------------------------

# Elim: 
# 1. why two Elim_t0 is still Elim_t1?. 
#    It means that they didn't eliminate the indiv!
#    Elim does not have data in the current year.    
# 2. Three instances Elim_t0 is D_t1. Cannot use it in year t1. it's going to be NA
trans_d %>% 
  dplyr::select(STAGE_t0, STAGE_t1, FID, NEWID) %>% 
  unique %>%
  subset( STAGE_t0 == "Elim") %>%
  arrange( STAGE_t0, STAGE_t1) %>%
  group_by( STAGE_t0, STAGE_t1 ) %>%
  summarise( rep = n())

# Implement # 2
trans_d <- mutate(trans_d,
                  STAGE_t1 = replace(STAGE_t1, STAGE_t0 == "Elim" & STAGE_t1 == "D", NA) )


# NF 1: individual, becomes NR, but no size information
# this plant won't give you information. So we have no transitions because of too many gaps
# Will leave this as is - no idiosync
trans_d %>% 
  dplyr::select(STAGE_t0, STAGE_t1, WID_t1, LEN_t1, id, FID, NEWID, year, LOCATION) %>% 
  unique %>%
  subset( STAGE_t0 == "NF" & STAGE_t1 == "NR") %>%
  arrange( STAGE_t0, STAGE_t1)
# Will remove this - I don't know how to use this information (it's not size-structured)
trans_d <- mutate(trans_d,
                  STAGE_t1 = replace(STAGE_t1, id == "745_4663" & STAGE_t0 == "NF", NA) )


# NF 1: NF becomes "A" without size information
# They pulled the tags the following year, because NOT IN INTERESTING LOCATION. 
trans_d %>% 
  dplyr::select(STAGE_t0, STAGE_t1, WID_t1, LEN_t1, FID, NEWID, year) %>% 
  unique %>%
  subset( STAGE_t0 == "NF" & STAGE_t1 == "A") %>%
  arrange( STAGE_t0, STAGE_t1)


# TAGLATER: TAGLATER_t0 to TAGLATER_t1. VERY COMMON
# OF COURSE, it's tagged LATER, anytime later
# also: TAGLATER_t0 to NA_t1.
trans_d %>%
  dplyr::select(STAGE_t0, STAGE_t1, FID, NEWID) %>% 
  unique %>%
  subset( STAGE_t0 == "TAGLATER") %>%
  arrange( STAGE_t0, STAGE_t1 ) %>%
  group_by( STAGE_t0, STAGE_t1 ) %>%
  summarise( rep = n() )

# NF: found in year t, but found in year t+1?
# no data for that transition. These are problem plants
trans_d %>% 
  dplyr::select(STAGE_t0, STAGE_t1, FID, NEWID) %>% 
  unique %>%
  subset( STAGE_t0 == "NF/E") %>%
  arrange( STAGE_t0, STAGE_t1) %>%
  group_by( STAGE_t0, STAGE_t1 ) %>%
  summarise( rep = n() )


# D: from D_t0 to R_t1 (6 cases), NR_t1 (3), Counting (4) ------------------
# dead to counting is an NA 
# 
trans_d %>% 
  dplyr::select(STAGE_t0, STAGE_t1, FID, NEWID) %>% 
  unique %>%
  subset( STAGE_t0 == "D") %>%
  arrange( STAGE_t0, STAGE_t1) %>%
  group_by( STAGE_t0, STAGE_t1 ) %>%
  summarise( rep = n() )

# NEWID 932 ELIMINATE. Not confident the plant is dormant
# NEWID 4441 eliminate, don't use it.
# NEWID 4469: a mistake!
subset(trans_d, STAGE_t0 == "D" & STAGE_t1 == "NR")
# Implement newid 932. I will eliminate this NEWID after 2010
trans_d <- subset(trans_d, !(NEWID == 932 & year > 2010) )
trans_d <- subset(trans_d, NEWID != 4441 )
trans_d <- subset(trans_d, NEWID != 4469 )


# NEWID 885: remove this plant
# NEWID 887: THIS IS A DORMANT in 2010!!!!
# NEWID 1129: Don't use it after 2012
# NEWID 2789: Don't use it! MORE LIKELY FOR THIS TO BE A NEW PLANT. 
# NEWID 4095: DORMANT IN 2012. This is a clear case
# NEWID 3534: DORMANT IN 2012. This is a clear case
subset(trans_d, STAGE_t0 == "D" & STAGE_t1 == "R")
trans_d <- subset(trans_d, NEWID != 885)
trans_d <- subset(trans_d, !(NEWID == 1129 & year > 2012) )
trans_d <- subset(trans_d, NEWID != 2789 )
# actually dormant individuals
trans_d <- mutate(trans_d,
                  STAGE_t0 = replace(STAGE_t0, NEWID==887 & year == 2010, "DORM"),
                  STAGE_t1 = replace(STAGE_t1, NEWID==887 & year == 2009, "DORM") )
trans_d <- mutate(trans_d,
                  STAGE_t0 = replace(STAGE_t0, NEWID==4095 & year == 2012, "DORM"),
                  STAGE_t1 = replace(STAGE_t1, NEWID==4095 & year == 2011, "DORM") )
trans_d <- mutate(trans_d,
                  STAGE_t0 = replace(STAGE_t0, NEWID==3534 & year == 2012, "DORM"),
                  STAGE_t1 = replace(STAGE_t1, NEWID==3534 & year == 2011, "DORM") )

# remove D_t0 to Counting_t0
trans_d <- subset(trans_d, !(STAGE_t0 == "D" & STAGE_t1 == "Counting") )

# Remove "taglater", and "dead" and "eliminate" in stage t0, taglater in stage t1
trans_d <- trans_d %>%
              subset( !(STAGE_t0 %in% c("TAGLATER", "D", "Elim","Counting","A")) ) %>%
              subset( !(STAGE_t1 %in% c("TAGLATER", "Elim","Counting","A")) ) 

# SIZE -------------------------------------------------------------------------------
trans_d <- trans_d %>%
              mutate( LEN_t0  =  as.numeric(LEN_t0),
                      LEN_t1  =  as.numeric(LEN_t1),
                      WID_t0  =  as.numeric(WID_t0),
                      WID_t1  =  as.numeric(WID_t1)  ) %>%
              mutate( area_t0 = LEN_t0 * WID_t0 * pi,
                      area_t1 = LEN_t1 * WID_t1 * pi ) %>%
              mutate( area_t0 = replace(area_t0, 
                                        STAGE_t0 %in% c("NF","DORM"), NA),
                      area_t1 = replace(area_t1, 
                                        STAGE_t1 %in% c("NF","DORM","D"), NA),
                      area_t0 = replace(area_t0, is.na(STAGE_t0),  NA),
                      area_t1 = replace(area_t1, is.na(STAGE_t1),  NA) )


# SURVIVAL ---------------------------------------------------------------------------

# Survival is NA UNLESS individuals have "D", "NR", "R", or "DORM" status
trans_d <- trans_d %>%
              mutate( surv_t1 = NA ) %>%
              mutate( surv_t1 = replace(surv_t1, STAGE_t1 == "D", 0) ) %>%
              mutate( surv_t1 = replace(surv_t1, STAGE_t1 %in% c("NR","R","DORM"), 1) ) %>%
              # NAs when no reliable information at STAGE_t0
              mutate( surv_t1 = replace(surv_t1, STAGE_t0 == "NF", NA) ) %>%
              mutate( surv_t1 = replace(surv_t1, is.na(STAGE_t0), NA) ) %>%
              # remove a handful of nonsensical transitions (probable mistakes)
              subset( !(is.na(STAGE_t0) & STAGE_t1 == "D") )


# FLOWERING --------------------------------------------------------------------------

# flowering at time t0 and t1 is NA UNLESS individuals have "R" and "NR" information 
trans_d <- trans_d %>%
              mutate( flow_t0 = NA, 
                      flow_t1 = NA ) %>%
              mutate( flow_t0 = replace(flow_t0, STAGE_t0 == "R",  1),
                      flow_t1 = replace(flow_t1, STAGE_t1 == "R",  1) ) %>%
              mutate( flow_t0 = replace(flow_t0, STAGE_t0 == "NR", 0),
                      flow_t1 = replace(flow_t1, STAGE_t1 == "NR", 0) ) 


# FERTILITY ------------------------------------------------------------------------
trans_d <- trans_d %>%
              rowwise() %>%
              mutate( NUMALL_t0 = NUMAB_t0 + NUMCL_t0 + NUMINT_t0 ) %>%
              rowwise %>%
              mutate( NUMALL_t1 = NUMAB_t1 + NUMCL_t1 + NUMINT_t1 ) %>%
              # last step for clarity
              arrange( LOCATION, year) 

# WRITE DOWN PROBLEM INDIVIDUALS -----------------------------------------------

fert0_t1 <- trans_d %>% 
              subset( STAGE_t1 == 'R' & NUMALL_t1 == 0 ) %>%
              select(year,LOCATION,NEWID) %>% 
              mutate( year = year + 1 ) %>% 
              unique
  
fert0_t0 <- trans_d %>% 
              subset( STAGE_t0 == 'R' & NUMALL_t0 == 0 ) %>%
              select(year,LOCATION,NEWID) %>% 
              unique

# store it
fert0_out <- rbind(fert0_t1, fert0_t0) %>% unique
write.csv(fert0_out, 'fert0_ids.csv', row.names=F)


# Look at raw data --------------------------------------------------
require(rgdal)
library(dplyr)
library(tidyr)
library(testthat)

# list of shape files (one for each population)
read_f    <- paste0('LUTI_0818_', c('AL', 'ATT', 'BR', 'BS', 'DR', 'NB', 'POP9'))

# read all shape files
shape_df  <- lapply(read_f, function(x) readOGR(dsn = "data/shapefiles", layer = x) ) %>% 
                # extract data frames from the shape files
                lapply(function(x) attributes(x)$data ) %>% 
                # bind it all together
                bind_rows      

# not needed columns
d_interim <- shape_df %>%
              # columns not useful
              dplyr::select( -c(GPS_DATE:FLWRLEN) ) %>% 
              # remove HT info: only for LUCH, not LUTI
              dplyr::select(-grep("[0-9][0-9]HT", names(.)) ) %>%  
              # discard other columns
              dplyr::select(-grep("SQUOOP|SQPCOMM", names(.)) ) %>%
              dplyr::select(-grep("CHCK|CHK|CHECK", names(.)) ) %>%          
              # remove successional information
              dplyr::select(-grep("SUCC|NUMFR|NEAR_|Y15HAB", names(.)) ) %>%
              dplyr::select( -c(GNTPYR, GENOTYPE, X11numrac,MATINGYST,POST2012) ) %>% 
              # No number of branches (this data + 2018 data is in data/lupine_area_vs_branches.csv
              dplyr::select( -c(Y08_NUMBRN, Y08_BRNCHS) ) %>% 
              # retai nonly LUTI individuals
              subset( SPECIES == "LUTI" ) %>% 
              # add location information at a finer scale ("sub_location")
              # for now we only split ATT in two sub-populations
              mutate( sub_loc = ATTSubpop ) %>% 
              mutate( sub_loc = replace(sub_loc, is.na(sub_loc), LOCATION[is.na(sub_loc)]) ) %>% 
              # link subpopulation names directly to ATT (8)
              mutate( sub_loc = replace(sub_loc, sub_loc=='beach',   'ATT (8)_beach')) %>% 
              mutate( sub_loc = replace(sub_loc, sub_loc=='blowout', 'ATT (8)_blowout')) %>% 
  
              # Individuals to delete/change (according to Eleanor's indications)
  
              # NEWID==99999 from BR, NOT from 2015, should be dropped!
              subset( !(NEWID %in% 99999 & LOCATION %in% "BR (6)") ) %>% 
              # NEWID %in% c(1053, 1055, 1056, 1057): delete IF FROM DR.
              subset( !(NEWID %in% 1053  & LOCATION %in% 'DR (3)') ) %>% 
              subset( !(NEWID %in% 1055  & LOCATION %in% 'DR (3)') ) %>% 
              subset( !(NEWID %in% 1056  & LOCATION %in% 'DR (3)') ) %>% 
              subset( !(NEWID %in% 1057  & LOCATION %in% 'DR (3)') ) %>% 
              # THESE PLANTS ARE DODGY, DO NOT TRUST DATA
              subset( !(NEWID %in% c(885,887)  & LOCATION %in% 'ATT (8)')) %>% 
              # NEWID == 4483: delete IF FROM ATT (8)
              subset( !(NEWID %in% 4483 & LOCATION %in% 'ATT (8)') ) %>% 
              # NEWID %in% 708
              subset( !(NEWID %in% 708) ) %>% 
              # change 10074 at ATT (8) to 10071: it's was a typo
              mutate( NEWID = replace(NEWID, NEWID == 10074 & LOCATION == 'ATT (8)', 10071))
              

# Introduce NEWID when needed
d_need <- d_interim %>% 
              # NEWID==0 to SL_08_1, SL_08_2, ect. IF SEEDLINGS            
              mutate( NEWID = replace(NEWID, NEWID == 0 & Y08STAGE == 'SL', 
                                      paste0('SL_08_',1:nrow(subset(.,NEWID == 0 &
                                                                      Y08STAGE == 'SL')))) ) %>% 
              # NEWID==0 to NR_08_1, NR_08_2, ect. IF Non reproductive
              mutate( NEWID = replace(NEWID, NEWID == 0 & Y08STAGE == 'NR', 
                                      paste0('NR_08_',1:nrow(subset(.,NEWID == 0 &
                                                                      Y08STAGE == 'NR')))
                                            ) 
                     ) %>% 
              # NEWID==99999 to FERT_15_1, FERT_15_2, ect. (ONLY FERTILITY DATA for these individuals)
              mutate( NEWID = replace(NEWID, NEWID == 99999,
                                      paste0('FERT_15_',1:nrow(subset(.,NEWID == 99999))) 
                                      ) 
                    )


d_need %>% 
  subset( NEWID    == fert0_out$NEWID[43] & 
          LOCATION == fert0_out$LOCATION[43])
fert0_out[43,]
