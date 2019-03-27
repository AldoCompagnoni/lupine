# Script to automatically format demographic data starting from the shapefile
# 1. read original files and pre-format 
# 2. 
setwd("C:/cloud/Dropbox/lupine")
options(stringsAsFactors = F)
require(rgdal)
library(dplyr)
library(tidyr)
library(testthat)


# 1. read original files and pre-format -------------------------------------------------------------

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
              # retain only LUTI individuals
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


# tests:

# dataset need have NO DOUBLES
expect_equal( nrow(d_need), unique(d_need) %>% nrow)

# No NEWID should be repeated more than once
expect_equal(d_need %>% count(NEWID) %>% subset(n > 1) %>% nrow, 0)


# Separate data frame per each variable 
by_variable <- function(var){
  
  d_need %>% 
    dplyr::select( c(NEWID,TAGYEAR,LOCATION, sub_loc, grep(paste0("Y[0-9][0-9]",var), names(.))) ) %>%
    gather( year, stage, -c(NEWID,TAGYEAR,LOCATION, sub_loc) ) %>%
    mutate( year = gsub("Y", "", year) ) %>%
    mutate( year = gsub(var, "", year) ) %>%
    mutate( year = as.numeric(paste0("20", year)) ) %>%
    rename_( .dots = setNames("stage", var)) 
    
}

# individual level measures
indiv_l <- lapply( c('STAGE','LEN','WID','NUMRAC','NUMAB','NUMCL','NUMINT'), by_variable) # c('WID','NUMAB','NUMCL','NUMINT')

# too many lines?
indiv_d <- Reduce( function(...) full_join(...), indiv_l)

# test: the indiv_d need have as many rows as the biggest dataframe in 
expect_equal( nrow(indiv_d), max(sapply(indiv_l,nrow)) )


# 3. introduce NAs ---------------------------------------------------------------------------

# specific variable where to introduce NAs
na_introduce <- function(var, vec){
  
  indiv_d[,var] <- gsub(paste(vec, collapse="|"), NA, indiv_d[,var] )
  return(indiv_d)
  
}

# STAGE (all these stages are supposed to be NA)
indiv_d <- na_introduce("STAGE", c("^999$","^999/Elim$","^BLNOVR$",
                                   "^NF/BLNOVR$","^NF/E$", "R/Elim",
                                   'N/A', 'D/Elim', 'Elim', 'NotRec',
                                   'NF', 'NF/Elim', 'NF/Bur', 'B/999', 'NF/999', 'NF/B/E',
                                   '[0-9]{4}','Counting','NA') )

# LEN and WID: Introduce NA in place of 9999 and 99
indiv_d <- na_introduce("LEN","9999|999")
indiv_d <- na_introduce("WID","9999|999")

# NUMCL, NUMAB, NUMINT: introuce NA where 9999 or 999
indiv_d <- na_introduce("NUMRAC","9999|999|1111|111")
indiv_d <- na_introduce("NUMCL","9999|999")
indiv_d <- na_introduce("NUMAB","9999|999")
indiv_d <- na_introduce("NUMINT","9999|999")

# make raceme numbers "numeric"
indiv_d <- indiv_d %>%
              mutate( NUMAB  = as.numeric(NUMAB),
                      NUMCL  = as.numeric(NUMCL),
                      NUMINT = as.numeric(NUMINT) ) %>%
              # remove useless data rows (all NAs)
              subset( !(is.na(STAGE) & is.na(LEN) & is.na(WID) & is.na(NUMAB) & is.na(NUMCL) & is.na(NUMINT)) )

# remove "99" in LEN, WID, NUMRAC, NUMCL, NUMAB, NUMINT when Stage == "SL" or "NR" ---

# function to remove instances of "99"
remove_99 <- function(var_in){
  
  # Seedling instances
  indiv_d[,var_in]  <- replace(indiv_d[,var_in],
                                which(indiv_d[,'STAGE'] == "SL" & indiv_d[,var_in] == "99"),
                                NA)

  # Non-reproductive individuals instances
  indiv_d[,var_in]  <- replace(indiv_d[,var_in],
                                which(indiv_d[,'STAGE'] == "NR" & indiv_d[,var_in] == "99"),
                                NA)
  
  return(indiv_d)
  
}

#LEN,WID, NUMCL, NUMAB, NUMINT
indiv_d <- remove_99("LEN")
indiv_d <- remove_99("WID")
indiv_d <- remove_99("NUMRAC")
indiv_d <- remove_99("NUMCL")
indiv_d <- remove_99("NUMAB")
indiv_d <- remove_99("NUMINT")

# create "NOTABORTED" column
indiv_d <- indiv_d %>% mutate( NOTAB = NUMCL + NUMINT )


# transitions ------------------------------------------------------------------------

# vital rates labels
vr      <- c("STAGE","LEN", "WID", 'NUMRAC', "NUMAB", "NUMCL", "NUMINT", "NOTAB")
vr_t0   <- paste0(vr, "_t0")
vr_t1   <- paste0(vr, "_t1")

# order and rename variables
indiv_t0  <- indiv_d %>%
                dplyr::select(NEWID, TAGYEAR, LOCATION, sub_loc, year, 
                              STAGE, LEN, WID, 
                              NUMRAC, NUMAB, NUMCL, NUMINT, NOTAB) %>%
                setNames( c('NEWID', 'TAGYEAR', 'LOCATION', 'sub_loc', 'year',
                            vr_t0) )

# year t+1
indiv_t1  <- indiv_t0 %>%
                mutate( year = year - 1 ) %>%
                setNames( c('NEWID', 'TAGYEAR','LOCATION', 'sub_loc', 'year',
                            vr_t1) )

# transition data frame
trans_d   <- full_join(indiv_t0, indiv_t1)


# remove "99" in LEN, WID, NUMRAC, NUMCL, NUMAB, NUMINT when Stage_t0/t1 == "SL" or "NR" ---

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
trans_d <- remove_99("NUMRAC")
trans_d <- remove_99("NUMCL")
trans_d <- remove_99("NUMAB")
trans_d <- remove_99("NUMINT")


# test that there are no 99
expect_equal( nrow(subset(trans_d, STAGE_t0 == "SL" & LEN_t0 =="99")), 0)
expect_equal( nrow(subset(trans_d, STAGE_t0 == "NR" & WID_t0 =="99")), 0)
expect_equal( nrow(subset(trans_d, STAGE_t1 == "NR" & NUMAB_t1 =="99")), 0)
expect_equal( nrow(subset(trans_d, STAGE_t0 == "NR" & NUMINT_t0 =="99")), 0)


# Remove or fix specific 'WEIRD' instances -----------------------------------------------
trans_d <- trans_d %>% 
              # 'A' in Stage_t0 is an NA
              mutate( STAGE_t0 = replace(STAGE_t0, STAGE_t0 == 'A', NA) ) %>% 
              # Also, data with Stage_t0 == 'NA' to Stage_t1 == 'A' is useless
              subset( !(is.na(STAGE_t0) & STAGE_t1 %in% 'A') ) %>%  
              # remove when Stage_t0 and Stage_t1 are both NAs: 
              # most these involve size information from individuals that have died (!!!)
              # Got to this conclusion by checking NEWID %in% c(885, 887, 425, 4650, 8388)
              subset( !(is.na(STAGE_t0) & is.na(STAGE_t1)) ) %>% 
              # remove when Stage_t0 is NA and Stage_t1 is == 'D'
              # Got to this conclusion by checking NEWID %in% c(3719, 3243, 3267, 8980, 8449)
              subset( !(is.na(STAGE_t0) & STAGE_t1 %in% 'D') ) %>% 
              # remove transitions from "D" to "D" (dead to dead)
              # these simply come because individuals were flagged "dead" for years after death 
              subset( !(STAGE_t0 %in% 'D' & STAGE_t1 %in% 'D') ) %>% 
              # remove dead/alive individuals at stage t0 (useless to us)
              subset( !(STAGE_t0 %in% c("D", "A")) ) #@%>%
              
              #rename( newid = NEWID )

# no instances of this should occurs
expect_equal(subset(trans_d, is.na(STAGE_t0) & STAGE_t1 == 'D') %>% nrow, 0)

              
# NEWID 932 ELIMINATE. Not confident the plant is dormant
# NEWID 4441 eliminate, don't use it.
# NEWID 4469: a mistake!
subset(trans_d, STAGE_t0 == "D" & STAGE_t1 == "NR")
# Implement newid 932. I will eliminate this NEWID after 2010
trans_d <- subset(trans_d, !(NEWID %in% 932 & year > 2010) )
trans_d <- subset(trans_d, !(NEWID %in% 4441) )
trans_d <- subset(trans_d, !(NEWID %in% 4469) )

# NEWID 885: remove this plant
# NEWID 887: THIS IS A DORMANT in 2010!!!!
# NEWID 1129: Don't use it after 2012
# NEWID 2789: Don't use it! MORE LIKELY FOR THIS TO BE A NEW PLANT. 
# NEWID 4095: DORMANT IN 2012. This is a clear case
# NEWID 3534: DORMANT IN 2012. This is a clear case
subset(trans_d, STAGE_t0 == "D" & STAGE_t1 == "R")
# trans_d <- subset(trans_d, !(NEWID %in% 885) )
trans_d <- subset(trans_d, !(NEWID %in% 1129 & year > 2012) )
trans_d <- subset(trans_d, !(NEWID %in% 2789) )

# Actually dormant individuals
# trans_d <- mutate(trans_d,
#                   STAGE_t0 = replace(STAGE_t0, NEWID==887 & year == 2010, "DORM"),
#                   STAGE_t1 = replace(STAGE_t1, NEWID==887 & year == 2009, "DORM") )
trans_d <- mutate(trans_d,
                  STAGE_t0 = replace(STAGE_t0, NEWID==4095 & year == 2012, "DORM"),
                  STAGE_t1 = replace(STAGE_t1, NEWID==4095 & year == 2011, "DORM") )
trans_d <- mutate(trans_d,
                  STAGE_t0 = replace(STAGE_t0, NEWID==3534 & year == 2012, "DORM"),
                  STAGE_t1 = replace(STAGE_t1, NEWID==3534 & year == 2011, "DORM") )

# the DORM/?? individual (3225) is clearly dormant (because subsequent year it was DORM)
trans_d <- mutate(trans_d,
                  STAGE_t0 = replace(STAGE_t0, STAGE_t0 == "DORM/??", "DORM"),
                  STAGE_t1 = replace(STAGE_t1, STAGE_t1 == "DORM/??", "DORM") )

# remove D_t0 to Counting_t0
# trans_d <- subset(trans_d, !(STAGE_t0 == "D" & STAGE_t1 == "Counting") )


# tests: catch suspicious transitions --------------------------------------------------------

# sensible transitions: copied and pasted from transition from 2008-2018.
sensible_trans <- data.frame( STAGE_t0 = c('DORM', 'DORM', 'DORM', 'NR', 'NR', 'NR', 'NR', 
                                           'NR', 'NR', 'R', 'R', 'R', 'R', 'R', 'R', 'SL', 
                                           'SL', 'SL', 'SL', 'SL', NA, NA, NA),
                              STAGE_t1 = c('DORM', 'R', NA, 'A', 'D', 'DORM', 'NR', 'R',
                                           NA, 'A', 'D', 'DORM', 'NR', 'R', NA, 
                                           'D', 'DORM', 'NR', 'R', NA, 'NR', 'R', 'SL'),
                              stringsAsFactors = F
                             ) %>% arrange( STAGE_t0, STAGE_t1)

# transitions that you have observed
obs_trans      <- trans_d %>% 
                    dplyr::select(STAGE_t0, STAGE_t1) %>% 
                    unique %>% 
                    arrange(STAGE_t0, STAGE_t1) 


# check whether transitions from STAGE_t0 to STAGE_t1 make sense, logically
check_transitions <- function(obs_df, sensib_df){

  obs_minus_sensib <- setdiff(obs_df, sensib_df)
  
  if( nrow(obs_minus_sensib) > 0 ){
    cat('This year you have transitions that MIGHT not make sense.\nReview the returned dataframe to make sure all suspicious transitions make sense.\nI suggest storing the output dataframe into an R object.\n\n')
    return( obs_minus_sensib )
  }
  
  if( nrow(obs_minus_sensib) == 0 ){
    print('All transitions from t0 to t1 make sense: you can go ahead and continue formatting')
  }
  
}

# test
check_transitions(obs_trans, sensible_trans)



# SIZE -------------------------------------------------------------------------------
trans_d <- trans_d %>%
              mutate( LEN_t0  =  as.numeric(LEN_t0),
                      LEN_t1  =  as.numeric(LEN_t1),
                      WID_t0  =  as.numeric(WID_t0),
                      WID_t1  =  as.numeric(WID_t1)  ) %>%
              mutate( area_t0 = (LEN_t0/2) * (WID_t0/2) * pi,
                      area_t1 = (LEN_t1/2) * (WID_t1/2) * pi ) %>%
              mutate( area_t0 = replace(area_t0, 
                                        STAGE_t0 %in% c('DORM'), NA),
                      area_t1 = replace(area_t1, 
                                        STAGE_t1 %in% c('DORM','D','A'), NA),
                      area_t0 = replace(area_t0, is.na(STAGE_t0),  NA),
                      area_t1 = replace(area_t1, is.na(STAGE_t1),  NA) )


# SURVIVAL ---------------------------------------------------------------------------

# Survival is NA UNLESS individuals have "D", "NR", "R", "DORM", or "A" status
trans_d <- trans_d %>%
              mutate( surv_t1 = NA ) %>%
              mutate( surv_t1 = replace(surv_t1, STAGE_t1 == "D", 0) ) %>%
              mutate( surv_t1 = replace(surv_t1, STAGE_t1 %in% c("NR","R","DORM","A"), 1) ) %>%
              # NAs when no reliable information at STAGE_t0
              mutate( surv_t1 = replace(surv_t1, is.na(STAGE_t0), NA) )


# FLOWERING --------------------------------------------------------------------------

# flowering at time t0 and t1 is NA UNLESS individuals have "R" and "NR" information 
trans_d <- trans_d %>%
              mutate( flow_t0 = NA, 
                      flow_t1 = NA ) %>%
              mutate( flow_t0 = replace(flow_t0, STAGE_t0 == "R",  1),
                      flow_t1 = replace(flow_t1, STAGE_t1 == "R",  1) ) %>%
              mutate( flow_t0 = replace(flow_t0, STAGE_t0 == "NR", 0),
                      flow_t1 = replace(flow_t1, STAGE_t1 == "NR", 0) ) %>% 
              # introduce NAs for 2015 fertility-only extra individual 
              mutate( flow_t0 = replace(flow_t0, grepl('FERT_15_',NEWID), NA),
                      flow_t1 = replace(flow_t1, grepl('FERT_15_',NEWID), NA) )
 

# FERTILITY ------------------------------------------------------------------------
numall_t0_v <- apply(dplyr::select(trans_d, NUMAB_t0, NUMCL_t0,NUMINT_t0),1, sum, na.rm=T)
numall_t1_v <- apply(dplyr::select(trans_d, NUMAB_t1, NUMCL_t1,NUMINT_t1),1, sum, na.rm=T)

# introduce numall (correctly calculated)
trans_d     <- trans_d %>%
                  mutate( NUMALL_t0 = numall_t0_v,
                          NUMALL_t1 = numall_t1_v ) 

# instances in which numall is unreliable (mostly due to NAs)
numall_misleading <- trans_d %>% 
                        subset( !(NUMALL_t0 == NUMRAC_t0) ) %>% 
                        subset( !(is.na(NUMAB_t0) & is.na(NUMCL_t0) & is.na(NUMINT_t0)) ) %>% 
                        dplyr::select(NEWID, year, NUMRAC_t0,NUMAB_t0,NUMCL_t0,NUMINT_t0) %>% 
                        unique


trans_d %>% subset(NUMAB_t0 == 99)

# tests on vial rate data ----------------------------------------------------------------

# Do all STAGE_t0 == NA correspond to area_t0 == NA?  
expect_true( subset(trans_d, is.na(STAGE_t0) ) %>% 
               .$area_t0 %>% 
               is.na %>% 
               all )

# Do all STAGE_t0 == NA correspond to surv_t1 == NA?  
expect_true( trans_d %>% 
               subset( is.na(STAGE_t0) ) %>%
               .$surv_t1 %>% 
               is.na %>% all )


# Visual tests --------------------------------------------------------------------------------

# Stage transitions with non-NA survival data
trans_d %>% 
  subset( !is.na(surv_t1) ) %>% 
  dplyr::select(STAGE_t0,STAGE_t1) %>% 
  unique %>% 
  arrange(STAGE_t0)  

# Seedlings for which there is no size information at time t0 
trans_d %>% 
  subset( STAGE_t0 == "SL" & is.na(area_t0)) %>%
  dplyr::select(LOCATION, year) %>% unique 

# seedlings for which there IS size information
trans_d %>% 
  subset( STAGE_t0 == "SL" & !is.na(area_t0)) %>% 
  nrow

#write final demography file! ------------------------------------------------------

trans_out <- trans_d %>% 
                # lower-case-it!
                setNames( tolower(names(.)) ) %>% 
                # remove outlier (probably not LUTI but different species)
                subset( newid != 905 ) %>% 
                # label year as "transition year"
                mutate( transition = paste(year, year+1, sep = '-') )

# 2008_2018 data
write.csv(trans_out, "data/lupine_08_18.csv", row.names=F)


# all data -------------------------------------------------------------------------

lupine_08 <- read.csv('data/lupine_05_08.csv') %>% 
                rename( numrac_t0 = notab_t0,
                        numrac_t1 = notab_t1 )
lupine_18 <- read.csv('data/lupine_08_18.csv')

# remove doubles: -----------------------------------------------------------------

# duplicates found by Valentin 
double_id <- read.csv( 'results/explorative_graphs/duplicated_ids.csv' )

# Same doubles in both 2007 and 2008
double_df <- intersect( select(lupine_08, year, newid, location),
                        select(lupine_18, year, newid, location) ) 
expect_equal( count(double_df, newid) %>% .$n %>% unique, 2)

# the differences are the individuals I remove below
setdiff( double_df$newid, 
         double_id$x )


# Doubles when year == 2008 should be removed from lupine_08

# select ids in 2008
select_id_08 <- function(x_df, ii){
  x_df %>% subset( newid == double_id$x[ii] & year == 2008 ) %>% 
    select(newid, year, location, stage_t0,area_t0, numrac_t0, flow_t0,
                                  stage_t1,area_t1, numrac_t1, flow_t1, surv_t1)
}

# test that all _t0 values correspond (they are the same)
compare_t0_08 <- function(ii){
  all( select_id_08(lupine_08, ii) %>% select(area_t0, numrac_t0, flow_t0) == 
       select_id_08(lupine_18, ii) %>% select(area_t0, numrac_t0, flow_t0) )
}
expect_true( sapply(1:length(double_id$x), compare_t0_08) %>% all )

# all stage_t1 info from _08 should be NA
compare_t1_08 <- function(ii){
  df_08 <- select_id_08(lupine_08, ii) %>% select(stage_t1, area_t1, numrac_t1, flow_t1, surv_t1)  
  is.na(df_08) %>% sum
}
expect_true( (sapply(1:length(double_id$x), compare_t1_08) == 5) %>% all )


# Doubles when year == 2007 should be removed from lupine_18

# select ids in 2007
select_id_07 <- function(x_df, ii){
  x_df %>% subset( newid == double_id$x[ii] & year == 2007 ) %>% 
    select(newid, year, location, stage_t0,area_t0, numrac_t0, flow_t0,
                                  stage_t1,area_t1, numrac_t1, flow_t1, surv_t1)
}
  

# test that all _t1 (but surv_t1) values correspond (they are the same)
compare_t1_07 <- function(ii){
  all( select_id_07(lupine_08, ii) %>% select(area_t1, numrac_t1, flow_t1) == 
       select_id_07(lupine_18, ii) %>% select(area_t1, numrac_t1, flow_t1) )
}
expect_true( sapply(1:length(double_id$x), compare_t1_07) %>% all )

# all "_t0" info AND surv_t1 from _18 should be NA (surv_t1 needs 2007 data...)
compare_t0_07 <- function(ii){
  df_18 <- select_id_07(lupine_18, ii) %>% select(stage_t0, area_t0, numrac_t0, flow_t0, surv_t1)  
  is.na(df_18) %>% sum
}
expect_true( (sapply(1:length(double_id$x), compare_t0_07) == 5) %>% all )


# put it all together ------------------------------------------------------------

# Remove doubles when year == 2008 should be removed from lupine_08
lupine_08 <- subset(lupine_08, !(newid %in% double_id$x & year == 2008) )

# Remove doubles when year == 2007 should be removed from lupine_18
lupine_18 <- subset(lupine_18, !(newid %in% double_id$x & year == 2007) )

# put i
trans_all <- bind_rows(lupine_08, lupine_18) %>% 
                # remove dodgy data points (adults in 2005, transform into seedling post 2008)
                subset( newid != 777 ) %>% 
                subset( newid != 746 ) %>% 
                subset( newid != 735 ) 


# Calculate AGE information 

# get birth days of all individuals which were observed as seedlings
births_df <- trans_all %>%
                mutate( birth_year = year ) %>%
                subset( stage_t0 == 'SL' ) %>% 
                dplyr::select( newid, birth_year ) %>% 
                unique

# introduce age information, and calculate log_area
trans_out <- trans_all %>% 
                left_join( births_df ) %>% 
                mutate( age_t0 = year - birth_year ) %>% 
                mutate( age_t0 = replace(age_t0, is.na(stage_t0), NA) ) %>% 
                mutate( area_t0 = replace(area_t0, area_t0 == 0, NA),
                        area_t1 = replace(area_t1, area_t1 == 0, NA) ) %>% 
                mutate( log_area_t0 = log(area_t0),
                        log_area_t1 = log(area_t1) ) %>% 
                # remove useless columns
                select(-UID, 
                       -len_t0, -wid_t0, size_t0, -numall_t0, 
                       -len_t1, -wid_t1, size_t1, -numall_t1 )

write.csv(trans_out, 'data/lupine_all.csv', row.names=F)


# plots ----------------------------------------------------------------------------

# files
surv <- subset(trans_d, !is.na(surv_t1) )
flow <- subset(trans_d, !is.na(flow_t1) )
fert <- subset(trans_d, !is.na(flow_t1) & flow_t1 == 1 )
grow <- subset(trans_d, !(STAGE_t0 %in% c("DORM", "NF")) & 
                        !(STAGE_t1 %in% c("D", "NF", "DORM")) ) %>%
        # remove zeroes from area_t0 and area_t1
        subset( area_t0 != 0) %>%
        subset( area_t1 != 0)


# raw scale
par(mfrow=c(2,2), mar = c(3,3, 0.5, 0.5), mgp = c(1.8, 0.7, 0) )
plot(jitter(surv_t1) ~ area_t0, pch=1, ylab = "Survival", data = surv )
plot(area_t1 ~ area_t0, pch=1, data = grow)
plot(jitter(flow_t1) ~ area_t0, pch=1, ylab = "Flowering", data = flow )
plot(NUMRAC_t1 ~ area_t0, pch=1, data = fert)

# log transformed
par(mfrow=c(2,2), mar = c(3,3, 0.5, 0.5), mgp = c(1.8, 0.7, 0) )
plot(jitter(surv_t1) ~ log(area_t0), pch=1, ylab = "Survival", 
     data = subset(surv, STAGE_t0 != 'SL') )
plot(log(area_t1) ~ log(area_t0), pch=1, data = grow)
plot(jitter(flow_t1) ~ log(area_t1), pch=1, ylab = "Flowering", 
     data = subset(flow, STAGE_t0 != 'SL') )
plot(NUMRAC_t0 ~ log(area_t1), pch=1, data = fert)


fert <- fert %>% 
            mutate(log_area_t0 = log(area_t0),
                   log_area_t1 = log(area_t1) ) %>% 
            subset( log_area_t0 != -Inf)
          

fert1 %>% 
  subset( NUMALL_t1 == 0) %>% 
  select( NUMALL_t1, NUMAB_t1, NUMCL_t1, NUMINT_t1)

yearz <- fert$year %>% unique %>% .[-7]
sitez <- fert$sub_loc %>% unique
par(mfrow=c(8,10))
for(ss in 1:length(sitez)){
  for(ii in 1:length(yearz)){
    plot(NUMRAC_t1 ~ log_area_t1, pch=1, ylab = "Survival", 
         data = subset(fert, year == yearz[ii] & sub_loc == sitez[ss]) )
  }
}

par(mfrow=c(1,1))
plot(NUMRAC_t1 ~ log_area_t1, pch=1, ylab = "fert", 
         data = fert )

par(mfrow=c(4,2))
for(ii in 1:length(sitez)){
  plot(NUMRAC_t1 ~ log_area_t1, pch=1, ylab = "Survival", 
       data = subset(fert, sub_loc == sitez[ii]) )
}

plot(NUMRAC_t1 ~ log_area_t1, pch=1, ylab = "Survival", 
       data = subset(fert, sub_loc == sitez[ii]) )


par(mfrow=c(1,1))
plot(NUMALL_t1 ~ log_area_t0, pch=1, ylab = "fert", data = fert1 )

par(mfrow=c(1,1))
plot(NUMALL_t1 ~ log_area_t1, pch=1, ylab = "fert", data = fert1 )



par(mfrow=c(8,10))
for(ss in 1:length(sitez)){
  for(ii in 1:length(yearz)){
    plot(jitter(surv_t1) ~ log(area_t0), pch=1, ylab = "Survival", 
         data = subset(surv, year == yearz[ii] & sub_loc == sitez[ss],
                             STAGE_t0 != 'SL') )
  }
}


subset(surv, year == yearz[9] & sub_loc == sitez[2]) %>% 
  as.data.frame %>% 
  subset( STAGE_t0 == 'SL') 
  
