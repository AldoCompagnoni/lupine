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


# Global changes -------------------------------------------------------------------
subset(trans_d, is.na(STAGE_t0) & STAGE_t1 == "NR")$area_t0
subset(trans_d, STAGE_t0 == "NF" )$surv_t1
subset(trans_d, is.na(STAGE_t0) ) %>% 
         dplyr::select(STAGE_t0, STAGE_t1) %>%
         unique

# automatic tests ----------------------------------------------------------------------------

# Do all is.na(STAGE_t0) correspond to surv_t1 == NA?  
expect_true( trans_d %>% 
               subset( is.na(STAGE_t0) ) %>%
               .$surv_t1 %>% 
               is.na %>% all )
# Do all STAGE_t0 == "NF" correspond to surv_t1 == NA?  
expect_true( trans_d %>% 
               subset( STAGE_t0 == "NF" ) %>%
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

trans_d %>% 
  subset( STAGE_t0 == "SL" & !is.na(area_t0)) %>% 
  nrow


# seedlings with area_t0 == 0
no_area_sl <- trans_d %>% 
                subset( STAGE_t0 == "SL" & is.na(area_t0)) %>%
                dplyr::select(FID, NEWID, LOCATION, year)
area_sl <- trans_d %>% 
                subset( STAGE_t0 == "SL" & !is.na(area_t0)) %>%
                dplyr::select(FID, NEWID, LOCATION, year, area_t0)
area_0_sl <- subset(area_sl, area_t0 == 0)
write.csv(no_area_sl, "no_area_seedlings_luti.csv", row.names=F)
write.csv(area_sl,    "area_seedlings_luti.csv",    row.names=F)
write.csv(area_0_sl,  "area_0_seedlings_luti.csv",  row.names=F)

subset(indiv_l[[3]], FID == 6690 & NEWID == 5986)
       
# missing size dat for year 12 
y12_cols <- names(dat)[grep("Y12",names(dat)) ]


# write final demography file! ------------------------------------------------------

trans_d <- trans_d %>% 
              # lower-case-it!
              setNames( tolower(names(.)) ) %>% 
              # remove outlier (probably not LUTI but different species)
              subset( newid != 905 ) %>% 
              # label year as "transition year"
              mutate( transition = paste(year, year+1, sep = '-') )

write.csv(trans_d, "data/lupine_08_17.csv", row.names=F)


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

# plots

# raw scale
par(mfrow=c(2,2), mar = c(3,3, 0.5, 0.5), mgp = c(1.8, 0.7, 0) )
plot(jitter(surv_t1) ~ area_t0, pch=1, ylab = "Survival", data = surv )
plot(area_t1 ~ area_t0, pch=1, data = grow)
plot(jitter(flow_t1) ~ area_t0, pch=1, ylab = "Flowering", data = flow )
plot(NUMALL_t1 ~ area_t0, pch=1, data = fert)

# log transformed
par(mfrow=c(2,2), mar = c(3,3, 0.5, 0.5), mgp = c(1.8, 0.7, 0) )
plot(jitter(surv_t1) ~ sqrt(area_t0), pch=1, ylab = "Survival", data = surv )
plot(sqrt(area_t1) ~ sqrt(area_t0), pch=1, data = grow)
plot(jitter(flow_t1) ~ sqrt(area_t0), pch=1, ylab = "Flowering", data = flow )
plot(NUMALL_t1 ~ sqrt(area_t0), pch=1, data = fert)



# models -----------------------------------------------------------------------------------

gr_mod <- lm(sqrt(area_t1) ~ sqrt(area_t0), data = grow)
fl_mod <- glm(flow_t0 ~ log(area_t0), family = "binomial", data = subset(flow, area_t0 !=0) )

# unequal variance?
pred_y  <- fitted(gr_mod)
resid_y <- resid(gr_mod)
plot(pred_y, resid_y)

# 
par(mfrow=c(1,1), mar = c(3,3, 0.5, 0.5), mgp = c(1.8, 0.7, 0) )
plot(jitter(flow_t0) ~ log(area_t0), pch=1, ylab = "Survival", data = surv )
x_seq  <- seq( min(log(flow$area_t0),na.rm=T),
               max(log(flow$area_t0),na.rm=T), length.out=100 )
lines(x_seq, inv.logit(coef(fl_mod)[[1]] + coef(fl_mod)[[2]]*x_seq) )




# checks/issues ----------------------------------------------------------------------------

# LAST CHECKS ----------------------------------------------------------------------

trans_d %>%
  dplyr::select(STAGE_t0,STAGE_t1,surv_t1) %>%
  unique %>%
  arrange(STAGE_t0, STAGE_t1) %>%
  as.data.frame

# # REMOVE THESE: few and we are not 100% sure these died that year.
# subset(trans_d, is.na(STAGE_t0) & STAGE_t1 == "D")
# subset(trans_d, FID == 198 & NEWID == 3719) %>% dplyr::select(year, STAGE_t0, STAGE_t1,area_t0,area_t1)
# subset(indiv_l[[1]], FID == 198 & NEWID == 3719)
# 
# # REMOVE THESE: all instances of gaps in data ("999")
# subset(trans_d, is.na(STAGE_t0) & STAGE_t1 == "R")
# 
# subset(trans_d, FID == 281 & NEWID == 1048) %>% dplyr::select(year, STAGE_t0, STAGE_t1,area_t0,area_t1)
# subset(indiv_l[[1]], FID == 281 & NEWID == 1048)
# subset(trans_d, FID == 11 & NEWID == 2584) %>% dplyr::select(year, STAGE_t0, STAGE_t1,area_t0,area_t1)
# subset(indiv_l[[1]], FID == 11 & NEWID == 2584)
# subset(trans_d, FID == 3393 & NEWID == 853) %>% dplyr::select(year, STAGE_t0, STAGE_t1)
# subset(indiv_l[[1]], FID == 3393 & NEWID == 853)
# subset(trans_d, FID == 6499 & NEWID == 3538) %>% dplyr::select(year, STAGE_t0, STAGE_t1)
# subset(indiv_l[[1]], FID == 6499 & NEWID == 3538)
# 
# subset(trans_d, is.na(STAGE_t0) & STAGE_t1 == "NR")
# subset(trans_d, FID ==  739 & NEWID == 4660) %>% dplyr::select(year, STAGE_t0, STAGE_t1,area_t0,area_t1)
# subset(indiv_l[[1]], FID == 739 & NEWID == 4660)
# subset(trans_d, FID == 825 & NEWID == 5006) %>% dplyr::select(year, STAGE_t0, STAGE_t1,area_t0,area_t1)
# subset(indiv_l[[1]], FID == 825 & NEWID == 5006)
# subset(trans_d, FID == 2041 & NEWID == 8741) %>% dplyr::select(year, STAGE_t0, STAGE_t1,area_t0,area_t1,flow_t0)
# subset(indiv_l[[1]], FID == 2041 & NEWID == 8741)
# 
# subset(trans_d, STAGE_t0 == "SL" & STAGE_t1 == "Elim") %>% dplyr::select(area_t0, area_t1)

# flowering status
trans_d %>% 
  dplyr::select(STAGE_t0, flow_t0) %>%
  unique
trans_d %>% 
  dplyr::select(STAGE_t1, flow_t1) %>%
  unique

# area and status
trans_d %>% 
  dplyr::select(STAGE_t0, area_t0) %>%
  subset( is.na(area_t0) ) %>%
  unique
trans_d %>% 
  dplyr::select(STAGE_t1, area_t1) %>%
  subset( is.na(area_t1) ) %>%
  unique

# total individuals censused
trans_d %>% 
  dplyr::select(LOCATION, id) %>%
  unique %>%
  group_by( LOCATION ) %>%
  summarize( reps = n() ) %>% 
  as.data.frame


# subset(trans_d, STAGE_t1 == "D") %>% 
#   dplyr::select(STAGE_t0, STAGE_t1) %>%
#   unique
# subset(trans_d, STAGE_t1 == "R") %>% 
#   dplyr::select(STAGE_t0, STAGE_t1) %>%
#   unique
# subset(trans_d, STAGE_t1 == "NR") %>% 
#   dplyr::select(STAGE_t0, STAGE_t1) %>%
#   unique
# subset(trans_d, STAGE_t0 == "SL") %>% 
#   dplyr::select(STAGE_t0, STAGE_t1) %>%
#   unique



# very large individual - is this "correct"?
subset(trans_d, area_t0 > 50000)
subset(trans_d, area_t1 > 50000)


trans_d %>% 
  select(STAGE_t0, STAGE_t1, flow_t1) %>% 
  arrange(STAGE_t0,STAGE_t1) %>% 
  unique

# other issues
trans_d %>% 
  dplyr::select(STAGE_t0, STAGE_t1) %>% 
  unique %>%
  subset( STAGE_t0 != "Elim" & STAGE_t0 != "NF" & STAGE_t0 != "TAGLATER" &
          STAGE_t0 != "D" & STAGE_t0 != "R" & STAGE_t0 != "R/Elim" &
          STAGE_t0 != "Counting" & STAGE_t0 != "DORM" ) %>%
  arrange( STAGE_t0, STAGE_t1)




## 
## CHECK OTHER ISSUES

# STAGE

# LEN: what is 999? (or viceversa, 9999?)

# 


# Metadata on stage codes
stage_codes <- d_need %>%
                  select( names(dat)[ grep("Y[0-9][0-9]STAGE", names(dat)) ] ) %>%
                  as.data.frame %>%
                  unlist %>%
                  unique

# 'NUMRAC' variable
# a. ONLY IN 2008: 999 are NAs; 111 reproductive but did not count.
# b. FROM 2009 ON: 9999 are NAs; 1111 reproductive but did not count. 

# Needed metadata
NUMAB   
NUMCL
NUMINT
SQUOOP
SQPCOMM
HAB
SUCCSTG
CHCK
CHK 

# stages codes:
SL
R
NR
TAGLATER
999
D
NF/E
NF
DORM     
BLNOVR
R/Elim
999/Elim
Elim
NF/BLNOVR
Counting
A
LUCH
