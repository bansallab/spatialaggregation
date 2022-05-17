## Name: Elizabeth Lee
## Date: 2/13/17
## Function: Functions for cleaning disease burden response data at the state level for INLA
## Filenames: 
## Data Source: 
## Notes: 
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### functions for model data cleaning ################################
require(dplyr); require(tidyr); require(readr); require(DBI); require(RMySQL)
require(igraph)

##### RESPONSE VARIABLES ##########################################

################################  

cleanR_iliEarly_shift1_st <- function(filepathList){
  # clean response variable: ilinDt.early
  print(match.call())
  
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean data
  iliEarly_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.early", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  print(filepathList$path_response_st)
  print(summary(iliEarly_data))
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there
  
  return_data <- full_join(iliEarly_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>% 
    mutate(y1 = y+1) %>%
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9)
  
  return(return_data)
}
################################

cleanR_iliPeak_shift1_st <- function(filepathList){
  # clean response variable: ilinDt.peak plus 1 (so it is comparable with iliSum); 9/15/17
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_st(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean burden data
  iliPeak_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.peak", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  # merge final data
  return_data <- full_join(iliPeak_data, pop_data, by = c("season", "abbr_st")) %>%
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>%
    mutate(y1 = y+1)  %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
################################  

cleanR_iliEarly_irDt_shift1_st <- function(filepathList){
  # clean response variable: irDt.early
  print(match.call())
  
  # grab disease burden metric (e.g., irDt): match "irDt" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean data
  iliEarly_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.early", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  print(filepathList$path_response_st)
  print(summary(iliEarly_data))
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there
  
  return_data <- full_join(iliEarly_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>% 
    mutate(y1 = y+1) %>%
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9)
  
  return(return_data)
}
################################

cleanR_iliPeak_irDt_shift1_st <- function(filepathList){
  # clean response variable: irDt.peak plus 1 (so it is comparable with iliSum); 9/15/17
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_st(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., irDt): match "irDt" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean burden data
  iliPeak_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.peak", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  # merge final data
  return_data <- full_join(iliPeak_data, pop_data, by = c("season", "abbr_st")) %>%
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>%
    mutate(y1 = y+1)  %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
################################  

cleanR_wksToEpi_st <- function(filepathList){
  # clean response variable: wks.to.epi; 7/5/17
  print(match.call())
  
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there

  # clean data
  wksToEpi_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == "wks.to.epi") %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)

  # merge final data  
  return_data <- full_join(wksToEpi_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>% 
    mutate(y1 = ifelse(y>0, y, NA)) %>%
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9)
  
  return(return_data)
}
################################  

cleanR_wksToPeak_st <- function(filepathList){
  # clean response variable: wks.to.peak from beginning of flu period; 7/14/17
  print(match.call())
  
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there

  # clean data
  wksToPeak_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == "wks.to.peak") %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)

  # merge final data  
  return_data <- full_join(wksToPeak_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>% 
    mutate(y1 = ifelse(y>0, y, NA)) %>%
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9)
  
  return(return_data)
}
################################  

cleanR_iliSum_shift1_st <- function(filepathList){
  # clean response variable: ilinDt.sum
  print(match.call())
  
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean data
  iliSum_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  print(filepathList$path_response_st)
  print(summary(iliSum_data))
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there
  
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>% 
    mutate(y1 = y+1) %>%
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9)
  
  return(return_data)
}
################################  

cleanR_iliSum_shift1_st_aggBias <- function(filepathList){
  # (original data) clean response variable: ilinDt.sum
  print(match.call())
  
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean data
  iliSum_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  print(filepathList$path_response_st)
  print(summary(iliSum_data))
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there
  
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop 
    mutate(logy1_st = log(y+1)) %>%
    filter(season >= 3 & season <= 9) %>%
    select(abbr_st, season, logy1_st) %>%
    rename(st = abbr_st)
  
  return(return_data)
}
################################  

cleanR_iliSum_irDt_shift1_st <- function(filepathList){
  # clean response variable: irDt.sum
  print(match.call())
  
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)
  # clean data
  iliSum_data <- read_csv(filepathList$path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden, abbr_st = state)
  
  print(filepathList$path_response_st)
  print(summary(iliSum_data))
  pop_data <- clean_pop_st(filepathList) # 4/12/16 all 51 pops are there
  
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "abbr_st")) %>% # 4/12/16 full_join so pops don't drop
    select(fips_st, abbr_st, state, lat, lon, season, year, pop, y, has.epi) %>% 
    mutate(y1 = y+1) %>%
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9)
  
  return(return_data)
}
################################  

##### SAMPLING EFFORT DATA ##########################################

cleanX_priorBurden_st <- function(filepathList){
  ## clean data variable: ilinDt.sum for previous year; new 11/10/16
  print(match.call())
  
  ## grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  ## 5/28/18 scales modesl all utilize irDt measure
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_st, "_")[[1]], value=T)

  # grab total pop response, not adult or child one
  total_path_response_st <- gsub("_child", "", gsub("_adult", "", filepathList$path_response_st))
  
  # crosswalk to fips_st data
  cw <- cw_zip3_st()

  # clean burden data into prior immunity
  output <- read_csv(total_path_response_st, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(priorBurden = burden, abbr_st = state) %>%
    mutate(season = season + 1) %>%
    left_join(cw, by = "abbr_st") %>%
    distinct(fips_st, season, priorBurden) %>%
    filter(season >= 3 & season <= 9)
  
  return(output)
}
##########################################

##### COVARIATE DATA RELYING ON RESPONSE VARIABLE ##########################################
cleanX_protectedFromPrevSeason_st <- function(filepathList){
  # clean variable indicating protection conferred from infection during previous flu season: incorporates previous and current season subtype/type distributions (H1/H3/B) and strain similarity
  print(match.call())
  
  priorburdenDat <- cleanX_priorBurden_st(filepathList) # fips_st, season, priorBurden
  protectedPropDat <- cleanX_multsrcSubtypeDistrStrainSim_reg() # season, region, estImmuneProp
  cw <- cleanX_cdcFluview_H3_region() %>% select(region, fips) %>% rename(fips_st = fips) %>% distinct # region, fips_st
  
  output <- left_join(priorburdenDat, cw, by = c("fips_st")) %>%
    left_join(protectedPropDat, by = c("season", "region")) %>%
    rowwise %>% 
    mutate(protectionPrevSeason = prod(priorBurden, estImmuneProp)) %>%
    mutate(protectionPrevSeason = ifelse(is.na(protectionPrevSeason), 0, protectionPrevSeason)) %>%
    select(fips_st, season, protectionPrevSeason) %>%
    ungroup
  
  return(output)
}

##### REFERENCE DATA ##########################################
cw_zip3_st <- function(){
  # spatial crosswalk for zip3-state
  print(match.call())

  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "spatialcw_zip3_state")
  # sel.statement <- "Select * from spatialcw_zip3_state limit 5"
  sel.statement <- "SELECT fips_st, zip3, abbr_st FROM spatialcw_zip3_state"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  output <- tbl_df(dummy) 
  return(output)
}

################################
clean_pop_st <- function(filepathList){
  # clean pop data with state abbr
  print(match.call())
  
  pop_data <- clean_pop_st_plain()
  abbr_data <- read_csv(filepathList$path_abbr_st, col_types = "ccc")
  coord_data <- read_csv(filepathList$path_latlon_st, col_types = "cdd", col_names = c("abbr_st", "lat", "lon"), skip = 1)
  
  dummy <- left_join(pop_data, abbr_data, by = c("fips_st" = "FIPS")) %>%
    rename(abbr_st = Abbreviation, state = State) %>%
    mutate(season = as.numeric(substring(year, 3, 4))) %>%
    select(fips_st, abbr_st, state, season, year, pop) %>%
    filter(season > 1)
  
  fulldata <- left_join(dummy, coord_data, by = "abbr_st")

  return(fulldata)
}

################################
clean_pop_st_plain <- function(){
  # clean pop data at state level
  print(match.call())
  
  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_state")
  # sel.statement <- "Select * from demog_Census_agePop_state limit 5"
  sel.statement <- "SELECT fips as fips_st, agegroup, year, pop FROM demog_Census_agePop_state WHERE scale = 'state' and agegroup = 'total'"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    select(fips_st, year, pop)
  return(output)
}

################################
clean_pop_adult_st_plain <- function(){
  # clean adult pop (20-69 yo) data at state level -- fips, all years, pop (for write_loess_fits_ILIn_age_cty.R)
  print(match.call())
  
  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_state")
  # sel.statement <- "Select * from demog_Census_agePop_state limit 5"
  sel.statement <- "SELECT fips as fips_st, agegroup, year, pop FROM demog_Census_agePop_state WHERE scale = 'state' and agegroup = 'adult'"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    select(fips_st, year, pop)
  return(output)
}

################################
clean_pop_child_st_plain <- function(){
  # clean child pop (5-19 yo) data at state level -- fips, all years, pop (for write_loess_fits_ILIn_age_cty.R)
  print(match.call())
  
  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_state")
  # sel.statement <- "Select * from demog_Census_agePop_state limit 5"
  sel.statement <- "SELECT fips as fips_st, agegroup, year, pop FROM demog_Census_agePop_state WHERE scale = 'state' and agegroup = 'child'"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    select(fips_st, year, pop)
  return(output)
}

################################

# #### testing area ################################
# dbCodeStr <- "_ilinDt_Octfit_span0.4_degree2"
# 
# setwd(dirname(sys.frame(1)$ofile))
# setwd('../reference_data')
# path_latlon_st <- paste0(getwd(), "/state_latlon.csv")
# path_abbr_st <- paste0(getwd(), "/state_abbreviations_FIPS.csv")
# setwd("../R_export")
# path_response_st <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_st.csv", dbCodeStr))
# 
# # put all paths in a list to pass them around in functions
# path_list <- list(path_abbr_st = path_abbr_st,
#                   path_latlon_st = path_latlon_st,
#                   path_response_st = path_response_st)
# setwd(dirname(sys.frame(1)$ofile))
# 
# db <- cleanR_iliSum_shift1_st(path_list)
# priordb <- cleanX_priorBurden_st(path_list)
# priorImm <- cleanX_protectedFromPrevSeason_st(path_list)
# cw <- cw_zip3_st()
# fullpop <- clean_pop_st(path_list)
# pop <- clean_pop_st_plain()
# aPop <- clean_pop_adult_st_plain()
# cPop <- clean_pop_child_st_plain()
