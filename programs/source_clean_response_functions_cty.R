## Name: Elizabeth Lee
## Date: 5/27/16
## Function: Functions for cleaning disease burden response data at the county level for INLA
## Filenames: 
## Data Source: 
## Notes: 7/28/16 downscaleDB suffix refers to downscaling procedure for the already-processed disease burden metrics; renamed functions from cleanR_iliSum_cty and cleanR_iliPeak_cty
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### functions for model data cleaning ################################
require(dplyr); require(tidyr); require(readr); require(DBI); require(RMySQL)
require(igraph)

##### COUNTY-LEVEL VARIABLES ##########################################

cleanR_iliEarly_shift1_cty <- function(filepathList){
  # clean response variable: ilinDt.early plus 1 (so it is comparable with iliSum); 9/14/17
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliEarly_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.early", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliEarly_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1) %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_iliPeak_shift1_cty <- function(filepathList){
  # clean response variable: ilinDt.peak plus 1 (so it is comparable with iliSum); 9/15/17
  print(match.call())
  
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliPeak_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.peak", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliPeak_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1)  %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_iliEarly_irDt_shift1_cty <- function(filepathList){
  # clean response variable: irDt.early plus 1 (so it is comparable with iliSum); 9/14/17
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliEarly_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.early", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliEarly_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1) %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_iliPeak_irDt_shift1_cty <- function(filepathList){
  # clean response variable: irDt.peak plus 1 (so it is comparable with iliSum); 10/15/17
  print(match.call())
  
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliPeak_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.peak", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliPeak_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1)  %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_wksToEpi_cty <- function(filepathList){
  # clean response variable: wks.to.epi; 3/31/17
  print(match.call())

  # pop data: fips, county, st, season, year, pop, lat, lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # clean burden data
  wksToEpi_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == "wks.to.epi") %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(wksToEpi_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = ifelse(y>0, y, NA)) %>% # 10/3/16
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16

  return(return_data)
}
##########################################

cleanR_wksToPeak_cty <- function(filepathList){
  # clean response variable: wks.to.peak (from epi onset); 7/14/17
  print(match.call())

  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # clean burden data
  wksToPeak_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == "wks.to.peak") %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(wksToPeak_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = ifelse(y>0, y, NA)) %>% # 10/3/16
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16

  return(return_data)
}
##########################################

cleanR_iliSum_shift1_cty <- function(filepathList){
  # clean response variable: ilinDt.sum plus 1; 12/15/16
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliSum_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1) %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_iliSum_irDt_shift1_cty <- function(filepathList){
  # clean response variable: irDt.sum plus 1; 5/27/18
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., irDt): match "irDt" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliSum_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1) %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_iliSum_shift1_cty_aggBias <- function(filepathList){
  # (original data) clean response variable: ilinDt.sum plus 1; 12/15/16 in preparation for aggBias modeling
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("ili+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)
  # clean burden data
  iliSum_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1) %>% # 12/15/16 add 1 to all seasonal intensity values
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

cleanR_iliRate_shift1_cty <- function(filepathList){
  # clean response variable: iliRate + 1; 7/6/17
  print(match.call())
  
  # pop data: fips, county, st, season, year, pop, lat lon
  pop_data <- clean_pop_cty(filepathList)
  
  # 7/18/16: add incl.analysis indicator
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- "iliRate"
  # clean burden data
  iliSum_data <- read_csv(filepathList$path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s", dbCode)) %>%
    select(-metric) %>%
    rename(y = burden)
  
  # merge final data
  return_data <- full_join(iliSum_data, pop_data, by = c("season", "fips")) %>%
    select(fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis) %>%
    mutate(y1 = y+1) %>% # 12/15/16 add 1 to all seasonal intensity values
    group_by(season) %>%
    mutate(E = weighted.mean(y1, pop, na.rm = TRUE)) %>%
    ungroup %>%
    filter(season >= 3 & season <= 9) # 12/12/16
  
  return(return_data)
}
##########################################

##### SAMPLING EFFORT DATA ##########################################
cleanO_imsCoverage_cty <- function(){
  # clean IMS Health adjusted physician coverage (database coverage) and visits per physician or visits per population (care-seeking behavior) from zip3 to county level, using overlapping pop bw zip3 & county as a weight for the weighted average
  # 1/5/17 rm careseek variables from this function, see source_clean_data_functions.R/cleanO_imsCareseekTot-Adult-Child_cty
  print(match.call())

  # spatial crosswalk: fips, zip3, proportion (of overlap in zip3 & fips population)
  cw <- cw_zip3_cty() 
  popDat <- clean_pop_cty_plain() # fips, year, pop
  
  # import physician coverage data
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "IMS_physicianCoverage_zip3")
  # sel.statement <- "Select * from IMS_physicianCoverage_zip3 limit 5"
  sel.statement <- "SELECT year, zip3, adjProviderCoverage, sampViz, sampProv FROM IMS_physicianCoverage_zip3"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)

  # clean zip3 coverage to county coverage
  covDat <- dummy %>%
    full_join(cw, by = "zip3") %>%
    group_by(fips, year) %>%
    summarise(zOverlaps = length(zip3), adjProviderCoverage = weighted.mean(adjProviderCoverage, proportion, na.rm = TRUE), sampViz = weighted.mean(sampViz, proportion, na.rm = TRUE), sampProv = weighted.mean(sampProv, proportion, na.rm = TRUE)) %>% 
    ungroup %>%
    filter(!is.na(fips)) %>% 
    mutate(adjProviderCoverage = ifelse(is.na(adjProviderCoverage), 0, adjProviderCoverage)) %>% # 9/27/16 for glm: 0 if NA
    mutate(visitsPerProvider = ifelse(is.na(sampViz), 0, sampViz/sampProv)) %>% # 9/27/16 for glm: 0 if NA
    left_join(popDat, by = c("fips", "year")) %>%
    mutate(visitsPerPop = ifelse(is.na(sampViz), 0, sampViz/pop)) %>% # 9/27/16 for glm: 0 if NA
    select(fips, year, adjProviderCoverage) %>% # 1/7/16 rm select for careseek terms
    arrange(fips, year)
  return(covDat)
}
##########################################

cleanX_priorBurden_cty <- function(filepathList){
  # clean data variable: ilinDt.sum for previous year; new 11/10/16
  print(match.call())
  
  # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
  dbCode <- grep("irDt+", strsplit(filepathList$path_response_cty, "_")[[1]], value=T)

  # grab total pop response, not adult or child one
  total_path_response_cty <- gsub("_child", "", gsub("_adult", "", filepathList$path_response_cty))
  
  # clean burden data into prior immunity
  output <- read_csv(total_path_response_cty, col_types = "icllcd") %>%
    filter(metric == sprintf("%s.sum", dbCode)) %>%
    select(-metric) %>%
    rename(priorBurden = burden) %>%
    mutate(season = season + 1) %>%
    select(fips, season, priorBurden) %>%
    filter(season >= 3 & season <= 9)
  
  return(output)
}
##########################################

##### COVARIATE DATA RELYING ON RESPONSE VARIABLE ##########################################
cleanX_protectedFromPrevSeason_cty <- function(filepathList){
  # clean variable indicating protection conferred from infection during previous flu season: incorporates previous and current season subtype/type distributions (H1/H3/B) and strain similarity
  print(match.call())
  
  priorburdenDat <- cleanX_priorBurden_cty(filepathList) %>% mutate(fips_st = substr(fips, 1, 2)) # fips, season, priorBurden
  protectedPropDat <- cleanX_multsrcSubtypeDistrStrainSim_reg() # season, region, estImmuneProp
  cw <- cleanX_cdcFluview_H3_region() %>% select(region, fips) %>% rename(fips_st = fips) %>% distinct # region, fips_st
  
  output <- left_join(priorburdenDat, cw, by = c("fips_st")) %>%
    left_join(protectedPropDat, by = c("season", "region")) %>%
    rowwise %>% 
    mutate(protectionPrevSeason = prod(priorBurden, estImmuneProp)) %>%
    mutate(protectionPrevSeason = ifelse(is.na(protectionPrevSeason), 0, protectionPrevSeason)) %>%
    select(fips, season, protectionPrevSeason) %>%
    ungroup
  
  return(output)
}

##### REFERENCE DATA ##########################################
clean_graphIDx <- function(filepathList, spatial_scale){
  # 10/30/16 import spatial crosswalk for fips-graph IDs
  print(match.call())
  
  if (spatial_scale == "county"){
    graphIdxDat <- read_csv(filepathList$path_graphIdx_cty) %>%
      select(fips, graphIdx)
  } else if (spatial_scale == "state"){
    graphIdxDat <- read_csv(filepathList$path_graphIdx_st) %>%
      select(fips_st, graphIdx_st)
  }
   
  return(graphIdxDat)
}

################################
clean_ctyCommmuter_stPassenger_graph <- function(filepathList, spatial_scale){
  # 12/19/16 import spatial dependence for commuting flows at county level or flight passenger flows at state level
  print(match.call())
  
  if (spatial_scale == "county"){
    g <- read_graph(filepathList$path_graphExport_cty, format = "edgelist", directed = FALSE)
  } else if (spatial_scale == "state"){
    g <- read_graph(filepathList$path_graphExport_st, format = "edgelist", directed = FALSE)
  }
  
  adjMx <- as_adjacency_matrix(g, type = "both")
  
  return(adjMx)
}

################################

cw_zip3_cty <- function(){
  # spatial crosswalk for zip3-county
  print(match.call())
  
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "spatialcw_zip3_county")
  # sel.statement <- "Select * from spatialcw_zip3_county limit 5"
  sel.statement <- "SELECT fips, zip3, zinf_prop as proportion FROM spatialcw_zip3_county"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  output <- tbl_df(dummy) 
  return(output)
}

################################
clean_pop_cty <- function(filepathList){
  # clean pop data at county level
  print(match.call())
  
  # read coord data by county: reference_data/cty_pop_latlon.csv
  coord_data <- read_csv(filepathList$path_latlon_cty , col_types = "cc__dd", col_names = c("st", "fips", "lat", "lon"), skip = 1) 
  # read state name data: reference_data/state_abbreviations_FIPS.csv
  abbr_data <- read_csv(filepathList$path_abbr_st, col_types = "cci", col_names = c("state", "st", "stateID"), skip = 1) 

  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_county")
  # sel.statement <- "Select * from demog_Census_agePop_county limit 5"
  sel.statement <- "SELECT fips, county, agegroup, year, pop FROM demog_Census_agePop_county WHERE scale = 'county' and agegroup = 'total' and year >= 2002 and year <= 2009"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    mutate(season = as.numeric(substring(year, 3, 4))) %>%
    left_join(coord_data, by = "fips") %>%
    left_join(abbr_data, by = "st") %>% 
    select(fips, county, st, state, stateID, season, year, pop, lat, lon)
  return(output)
}

################################
clean_pop_cty_plain <- function(){
  # clean pop data at county level -- fips, all years, pop (for write_loess_fits_ILIn.R)
  print(match.call())
  
  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_county")
  # sel.statement <- "Select * from demog_Census_agePop_county limit 5"
  sel.statement <- "SELECT fips, county, agegroup, year, pop FROM demog_Census_agePop_county WHERE scale = 'county' and agegroup = 'total'"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    select(fips, year, pop)
  return(output)
}

################################
clean_pop_adult_cty_plain <- function(){
  # clean adult pop (20-69 yo) data at county level -- fips, all years, pop (for write_loess_fits_ILIn_age_cty.R)
  print(match.call())
  
  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_county")
  # sel.statement <- "Select * from demog_Census_agePop_county limit 5"
  sel.statement <- "SELECT fips, county, agegroup, year, pop FROM demog_Census_agePop_county WHERE scale = 'county' and agegroup = 'adult'"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    select(fips, year, pop)
  return(output)
}

################################
clean_pop_child_cty_plain <- function(){
  # clean child pop (5-19 yo) data at county level -- fips, all years, pop (for write_loess_fits_ILIn_age_cty.R)
  print(match.call())
  
  # import population data from mysql
  con <- dbConnect(RMySQL::MySQL(), group = "rmysql-fludrivers")
  dbListTables(con)
  
  dbListFields(con, "demog_Census_agePop_county")
  # sel.statement <- "Select * from demog_Census_agePop_county limit 5"
  sel.statement <- "SELECT fips, county, agegroup, year, pop FROM demog_Census_agePop_county WHERE scale = 'county' and agegroup = 'child'"
  dummy <- dbGetQuery(con, sel.statement)
  
  dbDisconnect(con)
  
  # clean final dataset
  output <- tbl_df(dummy) %>% 
    select(fips, year, pop)
  return(output)
}

################################
