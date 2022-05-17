## Name: Elizabeth Lee
## Date: 11/14/19
## Function: main code to identify top 1% most probably flu season source (origin) location at the county level; export dataframes for our version and Charu2017 supplementary table
## Filenames: physicianCoverage_IMSHealth_state.csv, dbMetrics_periodicReg_ilinDt_Octfit_span0.4_degree2_analyzeDB_st.csv
## Data Source: IMS Health
## Notes: 
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### header #################################
rm(list = ls())
require(tidyverse) # clean_data_functions dependencies


#### set these! ################################

#### SOURCE: clean and import response data #################################
source("programs/functions_identify_origin.R")

#### FILEPATHS #################################


#### MAIN ############################
wksToEpiDat <- read_csv("R_export/inlaModelData_import/inlaImport_model10f_wksToEpi_irDt_v7.csv")
earlyLocDat <- subset_earliest_onset_locations_decile(wksToEpiDat) %>%
  rename(srcFips = fips_st)

seasons <- wksToEpiDat %>% distinct(season) %>% arrange(season) %>% unlist 
corrDat <- data.frame()

for (s in seasons){
  wksToEpiDat_seas <- wksToEpiDat %>% filter(season == s)
  earlyLocDat_seas <- earlyLocDat %>% filter(season == s)

  prepDat <- calculate_distance_from_sources(earlyLocDat_seas, wksToEpiDat_seas)

  dummyDat <- calculate_correlation_distance_onsetWeek(prepDat)
  dummyDf <- data.frame(season = s, srcID = names(dummyDat), corrCoef = t(dummyDat))
  corrDat <- tbl_df(bind_rows(corrDat, dummyDf))
}

fullDat <- left_join(earlyLocDat %>% select(srcID, srcFips), corrDat, by = "srcID")

srcLocDat_corr <- fullDat %>% 
  group_by(season) %>% 
  dplyr::top_n(2, corrCoef) %>% 
  ungroup %>%
  left_join(earlyLocDat %>% distinct(srcFips, srcLat, srcLon), by = c("srcFips")) %>% 
  dplyr::select(season, srcFips, srcLat, srcLon, corrCoef)

write_csv(srcLocDat_corr, "R_export/origin_locations/fluseason_source_states.csv")
