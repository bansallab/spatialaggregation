# 9/28/17
# Function: Plot response data in choropleths - by season and average across seasons
rm(list = ls())
require(tidyverse)
require(data.table)

setwd(dirname(sys.frame(1)$ofile))
source("source_import_modeldata.R")
source("source_response_data_explore.R")
####################################
# set these
dbCodeStr <- "_irDt_Octfit_span0.4_degree2"
w <- 5.5; h = 4
offset <- FALSE
incrementKm <- 10; resamp <- 500

####################################
# paths
setwd('../reference_data')
path_abbr_st <- paste0(getwd(), "/state_abbreviations_FIPS.csv")
path_latlon_cty <- paste0(getwd(), "/cty_pop_latlon.csv")
path_latlon_st <- paste0(getwd(), "/state_latlon.csv")
path_latlon_reg <- paste0(getwd(), "/region_latlon.csv")
path_region_cw <- paste0(getwd(), "/state_abbreviations_FIPS_region.csv")

setwd("../R_export")
path_response_cty <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_cty.csv", dbCodeStr))
path_fullIndic_cty <- paste0(getwd(), sprintf("/fullIndicAll_periodicReg%s_analyzeDB_cty.csv", dbCodeStr))
path_response_st <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_st.csv", dbCodeStr))
path_fullIndic_st <- paste0(getwd(), sprintf("/fullIndicAll_periodicReg%s_analyzeDB_st.csv", dbCodeStr))
path_response_reg <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_reg.csv", dbCodeStr))
path_fullIndic_reg <- paste0(getwd(), sprintf("/fullIndicAll_periodicReg%s_analyzeDB_reg.csv", dbCodeStr))


path_list <- list(path_abbr_st = path_abbr_st,
                  path_latlon_cty = path_latlon_cty,
                  path_region_cw = path_region_cw,
                  path_response_cty = path_response_cty,
                  path_fullIndic_cty = path_fullIndic_cty, 
                  path_latlon_st = path_latlon_st,
                  path_response_st = path_response_st,
                  path_fullIndic_st = path_fullIndic_st, 
                  path_latlon_reg = path_latlon_reg,
                  path_response_reg = path_response_reg,
                  path_fullIndic_reg = path_fullIndic_reg)


setwd("../graph_outputs/response_data_explore")
exportPath <- getwd()

# import cty & st data
wksToEpi_ctySt <- import_obs_wksToEpi_ctySt(offset, path_list)
wksToPeak_ctySt <- import_obs_wksToPeak_ctySt(offset, path_list)
iliEarly_ctySt <- import_obs_iliEarly_ctySt(offset, path_list)
iliPeak_ctySt <- import_obs_iliPeak_ctySt(offset, path_list)
# import cty & reg data (so it's plottable on a cty map)
wksToEpi_ctyReg <- import_obs_wksToEpi_ctyReg(offset, path_list)
wksToPeak_ctyReg <- import_obs_wksToPeak_ctyReg(offset, path_list)
iliEarly_ctyReg <- import_obs_iliEarly_ctyReg(offset, path_list)
iliPeak_ctyReg <- import_obs_iliPeak_ctyReg(offset, path_list)

# merge data for all three scales
wksToEpi <- merge_obs_ctyStReg(wksToEpi_ctySt, wksToEpi_ctyReg)
wksToPeak <- merge_obs_ctyStReg(wksToPeak_ctySt, wksToPeak_ctyReg)
iliEarly <- merge_obs_ctyStReg(iliEarly_ctySt, iliEarly_ctyReg)
iliPeak <- merge_obs_ctyStReg(iliPeak_ctySt, iliPeak_ctyReg)


#### MAIN CODE ################################
#### county data figures ####
format_wksToEpi <- list(w = w, h = h, measure = "wksToEpi", dataProcess = "irDt", legendStep = 5, offset_l = offset, exportPath = exportPath, dataScale = "cty", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(wksToEpi, format_wksToEpi))
# do.call(choro_obs_db_avgSeason, list(wksToEpi, format_wksToEpi))
correlog_wksToEpi_cty <- do.call(correlogStat_obs_allSeasons, list(wksToEpi, format_wksToEpi))

format_wksToPeak <- list(w = w, h = h, measure = "wksToPeak", dataProcess = "irDt", legendStep = 4, offset_l = offset, exportPath = exportPath, dataScale = "cty", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(wksToPeak, format_wksToPeak))
# do.call(choro_obs_db_avgSeason, list(wksToPeak, format_wksToPeak))
correlog_wksToPeak_cty <- do.call(correlogStat_obs_allSeasons, list(wksToPeak, format_wksToPeak))

format_iliEarly <- list(w = w, h = h, measure = "iliEarly", dataProcess = "irDt", legendStep = 0.5, offset_l = offset, exportPath = exportPath, dataScale = "cty", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(iliEarly, format_iliEarly))
# do.call(choro_obs_db_avgSeason, list(iliEarly, format_iliEarly))
correlog_iliEarly_cty <- do.call(correlogStat_obs_allSeasons, list(iliEarly, format_iliEarly))

format_iliPeak <- list(w = w, h = h, measure = "iliPeak", dataProcess = "irDt", legendStep = 0.5, offset_l = offset, exportPath = exportPath, dataScale = "cty", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(iliPeak, format_iliPeak))
# do.call(choro_obs_db_avgSeason, list(iliPeak, format_iliPeak))
correlog_iliPeak_cty <- do.call(correlogStat_obs_allSeasons, list(iliPeak, format_iliPeak))

################################
#### state data figures ####
format_wksToEpi <- list(w = w, h = h, measure = "wksToEpi", dataProcess = "irDt", legendStep = 5, offset_l = offset, exportPath = exportPath, dataScale = "st", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(wksToEpi, format_wksToEpi))
# do.call(choro_obs_db_avgSeason, list(wksToEpi, format_wksToEpi))
# correlog_wksToEpi_st <- do.call(correlogStat_obs_allSeasons, list(wksToEpi, format_wksToEpi))

format_wksToPeak <- list(w = w, h = h, measure = "wksToPeak", dataProcess = "irDt", legendStep = 4, offset_l = offset, exportPath = exportPath, dataScale = "st", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(wksToPeak, format_wksToPeak))
# do.call(choro_obs_db_avgSeason, list(wksToPeak, format_wksToPeak))
# correlog_wksToPeak_st <- do.call(correlogStat_obs_allSeasons, list(wksToPeak, format_wksToPeak))

format_iliEarly <- list(w = w, h = h, measure = "iliEarly", dataProcess = "irDt", legendStep = 50, offset_l = offset, exportPath = exportPath, dataScale = "st", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(iliEarly, format_iliEarly))
# do.call(choro_obs_db_avgSeason, list(iliEarly, format_iliEarly))
# correlog_iliEarly_st <- do.call(correlogStat_obs_allSeasons, list(iliEarly, format_iliEarly))

format_iliPeak <- list(w = w, h = h, measure = "iliPeak", dataProcess = "irDt", legendStep = 50, offset_l = offset, exportPath = exportPath, dataScale = "st", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(iliPeak, format_iliPeak))
# do.call(choro_obs_db_avgSeason, list(iliPeak, format_iliPeak))
# correlog_iliPeak_st <- do.call(correlogStat_obs_allSeasons, list(iliPeak, format_iliPeak))

################################
## region data figures ##
format_wksToEpi <- list(w = w, h = h, measure = "wksToEpi", dataProcess = "irDt", legendStep = 5, offset_l = offset, exportPath = exportPath, dataScale = "reg", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(wksToEpi, format_wksToEpi))
# do.call(choro_obs_db_avgSeason, list(wksToEpi, format_wksToEpi))
# correlog_wksToEpi_reg <- do.call(correlogStat_obs_allSeasons, list(wksToEpi, format_wksToEpi))

format_wksToPeak <- list(w = w, h = h, measure = "wksToPeak", dataProcess = "irDt", legendStep = 4, offset_l = offset, exportPath = exportPath, dataScale = "reg", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(wksToPeak, format_wksToPeak))
# do.call(choro_obs_db_avgSeason, list(wksToPeak, format_wksToPeak))
# correlog_wksToPeak_reg <- do.call(correlogStat_obs_allSeasons, list(wksToPeak, format_wksToPeak))

format_iliEarly <- list(w = w, h = h, measure = "iliEarly", dataProcess = "irDt", legendStep = 250, offset_l = offset, exportPath = exportPath, dataScale = "reg", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(iliEarly, format_iliEarly))
# do.call(choro_obs_db_avgSeason, list(iliEarly, format_iliEarly))
# correlog_iliEarly_reg <- do.call(correlogStat_obs_allSeasons, list(iliEarly, format_iliEarly))

format_iliPeak <- list(w = w, h = h, measure = "iliPeak", dataProcess = "irDt", legendStep = 500, offset_l = offset, exportPath = exportPath, dataScale = "reg", incrementKm = incrementKm, resamp = resamp)
# do.call(choro_obs_db_oneSeason, list(iliPeak, format_iliPeak))
# do.call(choro_obs_db_avgSeason, list(iliPeak, format_iliPeak))
# correlog_iliPeak_reg <- do.call(correlogStat_obs_allSeasons, list(iliPeak, format_iliPeak))