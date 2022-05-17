## Name: Elizabeth Lee
## Date: 10/19/17
## Function: 
## Data Source: IMS Health
## Notes: 
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### header #################################
rm(list = ls())
require(tidyverse)
require(INLA)

#### SOURCE: clean and import model data #################################
source("programs/source_import_modeldata.R")
source("programs/functions_tab1_SMtables.R")

#### FILEPATHS #################################
setwd('reference_data')
path_abbr_st <- paste0(getwd(), "/state_abbreviations_FIPS.csv")
path_latlon_cty <- paste0(getwd(), "/cty_pop_latlon.csv")
path_latlon_st <- paste0(getwd(), "/state_latlon.csv")
path_latlon_reg <- paste0(getwd(), "/region_latlon.csv")
path_region_cw <- paste0(getwd(), "/state_abbreviations_FIPS_region.csv")
path_adjMxExport_cty <- paste0(getwd(), "/UScounty_shapefiles/US_county_adjacency.graph")
path_graphIdx_cty <- paste0(getwd(), "/UScounty_shapefiles/US_county_graph_index.csv")


#### MODEL SETUP ###############################
formula_CAR <- Y ~ 1 + f(graphIdx, model = "besag", graph = path_adjMxExport_cty)
formula_iid <- Y ~ 1 + f(graphIdx, model = "iid")

st_datFormats <- list(variable1 = "obs_diff_stCty", variable2 = "obs_diff_stCty")
reg_datFormats <- list(variable1 = "obs_diff_regCty", variable2 = "obs_diff_regCty")
scale_datFormats <- list(variable1 = "obs_diff_regCty", variable2 = "obs_diff_stCty")

#### set these! ###############################
dbCodeStr <- "_irDt_Octfit_span0.4_degree2"
formula <- formula_iid
formulaType <- "iid"

#### MODEL INPUTS ###############################
setwd("R_export")
path_response_cty <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_cty.csv", dbCodeStr))
path_response_st <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_st.csv", dbCodeStr))
path_response_reg <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_reg.csv", dbCodeStr))

path_list <- list(path_abbr_st = path_abbr_st, 
                  path_graphIdx_cty = path_graphIdx_cty,
                  path_region_cw = path_region_cw,
                  path_latlon_st = path_latlon_st,
                  path_latlon_cty = path_latlon_cty,
                  path_latlon_reg = path_latlon_reg,
                  path_response_st = path_response_st,
                  path_response_reg = path_response_reg,
                  path_response_cty = path_response_cty)


#### IMPORT DATA ############################
offsetSetting <- FALSE
dataParams <- list(offset_l = offsetSetting, filepathList = path_list)

obs_wksToEpi_ctySt <- do.call(import_obs_wksToEpi_ctySt, c(dataParams))
obs_wksToEpi_ctyReg <- do.call(import_obs_wksToEpi_ctyReg, c(dataParams))
obs_wksToPeak_ctySt <- do.call(import_obs_wksToPeak_ctySt, c(dataParams))
obs_wksToPeak_ctyReg <- do.call(import_obs_wksToPeak_ctyReg, c(dataParams))

obs_iliEarly_ctySt <- do.call(import_obs_iliEarly_ctySt, c(dataParams))
obs_iliEarly_ctyReg <- do.call(import_obs_iliEarly_ctyReg, c(dataParams))
obs_iliPeak_ctySt <- do.call(import_obs_iliPeak_ctySt, c(dataParams))
obs_iliPeak_ctyReg <- do.call(import_obs_iliPeak_ctyReg, c(dataParams))

#### RUN MODELS #################################
#### Is aggregation bias more prevalent in early season or peak season measures? ################################# 
## State Timing ##
timing_st_data <- create_response_data(obs_wksToPeak_ctySt, obs_wksToEpi_ctySt, st_datFormats, path_list)
timing_st <- run_intercept_model(timing_st_data, formula)

## State Magnitude ##
mag_st_data <- create_response_data(obs_iliPeak_ctySt, obs_iliEarly_ctySt, st_datFormats, path_list)
mag_st <- run_intercept_model(mag_st_data, formula)

## Region Timing ##
timing_reg_data <- create_response_data(obs_wksToPeak_ctyReg, obs_wksToEpi_ctyReg, reg_datFormats, path_list)
timing_reg <- run_intercept_model(timing_reg_data, formula)

## Region Magnitude ##
mag_reg_data <- create_response_data(obs_iliPeak_ctyReg, obs_iliEarly_ctyReg, reg_datFormats, path_list)
mag_reg <- run_intercept_model(mag_reg_data, formula)


#### Does the magnitude of aggregation bias increase as the scale of aggregation increases from state to region? ################################# 
## Weeks to Epi ##
wksToEpi_data <- create_response_data(obs_wksToEpi_ctyReg, obs_wksToEpi_ctySt, scale_datFormats, path_list)
wksToEpi <- run_intercept_model(wksToEpi_data, formula)

## Weeks to Peak ##
wksToPeak_data <- create_response_data(obs_wksToPeak_ctyReg, obs_wksToPeak_ctySt, scale_datFormats, path_list)
wksToPeak <- run_intercept_model(wksToPeak_data, formula)

## Early ILI ##
earlyILI_data <- create_response_data(obs_iliEarly_ctyReg, obs_iliEarly_ctySt, scale_datFormats, path_list)
earlyILI <- run_intercept_model(earlyILI_data, formula)

## Peak ILI ##
peakILI_data <- create_response_data(obs_iliPeak_ctyReg, obs_iliPeak_ctySt, scale_datFormats, path_list)
peakILI <- run_intercept_model(peakILI_data, formula)


#### MODEL OUTPUTS #################################
setwd("./inla_aggBias_comparisons")
exportPath_cpo1 <- paste0(getwd(), "/modFit_timingMagnitude_", formulaType, ".csv")
exportPath_summ1 <- paste0(getwd(), "/summaryStats_timingMagnitude_", formulaType, ".csv")
exportPath_cpo2 <- paste0(getwd(), "/modFit_spatialScales_", formulaType, ".csv")
exportPath_summ2 <- paste0(getwd(), "/summaryStats_spatialScales_", formulaType, ".csv")


#### EXPORT DATA #################################
models1 <- list(timing_st = timing_st, mag_st = mag_st, timing_reg = timing_reg, mag_reg = mag_reg)
models2 <- list(wksToEpi = wksToEpi, wksToPeak = wksToPeak, earlyILI = earlyILI, peakILI = peakILI)

cpo1_ls <- lapply(models1, grab_dicCPO, formulaType = formulaType)
summ1_ls <- lapply(models1, grab_summaryStats, formulaType = formulaType)
cpo2_ls <- lapply(models2, grab_dicCPO, formulaType = formulaType)
summ2_ls <- lapply(models2, grab_summaryStats, formulaType = formulaType)

cpo1 <- rbindlist(cpo1_ls, use.names = TRUE, idcol = TRUE)
summ1 <- rbindlist(summ1_ls, use.names = TRUE, idcol = TRUE)
cpo2 <- rbindlist(cpo2_ls, use.names = TRUE, idcol = TRUE)
summ2 <- rbindlist(summ2_ls, use.names = TRUE, idcol = TRUE)

write_csv(cpo1, exportPath_cpo1)
write_csv(summ1, exportPath_summ1)
write_csv(cpo2, exportPath_cpo2)
write_csv(summ2, exportPath_summ2)