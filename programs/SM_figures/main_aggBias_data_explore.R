## Name: Elizabeth Lee
## Date: 7/10/17
## Function: main code to analyze aggregation bias from model outcomes for wksToEpi
## Filenames: 
## Data Source: 
## Notes: 
################################

require(tidyverse)
setwd(dirname(sys.frame(1)$ofile))
source("source_aggBias_data_explore_functions.R")

#### set these! ###############################
dbCodeStr <- "_irDt_Octfit_span0.4_degree2"
modules <- c("scatter_burden_aggBias") # "statistics", "lisa_oneSeason", "lisa_allSeasons", scatterplot", "choro", "choroAvg", "choroHotspot", "choroAvgHotspot", dataExport", "correlog_allSeasons", "scatter_burden_aggBias"

###############################
## PATHS ##
setwd('../reference_data')
path_abbr_st <- paste0(getwd(), "/state_abbreviations_FIPS.csv")
path_latlon_cty <- paste0(getwd(), "/cty_pop_latlon.csv")
path_latlon_st <- paste0(getwd(), "/state_latlon.csv")
path_latlon_reg <- paste0(getwd(), "/region_latlon.csv")
path_region_cw <- paste0(getwd(), "/state_abbreviations_FIPS_region.csv")

setwd("../R_export")
path_response_cty <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_cty.csv", dbCodeStr))
path_response_st <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_st.csv", dbCodeStr))
path_response_reg <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_reg.csv", dbCodeStr))
path_list <- list(path_abbr_st = path_abbr_st, 
                  path_region_cw = path_region_cw,
                  path_latlon_st = path_latlon_st,
                  path_latlon_cty = path_latlon_cty,
                  path_latlon_reg = path_latlon_reg,
                  path_response_st = path_response_st,
                  path_response_reg = path_response_reg,
                  path_response_cty = path_response_cty)

setwd("./test_dbVariance_aggBiasMagnitude")
path_exportData1 <- getwd()

################################
## MAIN ##
setwd(dirname(sys.frame(1)$ofile))

#### import aggBias data ############################
offsetSetting <- FALSE
dataParams <- list(offset_l = offsetSetting, filepathList = path_list)

obs_wksToEpi_ctySt <- do.call(import_obs_wksToEpi_ctySt, c(dataParams))
# obs_wksToEpi_ctyReg <- do.call(import_obs_wksToEpi_ctyReg, c(dataParams))
# obs_wksToPeak_ctySt <- do.call(import_obs_wksToPeak_ctySt, c(dataParams))
# obs_wksToPeak_ctyReg <- do.call(import_obs_wksToPeak_ctyReg, c(dataParams))

# obs_iliEarly_ctySt <- do.call(import_obs_iliEarly_ctySt, c(dataParams))
# obs_iliEarly_ctyReg <- do.call(import_obs_iliEarly_ctyReg, c(dataParams))
# obs_iliPeak_ctySt <- do.call(import_obs_iliPeak_ctySt, c(dataParams))
# obs_iliPeak_ctyReg <- do.call(import_obs_iliPeak_ctyReg, c(dataParams))

#### statistics ############################
if("statistics" %in% modules){
  # early and peak comparisons
  timeSt <- pairedTest_aggBias_timingMagnitude(obs_wksToEpi_ctySt, obs_wksToPeak_ctySt)
  timeReg <- pairedTest_aggBias_timingMagnitude(obs_wksToEpi_ctyReg, obs_wksToPeak_ctyReg)
  magSt <- pairedTest_aggBias_timingMagnitude(obs_iliEarly_ctySt, obs_iliPeak_ctySt)
  magReg <- pairedTest_aggBias_timingMagnitude(obs_iliEarly_ctyReg, obs_iliPeak_ctyReg)

  # scale comparisons
  wksToEpiScale <- pairedTest_aggBias_spatialScales(obs_wksToEpi_ctySt, obs_wksToEpi_ctyReg)
  wksToPeakScale <- pairedTest_aggBias_spatialScales(obs_wksToPeak_ctySt, obs_wksToPeak_ctyReg)
  iliEarlyScale <- pairedTest_aggBias_spatialScales(obs_iliEarly_ctySt, obs_iliEarly_ctyReg)
  iliPeakScale <- pairedTest_aggBias_spatialScales(obs_iliPeak_ctySt, obs_iliPeak_ctyReg)
  # list(histPlot=histPlot, absHistPlot=absHistPlot, dbPlot=dbPlot, absDbPlot=absDbPlot, ttest=ttest, absTtest=absTtest)
}
############################
if("correlog_allSeasons" %in% modules){
  staticFormats <- list(w = 6, h = 4, dataProcess = "irDt", incrementKm = 10, resamp = 500)
  dynFormatLs <- data.frame(dataProcess = staticFormats$dataProcess, measure = c(rep("wksToEpi", 2), rep("wksToPeak", 2), rep("iliEarly", 2), rep("iliPeak", 2)), scaleDiff = rep(c("stCty", "regCty"), 4)) 
  dynFormats <- split(dynFormatLs, seq(nrow(dynFormatLs)))

  correlogStat_aggBias_allSeasons(obs_wksToEpi_ctySt, staticFormats, dynFormats[[1]])
  correlogStat_aggBias_allSeasons(obs_wksToEpi_ctyReg, staticFormats, dynFormats[[2]])
  correlogStat_aggBias_allSeasons(obs_wksToPeak_ctySt, staticFormats, dynFormats[[3]])
  correlogStat_aggBias_allSeasons(obs_wksToPeak_ctyReg, staticFormats, dynFormats[[4]])
  correlogStat_aggBias_allSeasons(obs_iliEarly_ctySt, staticFormats, dynFormats[[5]])
  correlogStat_aggBias_allSeasons(obs_iliEarly_ctyReg, staticFormats, dynFormats[[6]])
  correlogStat_aggBias_allSeasons(obs_iliPeak_ctySt, staticFormats, dynFormats[[7]])
  correlogStat_aggBias_allSeasons(obs_iliPeak_ctyReg, staticFormats, dynFormats[[8]])
}
############################
if("lisa_oneSeason" %in% modules){
  plotFormatsDf <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    scaleDiff = rep(c("stCty", "regCty"), 4), 
    neighSize = 160, # km units
    resamp = 500,
    w = 6, 
    h = 4)) %>%
    mutate(pltVar = paste0("obs_diff_", scaleDiff))

  lisa_obs_wksToEpi_ctySt_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_wksToEpi_ctySt, as.list(plotFormatsDf[1,]))
  lisa_obs_wksToEpi_ctyReg_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_wksToEpi_ctyReg, as.list(plotFormatsDf[2,]))
  lisa_obs_wksToPeak_ctySt_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_wksToPeak_ctySt, as.list(plotFormatsDf[3,]))
  lisa_obs_wksToPeak_ctyReg_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_wksToPeak_ctyReg, as.list(plotFormatsDf[4,]))
  lisa_obs_iliEarly_ctySt_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_iliEarly_ctySt, as.list(plotFormatsDf[5,]))
  lisa_obs_iliEarly_ctyReg_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_iliEarly_ctyReg, as.list(plotFormatsDf[6,]))
  lisa_obs_iliPeak_ctySt_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_iliPeak_ctySt, as.list(plotFormatsDf[7,]))
  lisa_obs_iliPeak_ctyReg_oneSeas <- lisa_aggBias_timingMagnitude_oneSeason(obs_iliPeak_ctyReg, as.list(plotFormatsDf[8,]))
}
############################
if("lisa_allSeasons" %in% modules){
  plotFormatsDf <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    scaleDiff = rep(c("stCty", "regCty"), 4), 
    neighSize = 1000, # km units
    resamp = 500,
    w = 6, 
    h = 4)) %>%
    mutate(pltVar = paste0("obs_diff_", scaleDiff))

  lisa_obs_wksToEpi_ctySt <- lisa_aggBias_timingMagnitude_allSeasons(obs_wksToEpi_ctySt, as.list(plotFormatsDf[1,]))
  lisa_obs_wksToEpi_ctyReg <- lisa_aggBias_timingMagnitude_allSeasons(obs_wksToEpi_ctyReg, as.list(plotFormatsDf[2,]))
  lisa_obs_wksToPeak_ctySt <- lisa_aggBias_timingMagnitude_allSeasons(obs_wksToPeak_ctySt, as.list(plotFormatsDf[3,]))
  lisa_obs_wksToPeak_ctyReg <- lisa_aggBias_timingMagnitude_allSeasons(obs_wksToPeak_ctyReg, as.list(plotFormatsDf[4,]))
  lisa_obs_iliEarly_ctySt <- lisa_aggBias_timingMagnitude_allSeasons(obs_iliEarly_ctySt, as.list(plotFormatsDf[5,]))
  lisa_obs_iliEarly_ctyReg <- lisa_aggBias_timingMagnitude_allSeasons(obs_iliEarly_ctyReg, as.list(plotFormatsDf[6,]))
  lisa_obs_iliPeak_ctySt <- lisa_aggBias_timingMagnitude_allSeasons(obs_iliPeak_ctySt, as.list(plotFormatsDf[7,]))
  lisa_obs_iliPeak_ctyReg <- lisa_aggBias_timingMagnitude_allSeasons(obs_iliPeak_ctyReg, as.list(plotFormatsDf[8,]))
}

#### plots ############################
# scatterplot of state versus county/region
if("scatterplot" %in% modules){
  staticFormats <- list(w = 6, h = 4, offset_l = FALSE)
  dynFormatLs <- data.frame(measure = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), bigscale = rep(c("st", "reg"), 4))
  dynFormats <- split(dynFormatLs, seq(nrow(dynFormatLs)))

  scatter_obsCompare_aggBias_timingMagnitude(obs_wksToEpi_ctySt, staticFormats, dynFormats[[1]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_wksToEpi_ctyReg, staticFormats, dynFormats[[2]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_wksToPeak_ctySt, staticFormats, dynFormats[[3]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_wksToPeak_ctyReg, staticFormats, dynFormats[[4]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_iliEarly_ctySt, staticFormats, dynFormats[[5]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_iliEarly_ctyReg, staticFormats, dynFormats[[6]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_iliPeak_ctySt, staticFormats, dynFormats[[7]])
  scatter_obsCompare_aggBias_timingMagnitude(obs_iliPeak_ctyReg, staticFormats, dynFormats[[8]])
}
############################
# choropleth of magnitude of aggregation bias - one season
if("choro" %in% modules){
  plotFormatsDf <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    scaleDiff = rep(c("stCty", "regCty"), 4), 
    w = 6, 
    h = 4)) %>%
    mutate(pltVar = paste0("obs_diff_", scaleDiff))

  if(grepl("ilinDt", dbCodeStr)){
    breaksDf <- data.frame(
      wksToEpi_ctySt = c(-21, -8, -3, -1, 1, 3, 8, 16), 
      wksToEpi_ctyReg = c(-25, -15, -4, -1, 1, 4, 15),
      wksToPeak_ctySt = c(-20, -8, -4, -1, 1, 4, 8, 18),
      wksToPeak_ctyReg = c(-18, -8, -4, -1, 1, 4, 8, 20),
      iliEarly_ctySt = c(-80, -15, -5, -1, 1, 5, 25),
      iliEarly_ctyReg = c(-83, -13, -4, -1, 1, 4, 13),
      iliPeak_ctySt = c(-52, -14, -4, -1, 1, 4, 14),
      iliPeak_ctyReg = c(-60, -10, -4, -1, 1, 4, 10)
    )
    paletteDf <- data.frame(
      wksToEpi_ctySt = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41", "#09622a"), 
      wksToEpi_ctyReg = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41"),
      wksToPeak_ctySt = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41", "#09622a"),
      wksToPeak_ctyReg = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41", "#09622a"),
      iliEarly_ctySt = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41"),
      iliEarly_ctyReg = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41"),
      iliPeak_ctySt = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41"),
      iliPeak_ctyReg = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41")
    )
  } else if(grepl("irDt", dbCodeStr)){
    breaksDf <- data.frame(
      wksToEpi_ctySt = c(-22, -8, -3, -1, 1, 3, 8, 18), 
      wksToEpi_ctyReg = c(-22, -8, -3,-1, 1, 3, 8, 16),
      wksToPeak_ctySt = c(-20, -8, -3, -1, 1, 3, 8, 18),
      wksToPeak_ctyReg = c(-18, -8, -3, -1, 1, 3, 8, 18),
      iliEarly_ctySt = c(-22, -4, -2, -1, 1, 2, 4, 8),
      iliEarly_ctyReg = c(-26, -4, -2, -1, 1, 2, 4, 8),
      iliPeak_ctySt = c(-14, -4, -2, -1, 1, 2, 4, 6),
      iliPeak_ctyReg = c(-16, -4, -3, -1, 1, 2, 4, 6)
    )
    paletteDf <- data.frame(
      breaks8 = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41", "#09622a")
    )
  }

  choro_obs_aggBias_oneSeason_wrapper(obs_wksToEpi_ctySt, as.list(plotFormatsDf[1,]), breaksDf[,1], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_wksToEpi_ctyReg, as.list(plotFormatsDf[2,]), breaksDf[,2], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_wksToPeak_ctySt, as.list(plotFormatsDf[3,]), breaksDf[,3], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_wksToPeak_ctyReg, as.list(plotFormatsDf[4,]), breaksDf[,4], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_iliEarly_ctySt, as.list(plotFormatsDf[5,]), breaksDf[,5], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_iliEarly_ctyReg, as.list(plotFormatsDf[6,]), breaksDf[,6], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_iliPeak_ctySt, as.list(plotFormatsDf[7,]), breaksDf[,7], paletteDf[,1])
  choro_obs_aggBias_oneSeason_wrapper(obs_iliPeak_ctyReg, as.list(plotFormatsDf[8,]), breaksDf[,8], paletteDf[,1])

  # choro_obs_aggBias_stCty_wksToEpi_oneSeason(obs_wksToEpi_ctySt, plotFormats)
  # choro_obs_aggBias_regCty_wksToEpi_oneSeason(obs_wksToEpi_ctyReg, plotFormats)
  # choro_obs_aggBias_stCty_wksToPeak_oneSeason(obs_wksToPeak_ctySt, plotFormats)
  # choro_obs_aggBias_regCty_wksToPeak_oneSeason(obs_wksToPeak_ctyReg, plotFormats)
  # choro_obs_aggBias_stCty_iliEarly_oneSeason(obs_iliEarly_ctySt, plotFormats)
  # choro_obs_aggBias_regCty_iliEarly_oneSeason(obs_iliEarly_ctyReg, plotFormats)
  # choro_obs_aggBias_stCty_iliPeak_oneSeason(obs_iliPeak_ctySt, plotFormats)
  # choro_obs_aggBias_regCty_iliPeak_oneSeason(obs_iliPeak_ctyReg, plotFormats)

  # Interpretation: positive error (green) means that state model predicted a later epidemic onset than the county model
}
############################
# choropleth of magnitude of aggregation bias hotspots (identified with high spatial autocorrelation) - one season
if("choroHotspot" %in% modules){
  plotFormatsDf <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    scaleDiff = rep(c("stCty", "regCty"), 4), 
    w = 6, 
    h = 4)) 

  breaksDf <- data.frame(
    wksToEpi_ctySt = c(-22, -8, -3, -1, 1, 3, 8, 18), 
    wksToEpi_ctyReg = c(-22, -8, -3,-1, 1, 3, 8, 16),
    wksToPeak_ctySt = c(-20, -8, -3, -1, 1, 3, 8, 18),
    wksToPeak_ctyReg = c(-18, -8, -3, -1, 1, 3, 8, 18),
    iliEarly_ctySt = c(-22, -4, -2, -1, 1, 2, 4, 8),
    iliEarly_ctyReg = c(-26, -4, -2, -1, 1, 2, 4, 8),
    iliPeak_ctySt = c(-14, -4, -2, -1, 1, 2, 4, 6),
    iliPeak_ctyReg = c(-16, -4, -3, -1, 1, 2, 4, 6)
  )
  paletteDf <- data.frame(
    breaks8 = c("#10456a", "#1c73b1", "#67add4", "#cacaca", "#69a761", "#2f8e41", "#09622a") # 
  )
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[1,]), breaksDf[,1], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[2,]), breaksDf[,2], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[3,]), breaksDf[,3], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[4,]), breaksDf[,4], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[5,]), breaksDf[,5], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[6,]), breaksDf[,6], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[7,]), breaksDf[,7], paletteDf[,1])
  choro_obs_aggBias_hotspots_oneSeason_wrapper(as.list(plotFormatsDf[8,]), breaksDf[,8], paletteDf[,1])
}
############################
# choropleth of magnitude of aggregation bias - average seasons
if("choroAvg" %in% modules){
  plotFormatsDf <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    scaleDiff = rep(c("stCty", "regCty"), 4), 
    # breaks = ,
    w = 6, 
    h = 4)) %>%
    mutate(pltVar = paste0("obs_diff_", scaleDiff))

  if(grepl("ilinDt", dbCodeStr)){
    breaksDf <- data.frame(
      wksToEpi_ctySt = c(-14, -5, -3, -1, 1, 3, 5, 8), 
      wksToEpi_ctyReg = c(-15, -6, -3,-1, 1, 3, 6, 9),
      wksToPeak_ctySt = c(-10, -4, -2, -1, 1, 2, 4, 7),
      wksToPeak_ctyReg = c(-8, -3, -2, -1, 1, 2, 3, 7),
      iliEarly_ctySt = c(-20, -4, -2, -1, 1, 2, 4, 5),
      iliEarly_ctyReg = c(-21, -4, -2, -1, 1, 2, 4, 6),
      iliPeak_ctySt = c(-31, -4, -2, -1, 1, 2, 4, 6),
      iliPeak_ctyReg = c(-34, -6, -3, -1, 1, 2, 3.5, 5))
  } else if(grepl("irDt", dbCodeStr)){
    # breaks need to be updated for irDt?
    breaksDf <- data.frame(
      wksToEpi_ctySt = c(-14, -5, -3, -1, 1, 3, 5, 14), 
      wksToEpi_ctyReg = c(-9, -6, -3,-1, 1, 3, 6, 9),
      wksToPeak_ctySt = c(-7, -3, -2, -1, 1, 2, 3, 7),
      wksToPeak_ctyReg = c(-7, -3, -2, -1, 1, 2, 3, 7),
      iliEarly_ctySt = c(-20, -4, -2, -1, 1, 2, 4, 5),
      iliEarly_ctyReg = c(-21, -4, -2, -1, 1, 2, 4, 6),
      iliPeak_ctySt = c(-31, -4, -2, -1, 1, 2, 4, 6),
      iliPeak_ctyReg = c(-34, -6, -3, -1, 1, 2, 3.5, 5))
  }
 
  choro_obs_aggBias_avgSeason(obs_wksToEpi_ctySt, as.list(plotFormatsDf[1,]), breaksDf$wksToEpi_ctySt)
  choro_obs_aggBias_avgSeason(obs_wksToEpi_ctyReg, as.list(plotFormatsDf[2,]), breaksDf$wksToEpi_ctyReg)
  choro_obs_aggBias_avgSeason(obs_wksToPeak_ctySt, as.list(plotFormatsDf[3,]), breaksDf$wksToPeak_ctySt)
  choro_obs_aggBias_avgSeason(obs_wksToPeak_ctyReg, as.list(plotFormatsDf[4,]), breaksDf$wksToPeak_ctyReg)
  choro_obs_aggBias_avgSeason(obs_iliEarly_ctySt, as.list(plotFormatsDf[5,]), breaksDf$iliEarly_ctySt)
  choro_obs_aggBias_avgSeason(obs_iliEarly_ctyReg, as.list(plotFormatsDf[6,]), breaksDf$iliEarly_ctyReg)
  choro_obs_aggBias_avgSeason(obs_iliPeak_ctySt, as.list(plotFormatsDf[7,]), breaksDf$iliPeak_ctySt)
  choro_obs_aggBias_avgSeason(obs_iliPeak_ctyReg, as.list(plotFormatsDf[8,]), breaksDf$iliPeak_ctyReg)


  # Interpretation: positive error (green) means that state model predicted a later epidemic onset than the county model
}
############################
# choropleth of magnitude of aggregation bias - average seasons
if("choroAvgHotspot" %in% modules){
  plotFormatsDf <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    scaleDiff = rep(c("stCty", "regCty"), 4), 
    # breaks = ,
    w = 6, 
    h = 4)) %>%
    mutate(pltVar = paste0("obs_diff_", scaleDiff))
  
  breaksDf <- data.frame(
    wksToEpi_ctySt = c(-14, -5, -3, -1, 1, 3, 5, 14), 
    wksToEpi_ctyReg = c(-9, -6, -3,-1, 1, 3, 6, 9),
    wksToPeak_ctySt = c(-7, -3, -2, -1, 1, 2, 3, 7),
    wksToPeak_ctyReg = c(-7, -3, -2, -1, 1, 2, 3, 7),
    iliEarly_ctySt = c(-20, -4, -2, -1, 1, 2, 4, 5),
    iliEarly_ctyReg = c(-21, -4, -2, -1, 1, 2, 4, 6),
    iliPeak_ctySt = c(-31, -4, -2, -1, 1, 2, 4, 6),
    iliPeak_ctyReg = c(-34, -6, -3, -1, 1, 2, 3.5, 5))
 
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[1,]), breaksDf$wksToEpi_ctySt)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[2,]), breaksDf$wksToEpi_ctyReg)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[3,]), breaksDf$wksToPeak_ctySt)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[4,]), breaksDf$wksToPeak_ctyReg)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[5,]), breaksDf$iliEarly_ctySt)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[6,]), breaksDf$iliEarly_ctyReg)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[7,]), breaksDf$iliPeak_ctySt)
  choro_obs_aggBias_hotspots_avgSeason(as.list(plotFormatsDf[8,]), breaksDf$iliPeak_ctyReg)


  # Interpretation: positive error (green) means that state model predicted a later epidemic onset than the county model
}
############################
# export data to file
if("dataExport" %in% modules){
  # data export for explore_dbVariance_aggBiasMagnitude.R
  dataFormats <- tbl_df(data.frame(
    dbCode = c(rep("irDt_wksToEpi", 2), rep("irDt_wksToPeak", 2), rep("irDt_iliEarly", 2), rep("irDt_iliPeak", 2)), 
    pltVar = "obs_diff_stCty")) %>%
    mutate(exportPath = paste0(path_exportData1, "/aggBiasDiff_yVariance_", rep(c("st_", "cty_"), 4), dbCode, ".csv")) 

  write_st_aggBias_yVariance_mean(obs_wksToEpi_ctySt, as.list(dataFormats[1,]))
  # write_cty_aggBias_mean(obs_wksToEpi_ctySt, as.list(dataFormats[2,]))
  write_st_aggBias_yVariance_mean(obs_wksToPeak_ctySt, as.list(dataFormats[3,]))
  # write_cty_aggBias_mean(obs_wksToPeak_ctySt, as.list(dataFormats[4,]))
  write_st_aggBias_yVariance_mean(obs_iliEarly_ctySt, as.list(dataFormats[5,]))
  # write_cty_aggBias_mean(obs_iliEarly_ctySt, as.list(dataFormats[6,]))
  write_st_aggBias_yVariance_mean(obs_iliPeak_ctySt, as.list(dataFormats[7,]))
  # write_cty_aggBias_mean(obs_iliPeak_ctySt, as.list(dataFormats[8,]))

  # cty aggBias doesn't answer the correct question

}
############################
# export data to file
if("scatter_burden_aggBias" %in% modules){
   add code to import plot formats
  scatter_burden_vs_aggBias(obs_wksToEpi_ctySt)
  scatter_burden_vs_aggBias(obs_wksToPeak_ctySt)
  scatter_burden_vs_aggBias(obs_iliEarly_ctySt)
  scatter_burden_vs_aggBias(obs_iliPeak_ctySt)
}  