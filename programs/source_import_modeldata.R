## Name: Elizabeth Lee
## Date: 6/23/17
## Function: functions to import surveillance model data
## Filenames: 
## Data Source: 
## Notes: 
################################

require(tidyverse)
require(data.table)
library(maps)

source("programs/source_clean_response_functions_cty.R") # functions to clean response and IMS coverage data (cty)
source("programs/source_clean_response_functions_st.R") # functions to clean response (st)
source("programs/source_clean_response_functions_reg.R") # functions to clean response (reg)

#### group observed and fitted model data ################################
import_obs_allMeasures_cty <- function(filepathList){
  print(match.call())

  # import all 4 disease burden measures
  wksToEpiDat <- cleanR_wksToEpi_cty(filepathList) %>%
    mutate(y_wksToEpi = y1) %>%
    select(season, fips, y_wksToEpi)
  wksToPeakDat <- cleanR_wksToPeak_cty(filepathList) %>%
    mutate(y_wksToPeak = y1) %>%
    select(season, fips, y_wksToPeak)
  iliEarlyDat <- cleanR_iliEarly_irDt_shift1_cty(filepathList) %>%
    mutate(rr_iliEarly = y1/E) %>%
    select(season, fips, rr_iliEarly)
  iliPeakDat <- cleanR_iliPeak_irDt_shift1_cty(filepathList) %>%
    mutate(rr_iliPeak = y1/E) %>%
    select(season, fips, rr_iliPeak)

  # add lat/lon coords
  coordDat <- read_csv(filepathList$path_latlon_cty, col_types = "_c__dd")

  obsDat <- full_join(wksToEpiDat, wksToPeakDat, by = c("season", "fips")) %>%
    full_join(iliEarlyDat, by = c("season", "fips")) %>%
    full_join(iliPeakDat, by = c("season", "fips")) %>%
    left_join(coordDat, by = c("fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))

  return(obsDat)
}
################################
import_obs_aggBias_allMeasures <- function(filepathList, dataFormats){
  print(match.call())

  offset <- dataFormats$offset_l
  bigscale <- dataFormats$bigscale

  if(bigscale == "st"){
    wksToEpiDat <- import_obs_wksToEpi_ctySt(offset, filepathList) %>%
      rename(bias_wksToEpi = obs_diff_stCty) %>%
      select(season, fips, fips_st, latitude, longitude, bias_wksToEpi)
    wksToPeakDat <- import_obs_wksToPeak_ctySt(offset, filepathList) %>%
      rename(bias_wksToPeak = obs_diff_stCty) %>%
      select(season, fips, bias_wksToPeak)
    iliEarlyDat <- import_obs_iliEarly_ctySt(offset, filepathList) %>%
      rename(bias_iliEarly = obs_diff_stCty) %>%
      select(season, fips, bias_iliEarly)
    iliPeakDat <- import_obs_iliPeak_ctySt(offset, filepathList) %>%
      rename(bias_iliPeak = obs_diff_stCty) %>%
      select(season, fips, bias_iliPeak)
  } else if(bigscale == "reg"){
    wksToEpiDat <- import_obs_wksToEpi_ctyReg(offset, filepathList) %>%
      rename(bias_wksToEpi = obs_diff_regCty) %>%
      select(season, fips, fips_st, latitude, longitude, bias_wksToEpi)
    wksToPeakDat <- import_obs_wksToPeak_ctyReg(offset, filepathList) %>%
      rename(bias_wksToPeak = obs_diff_regCty) %>%
      select(season, fips, bias_wksToPeak)
    iliEarlyDat <- import_obs_iliEarly_ctyReg(offset, filepathList) %>%
      rename(bias_iliEarly = obs_diff_regCty) %>%
      select(season, fips, bias_iliEarly)
    iliPeakDat <- import_obs_iliPeak_ctyReg(offset, filepathList) %>%
      rename(bias_iliPeak = obs_diff_regCty) %>%
      select(season, fips, bias_iliPeak)
  }

  obsDat <- full_join(wksToEpiDat, wksToPeakDat, by = c("season", "fips")) %>%
    full_join(iliEarlyDat, by = c("season", "fips")) %>%
    full_join(iliPeakDat, by = c("season", "fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))

  return(obsDat)
}
################################
import_obs_correlog <- function(datFormats){
  print(match.call())

  dataScale <- datFormats$dataScale
  resamp <- datFormats$resamp

  fnames <- list.files(string_exportDat_response_data_folder(), pattern = paste0(dataScale, "_resamp", resamp), full.names = TRUE)
  correlogDat <- map_df(fnames, function(x){
    inDat <- read_csv(x)
    })

  return(correlogDat)
}
################################


#### wks.to.epi ################################
import_obsFit_wksToEpi <- function(modCodeStr, filepathList){
  print(match.call())
  # import observed and fitted data for weeks to epidemic onset at county level
  
  # import fitted data (on the scale of y)
  outDat <- read_csv(string_fit_fname(modCodeStr), col_types = "c_d_c_ddd_d___") %>%
    rename(fit_y = mean, fit_sd = sd) %>%
    select(modCodeStr, season, fips, fit_y, fit_sd, q_025, q_975)
  
  # import observed and expected wks to epi
  inDat <- cleanR_wksToEpi_cty(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips, obs_y, E)
  
  # prepare data for plotting breaks
  obsFitDat <- left_join(outDat, inDat, by = c("season", "fips")) %>%
    # mutate(obs_rr = obs_y/E, fit_rr = fit_y/E) %>%
    # mutate(cty_rr_q025 = q_025/E, cty_rr_q975 = q_975/E) %>%
    mutate(resid = (obs_y - fit_y)/fit_sd)
  
  return(obsFitDat)
}
################################
import_obsFit_wksToEpi_st <- function(modCodeStr, filepathList){
  print(match.call())
  # import observed and fitted data for weeks to epidemic onset at state level
  
  # import fitted data
  outDat <- read_csv(string_fit_fname(modCodeStr), col_types = "c_d_c_ddd_d___") %>%
    rename(fit_y = mean, fit_sd = sd) %>%
    select(modCodeStr, season, fips_st, fit_y, fit_sd, q_025, q_975)
  
  # import observed and expected wks to epi
  inDat <- cleanR_wksToEpi_st(filepathList) %>%
    mutate(obs_y = y1, E = E) %>%
    select(season, fips_st, obs_y, E)
  
  # prepare data for plotting breaks
  obsFitDat <- left_join(outDat, inDat, by = c("season", "fips_st")) %>%
    # mutate(obs_rr = obs_y/E, fit_rr = fit_y/E) %>% 
    # mutate(st_rr_q025 = q_025/E, st_rr_q975 = q_975/E) %>%
    mutate(resid = (obs_y - fit_y)/fit_sd)
  
  return(obsFitDat)
}
################################
import_obsFit_wksToEpi_ctySt <- function(modCodeStr_cty, modCodeStr_st, offset_l, filepathList){
  print(match.call())
  # import fitted values for county and state models: weeks to epidemic onset

  # import county and state data for models with offset
  if (offset_l){
    ctyDat <- import_obsFit_wksToEpi(modCodeStr_cty, filepathList) %>%
      mutate(fit_rr_cty = fit_y/E) %>%
      rename(fit_y_cty = fit_y, obs_y_cty = obs_y) %>%
      mutate(obs_rr_cty = obs_y_cty/E) %>%
      mutate(cty_LB = q_025/E, cty_UB = q_975/E) %>%
      select(season, fips, fit_rr_cty, fit_y_cty, cty_LB, cty_UB, obs_rr_cty, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obsFit_wksToEpi_st(modCodeStr_st, filepathList) %>%
      mutate(fit_rr_st = fit_y/E) %>%
      rename(fit_y_st = fit_y, obs_y_st = obs_y) %>%
      mutate(obs_rr_st = obs_y_st/E) %>%
      mutate(st_LB = q_025/E, st_UB = q_975/E) %>%
      select(season, fips_st, fit_rr_st, fit_y_st, st_LB, st_UB, obs_rr_st, obs_y_st)

    fullFitDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(fit_diff_stCty = fit_rr_st-fit_rr_cty) %>%
      mutate(obs_diff_stCty = obs_rr_st-obs_rr_cty) %>%
      select(season, fips, fips_st, fit_rr_cty, fit_rr_st, fit_diff_stCty, cty_LB, cty_UB, st_LB, st_UB, obs_rr_cty, obs_rr_st, obs_diff_stCty)

  } else{ # data without offset adjustment
    ctyDat <- import_obsFit_wksToEpi(modCodeStr_cty, filepathList) %>%
      rename(fit_y_cty = fit_y, obs_y_cty = obs_y) %>%
      rename(cty_LB = q_025, cty_UB = q_975) %>%
      select(season, fips, fit_y_cty, cty_LB, cty_UB, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obsFit_wksToEpi_st(modCodeStr_st, filepathList) %>%
      rename(fit_y_st = fit_y, obs_y_st = obs_y) %>%
      rename(st_LB = q_025, st_UB = q_975) %>%
      select(season, fips_st, fit_y_st, st_LB, st_UB, obs_y_st)

    fullFitDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(fit_diff_stCty = fit_y_st-fit_y_cty) %>%
      mutate(obs_diff_stCty = obs_y_st-obs_y_cty) %>%
      select(season, fips, fips_st, fit_y_cty, fit_y_st, fit_diff_stCty, cty_LB, cty_UB, st_LB, st_UB, obs_y_cty, obs_y_st, obs_diff_stCty)

  }
  
  return(fullFitDat) 
}
################################
import_obs_wksToEpi <- function(filepathList){
  print(match.call())
  # import observed for weeks to epidemic onset at county level
   
  # import observed and expected wks to epi
  inDat <- cleanR_wksToEpi_cty(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips, obs_y, E)
  # add lat/lon coords
  coordDat <- read_csv(filepathList$path_latlon_cty, col_types = "_c__dd")
  obsDat <- left_join(inDat, coordDat, by = c("fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))

  return(obsDat)
}
################################
import_obs_wksToEpi_st <- function(filepathList){
  print(match.call())
  # import observed for weeks to epidemic onset at state level
   
  # import observed and expected wks to epi
  obsDat <- cleanR_wksToEpi_st(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips_st, obs_y, E) %>%
        filter(!(fips_st %in% c("02", "15")))
   
  return(obsDat)
}
################################
import_obs_wksToEpi_ctySt <- function(offset_l, filepathList){
  print(match.call())
  # import fitted values for county and state models: weeks to epidemic onset

  # import county and state data for models with offset
  if (offset_l){
    ctyDat <- import_obs_wksToEpi(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      mutate(obs_rr_cty = obs_y_cty/E) %>%
      select(season, fips, latitude, longitude, obs_rr_cty, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_wksToEpi_st(modCodeStr_st, filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      mutate(obs_rr_st = obs_y_st/E) %>%
      select(season, fips_st, obs_rr_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_rr_st-obs_rr_cty) %>%
      mutate(obs_ratio_stCty = obs_y_st/obs_y_cty) %>%
      select(season, fips, fips_st, latitude, longitude, obs_rr_cty, obs_rr_st, obs_diff_stCty, obs_ratio_stCty) 

  } else{ # data without offset adjustment
    ctyDat <- import_obs_wksToEpi(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      select(season, fips, latitude, longitude, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_wksToEpi_st(filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      select(season, fips_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_y_st-obs_y_cty) %>%
      mutate(obs_ratio_stCty = obs_y_st/obs_y_cty) %>%
      select(season, fips, fips_st, latitude, longitude,  obs_y_cty, obs_y_st, obs_diff_stCty, obs_ratio_stCty) 
  }
  
  return(fullObsDat) 
}
################################
import_obs_wksToEpi_reg <- function(filepathList){
    print(match.call())
    # import observed data for weeks to epidemic onset at region level (no region level models)
    # acts as a wrapper for cleanR_wksToEpi_reg

    # import observed and expected wks to epi
    obsDat <- cleanR_wksToEpi_reg(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, regionID, obs_y, E) 
    
    return(obsDat)
}
################################
import_obs_wksToEpi_ctyReg <- function(offset_l, filepathList){
    print(match.call())
    # import observed values for county and region models: weeks to epidemic onset

    # import county and region data with offset
    if (offset_l){
        ctyDat <- import_obs_wksToEpi(filepathList) %>%
            mutate(obs_rr_cty = obs_y/E) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_wksToEpi_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            mutate(obs_rr_reg = obs_y/E) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_rr_reg-obs_rr_cty) %>%
          mutate(obs_ratio_regCty = obs_rr_reg/obs_rr_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_rr_cty, obs_rr_reg, obs_diff_regCty, obs_ratio_regCty)
    
    } else{ # data without offset adjustment
        ctyDat <- import_obs_wksToEpi(filepathList) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_wksToEpi_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_y_reg-obs_y_cty) %>%
          mutate(obs_ratio_regCty = obs_y_reg/obs_y_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_y_cty, obs_y_reg, obs_diff_regCty, obs_ratio_regCty)
    }
    
    return(fullObsDat)
}
#### wks.to.peak ################################
import_obs_wksToPeak <- function(filepathList){
  print(match.call())
  # import observed data for weeks to epidemic peak at county level
  # acts as a wrapper for cleanR_wksToPeak_cty
  
  # import observed and expected wks to epi
  inDat <- cleanR_wksToPeak_cty(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips, obs_y, E)
  # add lat/lon coords
  coordDat <- read_csv(filepathList$path_latlon_cty, col_types = "_c__dd")
  obsDat <- left_join(inDat, coordDat, by = c("fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_wksToPeak_st <- function(filepathList){
  print(match.call())
  # import observed data for weeks to epidemic peak at state level
  # acts as a wrapper for cleanR_wksToPeak_st
  
  # import observed and expected wks to epi
  obsDat <- cleanR_wksToPeak_st(filepathList) %>%
    mutate(obs_y = y1, E = E) %>%
    select(season, fips_st, obs_y, E) %>%
    filter(!(fips_st %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_wksToPeak_ctySt <- function(offset_l, filepathList){
  print(match.call())
  # import fitted values for county and state models: weeks to epidemic peak

  # import county and state data for models with offset
  if (offset_l){
    ctyDat <- import_obs_wksToPeak(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      mutate(obs_rr_cty = obs_y_cty/E) %>%
      select(season, fips, latitude, longitude, obs_rr_cty, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_wksToPeak_st(modCodeStr_st, filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      mutate(obs_rr_st = obs_y_st/E) %>%
      select(season, fips_st, obs_rr_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_rr_st-obs_rr_cty) %>%
      mutate(obs_ratio_stCty = obs_rr_st/obs_rr_cty) %>%
      select(season, fips, fips_st, latitude, longitude, obs_rr_cty, obs_rr_st, obs_diff_stCty, obs_ratio_stCty)

  } else{ # data without offset adjustment
    ctyDat <- import_obs_wksToPeak(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      select(season, fips, latitude, longitude,  obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_wksToPeak_st(filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      select(season, fips_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_y_st-obs_y_cty) %>%
      mutate(obs_ratio_stCty = obs_y_st/obs_y_cty) %>%
      select(season, fips, fips_st, latitude, longitude, obs_y_cty, obs_y_st, obs_diff_stCty, obs_ratio_stCty)
  }
  
  return(fullObsDat) 
}
################################
import_obs_wksToPeak_reg <- function(filepathList){
    print(match.call())
    # import observed data for weeks to peak at region level (no region level models)
    # acts as a wrapper for cleanR_wksToPeak_reg

    # import observed and expected wks to epi
    obsDat <- cleanR_wksToPeak_reg(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, regionID, obs_y, E)
    
    return(obsDat)
}
################################
import_obs_wksToPeak_ctyReg <- function(offset_l, filepathList){
    print(match.call())
    # import observed values for county and region models: weeks to peak

    # import county and region data with offset
    if (offset_l){
        ctyDat <- import_obs_wksToPeak(filepathList) %>%
            mutate(obs_rr_cty = obs_y/E) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_wksToPeak_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            mutate(obs_rr_reg = obs_y/E) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_rr_reg-obs_rr_cty) %>%
          mutate(obs_ratio_regCty = obs_rr_reg/obs_rr_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_rr_cty, obs_rr_reg, obs_diff_regCty, obs_ratio_regCty)
    
    } else{ # data without offset adjustment
        ctyDat <- import_obs_wksToPeak(filepathList) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_wksToPeak_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_y_reg-obs_y_cty) %>%
          mutate(obs_ratio_regCty = obs_y_reg/obs_y_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_y_cty, obs_y_reg, obs_diff_regCty, obs_ratio_regCty)
    }
    
    return(fullObsDat)
}
#### iliEarly ################################
import_obs_iliEarly <- function(filepathList){
  print(match.call())
  # import observed data for ili in early flu season at county level
  # acts as a wrapper for cleanR_iliEarly_irDt_shift1_cty
  
  # import observed and expected ili in early flu season
  inDat <- cleanR_iliEarly_irDt_shift1_cty(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips, obs_y, E)
  # add lat/lon coords
  coordDat <- read_csv(filepathList$path_latlon_cty, col_types = "_c__dd")
  obsDat <- left_join(inDat, coordDat, by = c("fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_iliEarly_st <- function(filepathList){
  print(match.call())
  # import observed data for ili in early flu season at state level
  # acts as a wrapper for cleanR_iliEarly_irDt_shift1_st
  
  # import observed and expected ili peak
  obsDat <- cleanR_iliEarly_irDt_shift1_st(filepathList) %>%
    mutate(obs_y = y1, E = E) %>%
    select(season, fips_st, obs_y, E) %>%
    filter(!(fips_st %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_iliEarly_ctySt <- function(offset_l, filepathList){
  print(match.call())
  # import observed values for county and state models: ili in early flu season

  # import county and state data for models with offset
  if (offset_l){
    ctyDat <- import_obs_iliEarly(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      mutate(obs_rr_cty = obs_y_cty/E) %>%
      select(season, fips, latitude, longitude, obs_rr_cty, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_iliEarly_st(filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      mutate(obs_rr_st = obs_y_st/E) %>%
      select(season, fips_st, obs_rr_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_rr_st-obs_rr_cty) %>%
      mutate(obs_ratio_stCty = obs_rr_st/obs_rr_cty) %>%
      select(season, fips, fips_st, obs_rr_cty, obs_rr_st, obs_diff_stCty, obs_ratio_stCty)

  } else{ # data without offset adjustment
    ctyDat <- import_obs_iliEarly(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      select(season, fips, latitude, longitude, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_iliEarly_st(filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      select(season, fips_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_y_st-obs_y_cty) %>%
      mutate(obs_ratio_stCty = obs_y_st/obs_y_cty) %>%
      select(season, fips, fips_st, latitude, longitude, obs_y_cty, obs_y_st, obs_diff_stCty, obs_ratio_stCty)

  }
  
  return(fullObsDat) 
}
################################
import_obs_iliEarly_reg <- function(filepathList){
    print(match.call())
    # import observed data for ili in early season at region level (no region level models)
    # acts as a wrapper for cleanR_iliEarly_shift1_reg

    # import observed and expected ili in early season
    obsDat <- cleanR_iliEarly_irDt_shift1_reg(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, regionID, obs_y, E)
    
    return(obsDat)
}
################################
import_obs_iliEarly_ctyReg <- function(offset_l, filepathList){
    print(match.call())
    # import observed values for county and region models: ili in early flu season

    # import county and region data with offset
    if (offset_l){
        ctyDat <- import_obs_iliEarly(filepathList) %>%
            mutate(obs_rr_cty = obs_y/E) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_iliEarly_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            mutate(obs_rr_reg = obs_y/E) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_rr_reg-obs_rr_cty) %>%
          mutate(obs_ratio_regCty = obs_rr_reg/obs_rr_cty) %>%
          select(season, fips, fips_st, regionID, obs_rr_cty, obs_rr_reg, obs_diff_regCty, obs_ratio_regCty)
    
    } else{ # data without offset adjustment
        ctyDat <- import_obs_iliEarly(filepathList) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_iliEarly_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_y_reg-obs_y_cty) %>%
          mutate(obs_ratio_regCty = obs_y_reg/obs_y_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_y_cty, obs_y_reg, obs_diff_regCty, obs_ratio_regCty)

    }
    
    return(fullObsDat)
}
#### iliPeak ################################
import_obs_iliPeak <- function(filepathList){
  print(match.call())
  # import observed data for peak ili at county level
  # acts as a wrapper for cleanR_iliPeak_irDt_shift1_cty
  
  # import observed and expected ili in early flu season
  inDat <- cleanR_iliPeak_irDt_shift1_cty(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips, obs_y, E)
  # add lat/lon coords
  coordDat <- read_csv(filepathList$path_latlon_cty, col_types = "_c__dd")
  obsDat <- left_join(inDat, coordDat, by = c("fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_iliPeak_st <- function(filepathList){
  print(match.call())
  # import observed data for ili peak at state level
  # acts as a wrapper for cleanR_iliPeak_irDt_shift1_st
  
  # import observed and expected ili peak
  obsDat <- cleanR_iliPeak_irDt_shift1_st(filepathList) %>%
    mutate(obs_y = y1, E = E) %>%
    select(season, fips_st, obs_y, E) %>%
    filter(!(fips_st %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_iliPeak_ctySt <- function(offset_l, filepathList){
  print(match.call())
  # import observed values for county and state models: peak ili

  # import county and state data for models with offset
  if (offset_l){
    ctyDat <- import_obs_iliPeak(filepathList) %>%
      rename(obs_y_cty = obs_y) %>%
      mutate(obs_rr_cty = obs_y_cty/E) %>%
      select(season, fips, latitude, longitude, obs_rr_cty, obs_y_cty) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_iliPeak_st(filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      mutate(obs_rr_st = obs_y_st/E) %>%
      select(season, fips_st, obs_rr_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_rr_st-obs_rr_cty) %>%
      mutate(obs_ratio_stCty = obs_rr_st/obs_rr_cty) %>%
      select(season, fips, fips_st, latitude, longitude, obs_rr_cty, obs_rr_st, obs_diff_stCty, obs_ratio_stCty)

  } else{ # data without offset adjustment
    ctyDat <- import_obs_iliPeak(filepathList) %>%
      select(season, fips, latitude, longitude,obs_y) %>%
      rename(obs_y_cty = obs_y) %>%
      mutate(fips_st = substring(fips, 1, 2))
  
    stDat <- import_obs_iliPeak_st(filepathList) %>%
      rename(obs_y_st = obs_y) %>%
      select(season, fips_st, obs_y_st)

    fullObsDat <- full_join(ctyDat, stDat, by = c("season", "fips_st")) %>%
      mutate(obs_diff_stCty = obs_y_st-obs_y_cty) %>%
      mutate(obs_ratio_stCty = obs_y_st/obs_y_cty) %>%
      select(season, fips, fips_st, latitude, longitude, obs_y_cty, obs_y_st, obs_diff_stCty, obs_ratio_stCty)

  }
  
  return(fullObsDat) 
}
################################
import_obs_iliPeak_reg <- function(filepathList){
    print(match.call())
    # import observed data for peak ili at region level (no region level models)
    # acts as a wrapper for cleanR_iliPeak_shift1_reg

    # import observed and expected peak ili
    obsDat <- cleanR_iliPeak_irDt_shift1_reg(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, regionID, obs_y, E)
    
    return(obsDat)
}
################################
import_obs_iliPeak_ctyReg <- function(offset_l, filepathList){
    print(match.call())
    # import observed values for county and region models: peak ili

    # import county and region data with offset
    if (offset_l){
        ctyDat <- import_obs_iliPeak(filepathList) %>%
            mutate(obs_rr_cty = obs_y/E) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_iliPeak_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            mutate(obs_rr_reg = obs_y/E) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_rr_reg-obs_rr_cty) %>%
          mutate(obs_ratio_regCty = obs_rr_reg/obs_rr_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_rr_cty, obs_rr_reg, obs_diff_regCty, obs_ratio_regCty)
    
    } else{ # data without offset adjustment
        ctyDat <- import_obs_iliPeak(filepathList) %>%
            rename(obs_y_cty = obs_y) %>%
            mutate(fips_st = substring(fips, 1, 2))

        reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1) %>%
          filter(!(fips_st %in% c("02", "15")))

        regDat <- import_obs_iliPeak_reg(filepathList) %>%
            left_join(reg_cw, by = c("regionID")) %>%
            rename(obs_y_reg = obs_y)

        fullObsDat <- full_join(ctyDat, regDat, by = c("season", "fips_st")) %>%
          mutate(obs_diff_regCty = obs_y_reg-obs_y_cty) %>%
          mutate(obs_ratio_regCty = obs_y_reg/obs_y_cty) %>%
          select(season, fips, fips_st, regionID, latitude, longitude, obs_y_cty, obs_y_reg, obs_diff_regCty, obs_ratio_regCty)

    }
    
    return(fullObsDat)
}
#### iliSum ################################
import_obs_iliSum <- function(filepathList){
  print(match.call())
  # import observed data for ili for total flu season at county level
  # acts as a wrapper for cleanR_iliSum_irDt_shift1_cty
  
  # import observed and expected ili for total flu season
  inDat <- cleanR_iliSum_irDt_shift1_cty(filepathList) %>%
        mutate(obs_y = y1, E = E) %>%
        select(season, fips, obs_y, E)
  # add lat/lon coords
  coordDat <- read_csv(filepathList$path_latlon_cty, col_types = "_c__dd")
  obsDat <- left_join(inDat, coordDat, by = c("fips")) %>%
    filter(!(substring(fips, 1, 2) %in% c("02", "15")))
  
  return(obsDat)
}
################################
import_obs_iliSum_st <- function(filepathList){
  print(match.call())
  # import observed data for ili for total flu season at state level
  # acts as a wrapper for cleanR_iliSum_irDt_shift1_st
  
  # import observed and expected ili for total flu season
  obsDat <- cleanR_iliSum_irDt_shift1_st(filepathList) %>%
    mutate(obs_y = y1, E = E) %>%
    select(season, fips_st, obs_y, E) %>%
    filter(!(fips_st %in% c("02", "15")))
  
  return(obsDat)
}

################################

#### data processing ################################
merge_obs_ctyStReg <- function(ctyStDat, ctyRegDat){
  print(match.call())

  fullDat <- full_join(ctyStDat, ctyRegDat %>% select(-fips_st, -contains("_cty"), -latitude, -longitude), by = c("season", "fips"))
  return(fullDat)
}
################################

#### paths  ################################
################################
string_fit_fname <- function(modCodeStr){
  searchDir <-  paste0(dirname(sys.frame(1)$ofile), "/../R_export/inlaModelData_export/", modCodeStr, "/")
  return(grep("summaryStatsFitted_", list.files(path = searchDir, full.names = TRUE), value = TRUE))
}
################################
string_exportFig_aggBias_model_folder <- function(){
  return(paste0(dirname(sys.frame(1)$ofile), "/../graph_outputs/aggBias_model_explore/"))
}
################################
string_exportFig_aggBias_data_folder <- function(){
  return(paste0(dirname(sys.frame(1)$ofile), "/../graph_outputs/aggBias_data_explore/"))
}
################################
string_msResults_folder <- function(){
    return(paste0(dirname(sys.frame(1)$ofile), "/../R_export/msResults/"))
}
################################
string_exportDat_response_data_folder <- function(){
  return(paste0(dirname(sys.frame(1)$ofile), "/../R_export/response_data_explore/"))
}
################################
string_exportDat_aggBias_data_folder <- function(){
  return(paste0(dirname(sys.frame(1)$ofile), "/../R_export/aggBias_data_explore/"))
}
################################
string_exportFig_wksToEpiAndPeak_folder <- function(){
  return(paste0(dirname(sys.frame(1)$ofile), "/../graph_outputs/wksToEpiAndPeak_explore/"))
}

#### plotting dependencies ################################
################################
substr.Right <- function(x, numchar){
  return(substr(x, nchar(x)-(numchar-1), nchar(x)))
}

import_county_geomMap <- function(){
  print(match.call())
  
  countyMap <- ggplot2::map_data("county")
  data(county.fips)
  polynameSplit <- tstrsplit(county.fips$polyname, ",")
  county_geomMap <- tbl_df(county.fips) %>%
    mutate(fips = substr.Right(paste0("0", fips), 5)) %>% 
    mutate(region = polynameSplit[[1]]) %>%
    mutate(subregion = polynameSplit[[2]]) %>%
    full_join(countyMap, by = c("region", "subregion")) %>%
    filter(!is.na(polyname) & !is.na(long)) %>%
    rename(state = region, county = subregion) %>%
    rename(region = fips) %>%
    select(-polyname)
  
  return(county_geomMap)
}

import_state_geomMap <- function(){
  print(match.call())
  
  stateMap <- ggplot2::map_data("state")
  data(state.fips)
  polynameSplit <- tstrsplit(state.fips$polyname, ":")
  state_geomMap <- tbl_df(state.fips) %>%
    mutate(fips = substr.Right(paste0("0", fips), 2)) %>% 
    mutate(region = polynameSplit[[1]]) %>%
    mutate(subregion = polynameSplit[[2]]) %>%
    full_join(stateMap, by = c("region", "subregion")) %>%
    filter(!is.na(polyname) & !is.na(long)) %>%
    rename(state = region) %>%
    rename(region = fips) %>%
    select(-polyname)
  
  return(state_geomMap)
}
