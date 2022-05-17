
## Name: Elizabeth Lee
## Date: 8/11/17
## Function: State-specific measures to see which states are more heterogeneous than others (variance in timing); Later, can use this to compare the variance within states (range of values or variance rank) to the aggregation effect

## Filenames: physicianCoverage_IMSHealth_state.csv, dbMetrics_periodicReg_ilinDt_Octfit_span0.4_degree2_analyzeDB_st.csv
## Data Source: IMS Health
## Notes: 
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### header #################################
rm(list = ls())
require(tidyverse); require(DBI); require(RMySQL) # clean_data_functions dependencies

#### set these! ################################
dbCodeStr <- "_irDt_Octfit_span0.4_degree2"

#### SOURCE: clean and import model data #################################
setwd(dirname(sys.frame(1)$ofile))
source("source_clean_response_functions_cty.R") # functions to clean response and IMS coverage data (cty)

#### FILEPATHS #################################
setwd('../reference_data')
path_abbr_st <- paste0(getwd(), "/state_abbreviations_FIPS.csv")
path_latlon_cty <- paste0(getwd(), "/cty_pop_latlon.csv")

setwd("../R_export")
path_response_cty <- paste0(getwd(), sprintf("/dbMetrics_periodicReg%s_analyzeDB_cty.csv", dbCodeStr))

path_list <- list(path_abbr_st = path_abbr_st,
                  path_latlon_cty = path_latlon_cty,
                  path_response_cty = path_response_cty)

setwd("../graph_outputs/explore_withinStateHeterogeneity/")
path_exportFig <- getwd()

setwd("../../R_export/test_dbVariance_aggBiasMagnitude/")
path_exportData <- getwd()

#### Functions #################################
scatter_variance <- function(prepData, plotFormats){
  # filter the data for a single season if needed
  print(match.call())
  
  ylab <- plotFormats$ylabScatter; dbCode <- plotFormats$dbCode

  pltData <- prepData %>%
    mutate(xAxis = factor(hetRank, labels = st))
  
  plt <- ggplot(pltData, aes(x = xAxis, y = variance)) +
    geom_point() +
    scale_y_continuous(ylab) +
    theme_bw() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1), axis.text=element_text(size=10), text = element_text(size = 10))
  
  ggsave(paste0(path_exportFig, "/", dbCode, "/scatter_variance_", dbCode, ".jpeg"), height = 4, width = 6, units = "in")
}
#################################
boxplot_response <- function(respData, plotFormats){
  # for all seasons
  print(match.call())

  ylab <- plotFormats$ylabBoxplot; dbCode <- plotFormats$dbCode

  pltData <- respData #%>%
    # mutate(xAxis = factor(hetRank, labels = st))
  
  plt <- ggplot(pltData, aes(st, y = y1)) +
    geom_boxplot() +
    # scale_x_continuous("", labels = )
    scale_y_continuous(ylab) +
    theme_bw() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1), axis.text=element_text(size=10), text = element_text(size = 10)) + 
    facet_grid(season~.)
  
  ggsave(paste0(path_exportFig, "/", dbCode, "/boxplot_response_", dbCode, ".jpeg"), height = 8, width = 8, units = "in")
}
#################################
boxplot_response_oneSeas <- function(respData, plotFormats){
  # for all seasons
  print(match.call())

  ylab <- plotFormats$ylabBoxplot; dbCode <- plotFormats$dbCode

  seasons <- respData %>% distinct(season) %>% unlist
  for (s in seasons){
    seasData <- respData %>% filter(season == s & !is.na(hetRank))
    labelsDf <- seasData %>%
      distinct(hetRank, st) %>% 
      arrange(hetRank)
    
    pltData <- seasData %>%
        mutate(xAxis = factor(hetRank, levels = labelsDf$hetRank, labels = labelsDf$st))
  
    plt <- ggplot(pltData, aes(xAxis, y = y1)) +
        geom_boxplot() +
        ggtitle(paste("Season", s)) +
        scale_y_continuous(ylab) +
        theme_bw() + 
        theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1), axis.text=element_text(size=10), text = element_text(size = 10))

    ggsave(paste0(path_exportFig, "/", dbCode, "/boxplot_response_", dbCode, "_S", s, ".jpeg"), height = 4, width = 6, units = "in")
  }
  
}

#### Process data #################################
# pool data across seasons at state level
process_pooledData <- function(yData){
  return(yData %>%
    mutate(fips_st = substr.Right(paste0("0", stateID), 2)) %>%
    group_by(fips_st) %>%
    mutate(forCount = ifelse(!is.na(y1), 1, 0)) %>%
    summarise(st = first(st), variance = var(y1, na.rm = TRUE), counted = sum(forCount), totCty = length(y1)) %>%
    filter(st != "DC" & counted > 0) %>%
    mutate(hetRank = data.table::frank(-variance)) %>%
    ungroup
  )
}
#################################
# pool data across seasons at county level
process_pooledData_cty <- function(yData){
  return(yData %>%
    group_by(fips) %>%
    mutate(forCount = ifelse(!is.na(y1), 1, 0)) %>%
    summarise(variance = var(y1, na.rm = TRUE), counted = sum(forCount), totCty = length(y1)) %>%
    mutate(hetRank = data.table::frank(-variance)) %>%
    ungroup
  )
}
#################################
# group data by season and link ranks and full response data
process_respData <- function(yData){
  varDat_bySeas <- yData %>%
    mutate(forCount = ifelse(!is.na(y1), 1, 0)) %>%
    mutate(fips_st = substr.Right(paste0("0", stateID), 2)) %>%
    group_by(fips_st, season) %>%
    summarise(st = first(st), variance = var(y1, na.rm = TRUE), counted = sum(forCount), totCty = length(y1)) %>%
    filter(st != "DC" & counted > 0) %>%
    arrange(season, desc(variance)) %>%
    group_by(season) %>%
    mutate(hetRank = data.table::frank(-variance)) %>%
    ungroup

  return(left_join(yData, varDat_bySeas %>% select(st, season, hetRank), by = c("season", "st")))
}
#################################
write_pooledData <- function(processedData, exportPath){
  writeData <- processedData %>%
    rename(contribCty = counted) 
  write_csv(writeData, exportPath) 
}
#################################


#### MAIN #################################
# import data
wksToEpi <- cleanR_wksToEpi_cty(path_list)
wksToEpi_pooled <- process_pooledData(wksToEpi)
wksToEpi_pooled_cty <- process_pooledData_cty(wksToEpi)
pltFormats_wksToEpi <- list(dbCode = "irDt_wksToEpi", ylabScatter = "Variance in epidemic onset", ylabBoxplot = "Weeks to epidemic onset")

wksToPeak <- cleanR_wksToPeak_cty(path_list)
wksToPeak_pooled <- process_pooledData(wksToPeak)
wksToPeak_pooled_cty <- process_pooledData_cty(wksToPeak)
pltFormats_wksToPeak <- list(dbCode = "irDt_wksToPeak", ylabScatter = "Variance in peak timing", ylabBoxplot = "Weeks to peak")

iliEarly <- cleanR_iliEarly_irDt_shift1_cty(path_list)
iliEarly_pooled <- process_pooledData(iliEarly)
iliEarly_pooled_cty <- process_pooledData_cty(iliEarly)
pltFormats_iliEarly <- list(dbCode = "irDt_iliEarly", ylabScatter = "Variance in early seasonal intensity", ylabBoxplot = "Early seasonal intensity")

iliPeak <- cleanR_iliPeak_irDt_shift1_cty(path_list)
iliPeak_pooled <- process_pooledData(iliPeak)
iliPeak_pooled_cty <- process_pooledData_cty(iliPeak)
pltFormats_iliPeak <- list(dbCode = "irDt_iliPeak", ylabScatter = "Variance in peak seasonal intensity", ylabBoxplot = "Peak seasonal intensity")

#### write pooled data ####
# state level
write_pooledData(wksToEpi_pooled, paste0(path_exportData, "/dbVariance_st_irDt_wksToEpi.csv"))
write_pooledData(wksToPeak_pooled, paste0(path_exportData, "/dbVariance_st_irDt_wksToPeak.csv"))
write_pooledData(iliEarly_pooled, paste0(path_exportData, "/dbVariance_st_irDt_iliEarly.csv"))
write_pooledData(iliPeak_pooled, paste0(path_exportData, "/dbVariance_st_irDt_iliPeak.csv"))
# county level
write_pooledData(wksToEpi_pooled_cty, paste0(path_exportData, "/dbVariance_cty_irDt_wksToEpi.csv"))
write_pooledData(wksToPeak_pooled_cty, paste0(path_exportData, "/dbVariance_cty_irDt_wksToPeak.csv"))
write_pooledData(iliEarly_pooled_cty, paste0(path_exportData, "/dbVariance_cty_irDt_iliEarly.csv"))
write_pooledData(iliPeak_pooled_cty, paste0(path_exportData, "/dbVariance_cty_irDt_iliPeak.csv"))

#### plot variance across states ####
# wks to epi
scatter_variance(wksToEpi_pooled, pltFormats_wksToEpi)
boxplot_response(process_respData(wksToEpi), pltFormats_wksToEpi)
boxplot_response_oneSeas(process_respData(wksToEpi), pltFormats_wksToEpi)
# wks to peak
scatter_variance(wksToPeak_pooled, pltFormats_wksToPeak)
boxplot_response(process_respData(wksToPeak), pltFormats_wksToPeak)
boxplot_response_oneSeas(process_respData(wksToPeak), pltFormats_wksToPeak)
# iliEarly
scatter_variance(iliEarly_pooled, pltFormats_iliEarly)
boxplot_response(process_respData(iliEarly), pltFormats_iliEarly)
boxplot_response_oneSeas(process_respData(iliEarly), pltFormats_iliEarly)
# iliPeak
scatter_variance(iliPeak_pooled, pltFormats_iliPeak)
boxplot_response(process_respData(iliPeak), pltFormats_iliPeak)
boxplot_response_oneSeas(process_respData(iliPeak), pltFormats_iliPeak)
