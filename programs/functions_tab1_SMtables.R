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
require(tidyverse)

#### INLA model and export functions #################################
run_intercept_model <- function(fullData, formula){
  print(match.call())

  modelOut <- inla(formula,
    family = "gaussian",
    data = list(Y = fullData$diff_aggBias,
                graphIdx = fullData$graphIdx),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE),
    verbose = TRUE,
    debug = TRUE)

  return(modelOut)
}
#################################
grab_dicCPO <- function(modelOut, formulaType){
  print(match.call())

  dic <- modelOut$dic$dic
  cpo <- sum(log(modelOut$cpo$cpo), na.rm = TRUE)
  cpoFail <- sum(modelOut$cpo$failure, na.rm = TRUE)

  dicCPO <- data.frame(formula = formulaType, exportDate = as.character(Sys.Date()), dic = dic, cpo = cpo, cpoFail = cpoFail)

  return(dicCPO)
}
#################################
grab_summaryStats <- function(modelOut, formulaType){
  print(match.call())
  
  columnNames <- c("mean", "sd", "q_025", "q_5", "q_975", "mode", "kld")
  summaryComplete <- data.frame()

  # clean fixed effects summary statistics output from INLA
  if (nrow(modelOut$summary.fixed)>0){
    names(modelOut$summary.fixed) <- columnNames
    summaryFixed <- tbl_df(modelOut$summary.fixed) %>%
      mutate(RV = rownames(modelOut$summary.fixed)) %>%
      mutate(effectType = "fixed") %>%
       select(RV, effectType, mean, sd, q_025, q_5, q_975, mode, kld)
    summaryComplete <- bind_rows(summaryComplete, summaryFixed)
  }
  # clean hyperpar summary statistics output from INLA
  if (!is.null(modelOut$summary.hyperpar)){
    names(modelOut$summary.hyperpar) <- columnNames[1:6] 
    summaryHyperpar <- tbl_df(modelOut$summary.hyperpar) %>%
      mutate(RV = rownames(modelOut$summary.hyperpar)) %>%
      mutate(effectType = "hyperpar", kld = NA) %>%
      select(RV, effectType, mean, sd, q_025, q_5, q_975, mode, kld)
    summaryComplete <- bind_rows(summaryComplete, summaryHyperpar)
  }
  if (!is.null(modelOut$summary.random$graphIdx)){
    names(modelOut$summary.random$graphIdx) <- c("RV", columnNames)
    summaryRandomGraphid <- modelOut$summary.random$graphIdx %>% 
      mutate(RV = as.character(paste0("error", RV))) %>%
      mutate(effectType = "random") %>%
      select(RV, effectType, mean, sd, q_025, q_5, q_975, mode, kld)
    summaryComplete <- bind_rows(summaryComplete, summaryRandomGraphid)
  }
  
  # bind data together
  summaryStats <- summaryComplete %>%
    mutate(formula = formulaType) %>%
    mutate(exportDate = as.character(Sys.Date())) %>%
    select(formula, exportDate, RV, effectType, mean, sd, q_025, q_5, q_975, mode, kld)
  
  return(summaryStats)  
}
#################################

#### data cleaning functions #################################
create_response_data <- function(data1, data2, dataFormats, path_list){
  print(match.call())
  # create new response variable y, which represents the difference in aggregation bias between the two datasets

  variable1 <- dataFormats$variable1; variable2 <- dataFormats$variable2

  newdata1 <- data1 %>%
    rename_("var1" = variable1) %>%
    select(-fips_st, -latitude, -longitude, -contains("obs_y"))
  newdata2 <- data2 %>%
    rename_("var2" = variable2) %>%
    select(-contains("obs_y"))

  graphIdx_df <- clean_graphIDx(path_list, "county")

  fullData <- full_join(newdata1, newdata2, by = c("season", "fips")) %>%
    mutate(diff_aggBias = var1-var2) %>%
    left_join(graphIdx_df, by = c("fips"))

  return(fullData)
}
#################################
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
#################################