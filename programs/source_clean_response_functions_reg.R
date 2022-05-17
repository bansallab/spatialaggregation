## Name: Elizabeth Lee
## Date: 9/15/17
## Function: Functions for cleaning disease burden response data at the region level 
## Filenames: 
## Data Source: 
## Notes: 
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### functions for model data cleaning ################################
require(tidyverse)
source("programs/source_clean_response_functions_st.R")

##### REGION-LEVEL VARIABLES ##########################################

cleanR_iliEarly_shift1_reg <- function(filepathList){
    # clean response variable: ilinDt.early plus 1 (so it is comparable with iliSum)
    print(match.call())

    # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
    dbCode <- grep("ili+", strsplit(filepathList$path_response_reg, "_")[[1]], value=T)
    # clean data
    iliEarly_data <- read_csv(filepathList$path_response_reg, col_types = "icllcd") %>%
        filter(metric == sprintf("%s.early", dbCode)) %>%
        select(-metric) %>%
        mutate(regionID = as.numeric(substring(region, 2, nchar(region)))) %>%
        rename(y = burden)

    print(filepathList$path_response_reg)
    print(summary(iliEarly_data))
    pop_data <- clean_pop_reg(filepathList)

    return_data <- full_join(iliEarly_data, pop_data, by = c("season", "regionID")) %>% 
        mutate(year = 2000 + season) %>%
        select(regionID, lat ,lon, season, year, pop, y, has.epi) %>%
        mutate(y1 = y+1) %>%
        group_by(season) %>%
        mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
        ungroup %>%
        filter(season >= 3 & season <= 9)

    return(return_data)
}
################################

cleanR_iliPeak_shift1_reg <- function(filepathList){
    # clean response variable: ilinDt.peak plus 1 (so it is comparable with iliSum)
    print(match.call())

    # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
    dbCode <- grep("ili+", strsplit(filepathList$path_response_reg, "_")[[1]], value=T)
    # clean data
    iliPeak_data <- read_csv(filepathList$path_response_reg, col_types = "icllcd") %>%
        filter(metric == sprintf("%s.peak", dbCode)) %>%
        select(-metric) %>%
        mutate(regionID = as.numeric(substring(region, 2, nchar(region)))) %>%
        rename(y = burden)

    print(filepathList$path_response_reg)
    print(summary(iliPeak_data))
    pop_data <- clean_pop_reg(filepathList)

    return_data <- full_join(iliPeak_data, pop_data, by = c("season", "regionID")) %>% 
        mutate(year = 2000 + season) %>%
        select(regionID, lat ,lon, season, year, pop, y, has.epi) %>%
        mutate(y1 = y+1) %>%
        group_by(season) %>%
        mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
        ungroup %>%
        filter(season >= 3 & season <= 9)

    return(return_data)
}
################################

cleanR_iliEarly_irDt_shift1_reg <- function(filepathList){
    # clean response variable: irDt.early plus 1 (so it is comparable with iliSum)
    print(match.call())

    # grab disease burden metric (e.g., irDt): match "ili" 1+ times
    dbCode <- grep("irDt+", strsplit(filepathList$path_response_reg, "_")[[1]], value=T)
    # clean data
    iliEarly_data <- read_csv(filepathList$path_response_reg, col_types = "icllcd") %>%
        filter(metric == sprintf("%s.early", dbCode)) %>%
        select(-metric) %>%
        mutate(regionID = as.numeric(substring(region, 2, nchar(region)))) %>%
        rename(y = burden)

    print(filepathList$path_response_reg)
    print(summary(iliEarly_data))
    pop_data <- clean_pop_reg(filepathList)

    return_data <- full_join(iliEarly_data, pop_data, by = c("season", "regionID")) %>% 
        mutate(year = 2000 + season) %>%
        select(regionID, lat ,lon, season, year, pop, y, has.epi) %>%
        mutate(y1 = y+1) %>%
        group_by(season) %>%
        mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
        ungroup %>%
        filter(season >= 3 & season <= 9)

    return(return_data)
}
################################

cleanR_iliPeak_irDt_shift1_reg <- function(filepathList){
    # clean response variable: irDt.peak plus 1 (so it is comparable with iliSum)
    print(match.call())

    # grab disease burden metric (e.g., ilinDt): match "ili" 1+ times
    dbCode <- grep("irDt+", strsplit(filepathList$path_response_reg, "_")[[1]], value=T)
    # clean data
    iliPeak_data <- read_csv(filepathList$path_response_reg, col_types = "icllcd") %>%
        filter(metric == sprintf("%s.peak", dbCode)) %>%
        select(-metric) %>%
        mutate(regionID = as.numeric(substring(region, 2, nchar(region)))) %>%
        rename(y = burden)

    print(filepathList$path_response_reg)
    print(summary(iliPeak_data))
    pop_data <- clean_pop_reg(filepathList)

    return_data <- full_join(iliPeak_data, pop_data, by = c("season", "regionID")) %>% 
        mutate(year = 2000 + season) %>%
        select(regionID, lat ,lon, season, year, pop, y, has.epi) %>%
        mutate(y1 = y+1) %>%
        group_by(season) %>%
        mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
        ungroup %>%
        filter(season >= 3 & season <= 9)

    return(return_data)
}
################################

cleanR_wksToEpi_reg <- function(filepathList){
    # clean response variable: wks.to.epi from beginning of flu period
    print(match.call())

    # clean data
    wksToEpi_data <- read_csv(filepathList$path_response_reg, col_types = "icllcd") %>%
        filter(metric == "wks.to.epi") %>%
        select(-metric) %>%
        mutate(regionID = as.numeric(substring(region, 2, nchar(region)))) %>%
        rename(y = burden)

    print(filepathList$path_response_reg)
    pop_data <- clean_pop_reg(filepathList)

    return_data <- full_join(wksToEpi_data, pop_data, by = c("season", "regionID")) %>% 
        mutate(year = 2000 + season) %>%
        select(regionID, lat ,lon, season, year, pop, y, has.epi) %>%
        mutate(y1 = ifelse(y>0, y, NA)) %>%
        group_by(season) %>%
        mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
        ungroup %>%
        filter(season >= 3 & season <= 9)

    return(return_data)
}
################################

cleanR_wksToPeak_reg <- function(filepathList){
    # clean response variable: wks.to.peak from beginning of flu period
    print(match.call())

    # clean data
    wksToPeak_data <- read_csv(filepathList$path_response_reg, col_types = "icllcd") %>%
        filter(metric == "wks.to.peak") %>%
        select(-metric) %>%
        mutate(regionID = as.numeric(substring(region, 2, nchar(region)))) %>%
        rename(y = burden)

    print(filepathList$path_response_reg)
    pop_data <- clean_pop_reg(filepathList)

    return_data <- full_join(wksToPeak_data, pop_data, by = c("season", "regionID")) %>% 
        mutate(year = 2000 + season) %>%
        select(regionID, lat ,lon, season, year, pop, y, has.epi) %>%
        mutate(y1 = ifelse(y>0, y, NA)) %>%
        group_by(season) %>%
        mutate(E = weighted.mean(y1, pop, na.rm=TRUE)) %>%
        ungroup %>%
        filter(season >= 3 & season <= 9)

    return(return_data)
}
################################

##########################################

clean_pop_reg <- function(filepathList){
    # clean pop data at region level
    print(match.call())

    st_pop_data <- clean_pop_st(filepathList) %>%
        select(fips_st, season, year, pop)
    reg_cw <- read_csv(filepathList$path_region_cw, col_types = "c__i", col_names = c("fips_st", "regionID"), skip = 1)
    coord_data <- read_csv(filepathList$path_latlon_reg, col_types = "cdd", col_names = c("regionID", "lat", "lon"), skip = 1) %>%
        mutate(regionID = as.numeric(substring(regionID, 2, nchar(regionID))))

    dummy <- full_join(st_pop_data, reg_cw, by = c("fips_st")) %>%
        group_by(regionID, season) %>%
        summarise(pop = sum(pop)) %>%
        filter(season > 1)

    fulldata <- left_join(dummy, coord_data, by = "regionID")

    return(fulldata)
}