## 8/1/2021
## Identify overlap in % of population and % of counties with overlap between county and state source locations by season

library(tidyverse)
library(tidycensus)
# source("programs/census_api_key.R")

ctyOrig <- read_csv("R_export/origin_locations/fluseason_source_counties.csv") %>%
    dplyr::rename(ctySrcFips = srcFips) %>%
    dplyr::select(season, ctySrcFips)
stOrig <- read_csv("R_export/origin_locations/fluseason_source_states.csv") %>%
    dplyr::rename(stSrcFips = srcFips) %>%
    dplyr::select(season, stSrcFips)

overlap_df <- dplyr::mutate(ctyOrig, st_ctySrcFips = str_sub(ctySrcFips, 1, 2)) %>%
    dplyr::full_join(stOrig, by = c("season")) %>%
    dplyr::mutate(overlap = ifelse(st_ctySrcFips == stSrcFips, TRUE, FALSE))

ctyPop <- read_csv("R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_cty.csv") %>%
    dplyr::distinct(fips, year, pop) %>%
    dplyr::rename(ctyPop = pop)
stPop <- read_csv("R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_st.csv") %>%
    dplyr::distinct(state, fips_st, year, pop) %>%
    dplyr::rename(stPop = pop)

# census_api_key(census_api_key)

# stPop <- get_estimates(geography = "state",
#                         product = "population",
#                         time_series = TRUE)


only_overlap <- dplyr::filter(overlap_df, overlap) %>%
    dplyr::mutate(year = 2000+season) %>%
    dplyr::left_join(ctyPop, by = c("year", "ctySrcFips"="fips")) %>%
    dplyr::left_join(stPop, by = c("year", "stSrcFips" = "fips_st")) %>%
    dplyr::mutate(prop_overlap = ctyPop/stPop) %>%
    dplyr::group_by(season, stSrcFips) %>%
    dplyr::summarise(prop_overlap = sum(prop_overlap), numCty_overlap = n(), ctyPop_overlap = sum(ctyPop), stPop = first(stPop))

## results stated under "Probable epidemic source locations" results section