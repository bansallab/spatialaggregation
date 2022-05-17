library(tidyverse)
library(usmap)

source("programs/source_import_modeldata.R")

ctymap <- import_county_geomMap()
stmap <- import_state_geomMap()

#############################
## COUNTY ORIGINS
uqfips <- us_map(regions = "counties") %>% distinct(fips)
origin_cty <- read_csv("R_export/origin_locations/fluseason_source_counties.csv") %>% ## top 20 counties
  dplyr::rename(fips = srcFips) %>%
  right_join(uqfips) %>%
  tidyr::complete(season, fips) %>%
  dplyr::filter(!is.na(season)) %>%
  dplyr::select(season, fips, corrCoef)

seasons <- unique(origin_cty$season)

for(s in seasons){
  pltDat <- origin_cty %>% dplyr::filter(season == s) %>%
    dplyr::mutate(is_origin = ifelse(!is.na(corrCoef), TRUE, FALSE))
  cols <- c('TRUE' = "#B22222", 'FALSE' = "white")
  mapdf_cty <- ggplot() +
    geom_map(data = pltDat, map = ctymap, aes(fill = is_origin, map_id = fips), colour = "white", size = .025) +
     geom_map(data = stmap, map = stmap, aes(x = long, y = lat, map_id = region), colour = "grey50", fill = "transparent", size = .1) +
    scale_fill_manual(name = "Top 20 origin", values = cols, na.value = "white", drop = TRUE) +
    expand_limits(x = ctymap$long, y = ctymap$lat) +
    guides(fill = "none") +
    theme_minimal() +
    theme(text = element_text(size = 9), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.key.size = unit(.35, "cm"), legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), legend.title = element_text(size = 8))

  fname <- paste0("graph_outputs/explore_origins/season_origins_cty_S", s, ".png")
  ggsave(fname, mapdf_cty, width = 4, height = 2.75, units = "in")
}


#############################
## STATE ORIGINS
uqstate <- us_map(regions = "states") %>% distinct(fips) %>% dplyr::rename(fips_st = fips)
origin_st <- read_csv("R_export/origin_locations/fluseason_source_states.csv") %>% ## top 2 states
  dplyr::rename(fips_st = srcFips) %>%
  right_join(uqstate) %>%
  tidyr::complete(season, fips_st) %>%
  dplyr::select(season, fips_st, corrCoef)

seasons <- unique(origin_st$season)[which(!is.na(unique(origin_st$season)))]

for(s in seasons){
  pltDat <- origin_st %>% dplyr::filter(season == s) %>%
    dplyr::mutate(is_origin = ifelse(!is.na(corrCoef), TRUE, FALSE))
  cols <- c('TRUE' = "#B22222", 'FALSE' = "white")
  mapdf_st <- ggplot() +
    geom_map(data = stmap, map = stmap, aes(x = long, y = lat, map_id = region), colour = "grey50", size = .025) +
    geom_map(data = pltDat, map = stmap, aes(fill = is_origin, map_id = fips_st), colour = "grey50", size = .1) +
    scale_fill_manual(name = "Top 2 origin", values = cols, na.value = "white", drop = TRUE) +
    expand_limits(x = ctymap$long, y = ctymap$lat) +
    guides(fill = "none") +
    theme_minimal() +
    theme(text = element_text(size = 9), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.key.size = unit(.35, "cm"), legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), legend.title = element_text(size = 8))

  fname <- paste0("graph_outputs/explore_origins/season_origins_st_S", s, ".png")
  ggsave(fname, mapdf_st, width = 4, height = 2.75, units = "in")
}


#############################
## Combined county state origins
for(s in seasons){
  stDat <- origin_st %>% dplyr::filter(season == s) %>%
    dplyr::mutate(is_origin = ifelse(!is.na(corrCoef), "st", "no"))
  ctyDat <- origin_cty %>% dplyr::filter(season == s) %>%
    dplyr::mutate(is_origin = ifelse(!is.na(corrCoef), "cty", "no"))
  cols <- c('cty' = "#B22222", 'st' = "#d89090", 'no' = "transparent")
  bcols <- c('cty' = "transparent", 'st' = "transparent", 'no' = "transparent")
  mapdf_comb <- ggplot() +
    geom_map(data = stmap, map = stmap, aes(x = long, y = lat, map_id = region), colour = "grey50", size = .025, fill = "white") +
    geom_map(data = stDat, map = stmap, aes(fill = is_origin, map_id = fips_st, colour = is_origin), size = 0.3) +
    geom_map(data = ctyDat, map = ctymap, aes(fill = is_origin, map_id = fips), size = 0.3, colour = "transparent") +
    geom_map(data = stmap, map = stmap, aes(x = long, y = lat, map_id = region), colour = "grey50", fill = "transparent", size=0.5) +
    scale_colour_manual(values = bcols, drop = TRUE, na.value = "transparent") +
    scale_fill_manual(values = cols, drop = TRUE, na.value = "transparent") +
    expand_limits(x = ctymap$long, y = ctymap$lat) +
    guides(fill = "none", colour = "none") +
    theme_minimal() +
    theme(text = element_text(size = 9), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom", legend.key.size = unit(.35, "cm"), legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt"), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), legend.title = element_text(size = 8))

  fname <- paste0("graph_outputs/explore_origins/season_origins_comb_S", s, ".png")
  ggsave(fname, mapdf_comb, width = 4, height = 2.75, units = "in")
}