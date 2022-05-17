
## Name: Elizabeth Lee
## Date: 7/31/17
## Function: Functions for identifying the most probable origin location for an epidemic
## Filenames: 
## Data Source: 
## Notes: main program needs to import source_clean_response_functions.. for cleanR_wksToEpi functions (both cty and state level)
## 
## useful commands:
## install.packages("pkg", dependencies=TRUE, lib="/usr/local/lib/R/site-library") # in sudo R
## update.packages(lib.loc = "/usr/local/lib/R/site-library")

#### functions for model data cleaning ################################
require(tidyverse)

################################

subset_earliest_onset_locations_decile <- function(wksToEpi_data){
	# identify top 10% of locations with earliest epidemic onset in each season
	print(match.call())

	# wksToEpi_data columns: fips, county, st, stateID, lat, lon, season, year, pop, y, has.epi, incl.analysis
	quants <- wksToEpi_data %>%
		group_by(season) %>%
		summarise(cutoff = quantile(y1, probs = c(.1), names = FALSE, na.rm = TRUE)) %>%
		select(season, cutoff)

	earlyLocData <- left_join(wksToEpi_data, quants, by = c("season")) %>%
		filter(y1 <= cutoff) %>%
		mutate(srcID = paste0("srcID", seq_along(y1))) %>%
		select(srcID, season, contains("fips"), lat, lon) %>%
		rename(srcLat = lat, srcLon = lon) 

	return(earlyLocData)
}
################################
distance_function <- function(srcPt1, srcPt2, pt1, pt2){
	return(sqrt(((srcPt1-pt1)^2) + ((srcPt2-pt2)^2)))
}
################################

calculate_distance_from_sources <- function(earlyLocData, wksToEpi_data){
	# calculate distance between potential source locations and all counties for data from a single season
	print(match.call())

	dummy_distPts <- list()

	# calculate distances between source locations (season-specific) and all fips in a given season
	for (i in 1:nrow(earlyLocData)){
		srcLatPts <- rep(earlyLocData[i,]$srcLat, nrow(wksToEpi_data))
		srcLonPts <- rep(earlyLocData[i,]$srcLon, nrow(wksToEpi_data))

		distPts <- distance_function(srcLatPts, srcLonPts, wksToEpi_data$lat, wksToEpi_data$lon) 
		dummy_distPts[[earlyLocData[i,]$srcID]] <- distPts
	}
	
	# merge distance information for each potential source location with onset week 
	distData <- bind_cols(wksToEpi_data %>% select(contains("fips"), season, y1), as.data.frame(dummy_distPts))

	return(distData)

}
################################

calculate_correlation_distance_onsetWeek <- function(distData){
	# 
	print(match.call())

	onsetVec <- distData$y1
	evalCols <- names(distData %>% select(matches("srcID")))
	corrData <- distData %>% 
		summarise_at(vars(evalCols), cor, y = onsetVec, method = "pearson", use = "complete.obs")

	return(corrData)
}