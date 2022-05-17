## Name: Elizabeth Lee
## Date: 10/25/17
## Function: functions to analyze response data
## Filenames: 
## Data Source: 
## Notes: 
################################

require(tidyverse)
require(lazyeval)
require(ncf)
setwd(dirname(sys.frame(1)$ofile))
# source("source_import_modeldata.R")

#### STATISTICS ################################
correlogStat_obs_allSeasons <- function(prepDat, datFormats){
  print(match.call())

  # data formatting
  dataScale <- datFormats$dataScale
  datFormats$statVar <- ifelse(datFormats$offset_l, paste0("obs_rr_", dataScale), paste0("obs_y_", dataScale))
  incrementKm <- datFormats$incrementKm
  resamp <- datFormats$resamp

  # plot formatting
  measure <- datFormats$measure
  dataProcess <- datFormats$dataProcess
  w <- datFormats$w; h <- datFormats$h; dp <- 300
  exportPath <- datFormats$exportPath

  # clean data
  statDat <- prepDat %>%
    rename_("statVar" = datFormats$statVar) %>%
    mutate(season = paste0("S", season)) %>%
    select(season, fips, latitude, longitude, statVar) %>%
    spread(season, statVar)

  # columns represent data over time
  statMx <- as.matrix(statDat %>% select(num_range("S", 3:9)))

  # calculate correlogram
  correlogOut <- correlog(x = statDat$longitude, y = statDat$latitude, z = statMx, na.rm = TRUE, increment = incrementKm, resamp = resamp, latlon = TRUE)
  # seems like increment has the unit km: https://stat.ethz.ch/pipermail/r-sig-geo/2010-October/009506.html

  # plot and export correlogram
  exportFname <- paste0(exportPath, "/", dataProcess, "_", measure, "/correlog_obs_", dataProcess, "_", measure, "_", dataScale, "_resamp", resamp, ".png")

  png(exportFname, units = "in", width = w, height = h, res = dp)
  plot.correlogMod(correlogOut, datFormats)
  dev.off()

  # write correlogram data
  exportDatname <- paste0(string_exportDat_response_data_folder(), "correlog_obs_", dataProcess, "_", measure, "_", dataScale, "_resamp", resamp, ".csv")

  exportCorrelog <- data.frame(dataProcess = dataProcess, measure = measure, dataScale = dataScale, xIntercept = correlogOut$x.intercept, correlation = correlogOut$correlation, meanOfClass = correlogOut$mean.of.class, numPairs = correlogOut$n, pValue = correlogOut$p)
  write_csv(exportCorrelog, exportDatname)

  return(correlogOut)
}

################################
plot.correlogMod <- function (x, datFormats){
    obj <- x
    plot(obj$mean.of.class, obj$correlation, ylab = "correlation", 
        xlab = "distance (mean-of-class, km)")
    lines(obj$mean.of.class, obj$correlation)
    abline(h = 0, col = "red")
    if (!is.null(obj$p)) {
        points(obj$mean.of.class[obj$p < 0.025], obj$correlation[obj$p < 
            0.025], pch = 21, bg = "black")
    }
    title(paste("Correlogram:", datFormats$dataProcess, datFormats$measure))
}

#### FIGURES ################################
################################
choro_obs_db_oneSeason <- function(prepDat, pltFormats){
  print(match.call())
  # plot single season choropleth for disease burden measure

  # plot formatting
  w <- pltFormats$w; h <- pltFormats$h; dp <- 300
  if (is.null(pltFormats$legendStep)){
    legendStep <- 1
  } else{
    legendStep <- pltFormats$legendStep
  }
  offset_l <- pltFormats$offset_l
  measure <- pltFormats$measure
  dataProcess <- pltFormats$dataProcess
  dataScale <- pltFormats$dataScale
  exportPath <- datFormats$exportPath

  names(prepDat) <- gsub(paste0("y_", dataScale), "y", names(prepDat)) 

  # set breaks based on distribution of observed data
  if(!offset_l){
    breaks <- seq(floor(min(prepDat$obs_y, na.rm = TRUE)), ceiling(max(prepDat$obs_y, na.rm = TRUE)), by = legendStep)
    prepDat2 <- prepDat %>%
      mutate(Observed = cut(obs_y, breaks, right = TRUE, include.lowest = TRUE, ordered_result = TRUE)) 
  } else{
    breaks <- seq(floor(min(prepDat$obs_rr, na.rm = TRUE)), ceiling(max(prepDat$obs_rr, na.rm = TRUE)), by = legendStep)
    prepDat2 <- prepDat %>%
      mutate(Observed = cut(obs_rr, breaks, right = TRUE, include.lowest = TRUE, ordered_result = TRUE)) 
  }
  
  factorlvls <- levels(prepDat2$Observed)
  plotDat <- prepDat2 %>%
    select(season, fips, Observed) %>%
    gather(fig, bin, Observed) %>%
    filter(fig == "Observed") %>%
    mutate(fig = factor(fig, levels = c("Observed"))) %>%
    mutate(bin = factor(bin, levels = factorlvls, labels = factorlvls, ordered = TRUE)) 
  print(levels(plotDat$bin))
 
  seasLs <- plotDat %>% distinct(season) %>% unlist

  
  plotChoro <- function(x){
    exportFname <- paste0(exportPath, "/", dataProcess, "_", measure, "/choro_obs_", dataProcess, "_", measure, "_", dataScale, "_S", x, ".png")
    pltDat <- plotDat %>% filter(season == x)

    # import county mapping info
    ctyMap <- import_county_geomMap()
    
    # plot
    choro <- ggplot() +
      geom_map(data = ctyMap, map = ctyMap, aes(x = long, y = lat, map_id = region)) +
      geom_map(data = pltDat, map = ctyMap, aes(fill = bin, map_id = fips), color = "grey25", size = 0.025) +
      scale_fill_brewer(name = paste0("Observed ", measure), palette = "OrRd", na.value = "grey60", drop = FALSE) +
      expand_limits(x = ctyMap$long, y = ctyMap$lat) +
      theme_minimal() +
      theme(text = element_text(size = 10), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
    
    ggsave(exportFname, choro, height = h, width = w, dpi = dp)
  } 
  
  purrr::map(seasLs, plotChoro)
     
}
################################
choro_obs_db_avgSeason <- function(inDat, pltFormats){
  # plot choropleth across seasons for disease burden measure
  
  # plot formatting
  w <- pltFormats$w; h <- pltFormats$h; dp <- 300
  if (is.null(pltFormats$legendStep)){
    legendStep <- 1
  } else{
    legendStep <- pltFormats$legendStep
  }
  offset_l <- pltFormats$offset_l
  measure <- pltFormats$measure
  dataProcess <- pltFormats$dataProcess
  dataScale <- pltFormats$dataScale
  exportPath <- pltFormats$exportPath

  names(inDat) <- gsub(paste0("y_", dataScale), "y", names(inDat)) 

  prepDat <- inDat %>%
    group_by(fips) %>%
    summarise(obs_y = mean(obs_y, na.rm = TRUE), E = mean(E, na.rm = TRUE)) %>%
    mutate(obs_rr = obs_y/E) %>%
    ungroup

  if(!offset_l){
    breaks <- seq(floor(min(prepDat$obs_y, na.rm = TRUE)), ceiling(max(prepDat$obs_y, na.rm = TRUE)), by = legendStep)
    prepDat2 <- prepDat %>%
      mutate(Observed = cut(obs_y, breaks, right = TRUE, include.lowest = TRUE, ordered_result = TRUE)) 
  } else{
    breaks <- seq(floor(min(prepDat$obs_rr, na.rm = TRUE)), ceiling(max(prepDat$obs_rr, na.rm = TRUE)), by = legendStep)
    prepDat2 <- prepDat %>%
      mutate(Observed = cut(obs_rr, breaks, right = TRUE, include.lowest = TRUE, ordered_result = TRUE)) 
  }
  factorlvls <- levels(prepDat2$Observed)
  pltDat <- prepDat2 %>%
    select(fips, Observed) %>%
    gather(fig, bin, Observed) %>%
    filter(fig == "Observed") %>%
    mutate(fig = factor(fig, levels = c("Observed"))) %>%
    mutate(bin = factor(bin, levels = factorlvls, labels = factorlvls, ordered = TRUE)) 
  print(levels(pltDat$bin))
 
  exportFname <- paste0(exportPath, "/", dataProcess, "_", measure, "/choro_obs_", dataProcess, "_", measure, "_", dataScale, "_avg.png")

  # import county mapping info
  ctyMap <- import_county_geomMap()
  
  # plot
  choro <- ggplot() +
    geom_map(data = ctyMap, map = ctyMap, aes(x = long, y = lat, map_id = region)) +
    geom_map(data = pltDat, map = ctyMap, aes(fill = bin, map_id = fips), color = "grey25", size = 0.025) +
    scale_fill_brewer(name = paste0("Observed ", measure), palette = "OrRd", na.value = "grey60", drop = FALSE) +
    expand_limits(x = ctyMap$long, y = ctyMap$lat) +
    theme_minimal() +
    theme(text = element_text(size = 10), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), legend.position = "bottom")
  
  ggsave(exportFname, choro, height = h, width = w, dpi = dp)
     
}
################################
