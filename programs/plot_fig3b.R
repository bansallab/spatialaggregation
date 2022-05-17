library(tidyverse)
library(ggthemes)

ctyDat <- read_csv("R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_cty.csv")
stDat <- read_csv("R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_st.csv")
regDat <- read_csv("R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_reg.csv")
refs <- read_csv("reference_data/state_abbreviations_FIPS_region.csv") %>%
  dplyr::mutate(region = paste0("R", Region)) %>%
  dplyr::rename(fips_st = fips, abbr_st = Abbreviation) %>%
  dplyr::select(fips_st, region, abbr_st)
refs_reg <- refs %>% dplyr::select(fips_st, region)
refs_st <- refs %>% dplyr::select(fips_st, abbr_st)

stpops <- stDat %>%
  dplyr::filter(lubridate::month(Thu.week) < 9) %>%
  distinct(fips_st, season, pop) %>%
  dplyr::rename(pop_st = pop)

####################################
## Process peak data ##

## process cty data and dummy fill in intermediate weeks where there was no change in cumpop_prop
weeksTemplate <- ctyDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & flu.week) %>%
  dplyr::mutate(fips_st = str_sub(fips, 1, 2)) %>%
  dplyr::distinct(fips_st, season, Thu.week) %>%
  dplyr::rename(peakwk_cty = Thu.week) %>%
  dplyr::mutate(uqid = paste(season, fips_st, peakwk_cty, by = "_"))
cty_peak_1 <- ctyDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  group_by(fips, season) %>%
  dplyr::filter(ir.dt == max(ir.dt)) %>%
  ungroup %>%
  dplyr::mutate(fips_st = str_sub(fips, 1, 2)) %>%
  left_join(refs_reg, by = c("fips_st")) %>%
  left_join(stpops, by = c("fips_st", "season")) %>%
  group_by(Thu.week, fips_st) %>%
  summarise(region = first(region), season = first(season), pop_cty = sum(pop), num_cty = n(), pop_st = first(pop_st)) %>%
  ungroup %>%
  dplyr::rename(peakwk_cty = Thu.week) %>%
  dplyr::mutate(uqid = paste(season, fips_st, peakwk_cty, by = "_"))
dummyWeekTemplate <- weeksTemplate %>%
  dplyr::filter(!(uqid %in% cty_peak_1$uqid)) %>%
  dplyr::select(-uqid) %>%
  left_join(refs_reg, by = c("fips_st")) %>%
  left_join(stpops, by = c("fips_st", "season")) %>%
  dplyr::mutate(pop_cty = 0)
cty_peak <- cty_peak_1 %>%
  dplyr::select(-uqid) %>%
  bind_rows(dummyWeekTemplate) %>%
  group_by(season, fips_st) %>%
  dplyr::arrange(fips_st, peakwk_cty) %>% 
  dplyr::mutate(cumpop_prop_cty = cumsum(pop_cty)/pop_st) %>%
  ungroup %>%
  dplyr::mutate(cumpop_prop_cty = ifelse(cumpop_prop_cty>1, 1, cumpop_prop_cty)) ## a few were slightly over 1 just due to pop rounding errors

## dummy fill in cumpop_prop_cty as 0 for weeks before seasonp peak
inRealData <- cty_peak %>%
  distinct(peakwk_cty, fips_st) %>%
  dplyr::mutate(uqid = paste(peakwk_cty, fips_st, by = "_"))
realPeak <- cty_peak %>%
  group_by(fips_st, season) %>%
  dplyr::filter(peakwk_cty == min(peakwk_cty)) %>%
  ungroup %>%
  dplyr::select(fips_st, season, peakwk_cty) %>%
  dplyr::rename(realPeakWk = peakwk_cty)
dummyStartEnd_st <- stDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & flu.week) %>%
  group_by(fips_st, season) %>%
  dplyr::mutate(cumpop_prop_cty = 0) %>%
  ungroup %>%
  dplyr::rename(peakwk_cty = Thu.week) %>%
  dplyr::select(fips_st, season, peakwk_cty, cumpop_prop_cty) %>%
  left_join(refs_reg, by = c("fips_st")) %>%
  left_join(realPeak, by = c("fips_st", "season")) %>%
  dplyr::filter(peakwk_cty < realPeakWk) %>%
  dplyr::mutate(uqid = paste(peakwk_cty, fips_st, by = "_")) %>%
  dplyr::filter(!(uqid %in% inRealData$uqid)) %>%
  dplyr::select(-uqid, -realPeakWk)

st_peak <- stDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  group_by(fips_st, season) %>%
  dplyr::filter(ir.dt == max(ir.dt)) %>%
  ungroup %>%
  dplyr::rename(peakwk_st = Thu.week) %>%
  dplyr::select(fips_st, season, peakwk_st) 

reg_peak <- regDat %>% 
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  group_by(region, season) %>%
  dplyr::filter(ir.dt == max(ir.dt)) %>%
  ungroup %>%
  dplyr::rename(peakwk_reg = Thu.week) %>%
  dplyr::select(region, season, peakwk_reg)

state_colors <- refs %>%
  group_by(region) %>%
  arrange(fips_st) %>%
  dplyr::mutate(plt.st.color = factor(seq_along(fips_st))) %>%
  ungroup %>%
  dplyr::select(fips_st, plt.st.color)


peakDat <- bind_rows(cty_peak, dummyStartEnd_st) %>%
  dplyr::arrange(fips_st, peakwk_cty) %>%
  full_join(st_peak, by = c("fips_st", "season")) %>%
  full_join(reg_peak, by = c("region", "season")) %>%
  left_join(state_colors, by = c("fips_st")) %>%
  dplyr::mutate(plt.date_cty = ifelse(lubridate::month(peakwk_cty)>=11, paste0("2005", substring(peakwk_cty, 5, 10)), paste0("2006", substring(peakwk_cty, 5, 10)))) %>%
  dplyr::mutate(plt.date_cty = lubridate::as_date(plt.date_cty)) %>%
  dplyr::mutate(plt.date_st = ifelse(lubridate::month(peakwk_st)>=11, paste0("2005", substring(peakwk_st, 5, 10)), paste0("2006", substring(peakwk_st, 5, 10)))) %>%
  dplyr::mutate(plt.date_st = lubridate::as_date(plt.date_st)) %>%
  dplyr::mutate(plt.date_reg = ifelse(lubridate::month(peakwk_reg)>=11, paste0("2005", substring(peakwk_reg, 5, 10)), paste0("2006", substring(peakwk_reg, 5, 10)))) %>%
  dplyr::mutate(plt.date_reg = lubridate::as_date(plt.date_reg)) %>%
  dplyr::mutate(plt.season = factor(season, levels = 3:9, labels = c("2002-2003", "2003-2004", "2004-2005", "2005-2006", "2006-2007", "2007-2008", "2008-2009"))) %>%
  dplyr::mutate(plt.region = factor(region, levels = paste0("R", 1:10), labels = c("Region 1:\nCT,ME,MA,\nNH,RI,VT", "Region 2:\nNJ,NY", "Region 3:\nDE,MD,PA,\nVA,WV", "Region 4:\nAL,FL,GA,\nKY,MS,NC,\nSC,TN", "Region 5:\nIL,IN,MI,\nMN,OH,WI", "Region 6:\nAR,LA,NM,\nOK,TX", "Region 7:\nIA,KS,MO,NE", "Region 8:\nCO,MT,ND,\nSD,UT,WY", "Region 9:\nAZ,CA,NV", "Region 10:\nAK,ID,OR,WA"))) %>%
  dplyr::mutate(dummy1 = -.03, dummy2=-.1, dummy3 = -.07)

####################################
## Plot peak vs cumulative pop ##

plt <- ggplot(peakDat, aes(x = plt.date_cty, y = cumpop_prop_cty, group = fips_st)) +
  geom_line(aes(colour = plt.st.color)) +
  scale_x_date("Peak Week", date_labels = "%b", date_breaks = "2 months") +
  scale_y_continuous("Cumulative proportion of the population", breaks = c(.1, .5, .9)) +
  scale_colour_gdocs() +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size=10), strip.text.y = element_text(size=8, angle=360), panel.spacing = unit(0, "lines")) +
  guides(colour = "none") +
  facet_grid(plt.region~plt.season)

ggsave("graph_outputs/explore_irDt_cumpop/peak_cumpop_region_season_grid.png", plt, width = 7, height = 9)



plt <- ggplot(peakDat, aes(x = plt.date_cty, y = cumpop_prop_cty, group = fips_st)) +
  geom_line(aes(colour = plt.st.color)) +
  geom_segment(aes(x = plt.date_st, xend = plt.date_st, y = dummy1, yend = dummy2, colour = plt.st.color)) +
  geom_point(aes(x = plt.date_reg, y = dummy3), colour = "black", shape = 1) +
  scale_x_date("Peak Week", date_labels = "%b", date_breaks = "2 months") +
  scale_y_continuous("Cumulative proportion of the population", breaks = c(.1, .5, .9), limits = c(-.1,1)) +
  scale_colour_gdocs() +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size=10), strip.text.y = element_text(size=8, angle=360), panel.spacing = unit(0, "lines")) +
  guides(colour = "none") +
  facet_grid(plt.region~plt.season)

ggsave("graph_outputs/explore_irDt_cumpop/peak_cumpop_region_season_grid_wScales.png", plt, width = 7, height = 9)

####################################
## % of county population with flu outbreaks at the time when state outbreak is identified

statePeakDat <- peakDat %>%
  dplyr::filter(peakwk_cty == peakwk_st) %>%
  left_join(refs_st, by = c("fips_st"))
seasons <- sort(unique(statePeakDat$season))

for(s in seasons){

  dummyDat <- statePeakDat %>% 
    dplyr::filter(season == s) %>%
    dplyr::arrange(desc(cumpop_prop_cty))
  pltdat <- dummyDat %>%
    dplyr::mutate(abbr_st = factor(abbr_st, levels = dummyDat$abbr_st))
  plt <- ggplot(pltdat, aes(x = abbr_st, y = cumpop_prop_cty)) +
    geom_col() +
    theme_bw() +
    theme(legend.position = "bottom", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous("Cumulative county population (%)\npast peak season at state peak")
  ggsave(paste0("graph_outputs/explore_irDt_cumpop/peak_cumpop_state_S", s, ".png"), plt, width = 7, height = 3)
}

statePeaksummaryDat <- statePeakDat %>%
  group_by(abbr_st) %>%
  summarise(cumpop_prop_cty_mn = mean(cumpop_prop_cty), cumpop_prop_cty_sd = sd(cumpop_prop_cty)) %>%
  dplyr::mutate(cumpop_min = cumpop_prop_cty_mn-cumpop_prop_cty_sd,
                cumpop_max = cumpop_prop_cty_mn+cumpop_prop_cty_sd ) %>%
  dplyr::arrange(desc(cumpop_prop_cty_mn))
pltDat2 <- statePeaksummaryDat %>%
  dplyr::mutate(abbr_st = factor(abbr_st, levels = statePeaksummaryDat$abbr_st))

plt <- ggplot(pltDat2, aes(x = abbr_st, y = cumpop_prop_cty_mn)) +
    geom_linerange(aes(ymin = cumpop_min, ymax = cumpop_max)) +
    geom_point() +
    geom_hline(aes(yintercept = mean(cumpop_prop_cty_mn)), colour = "red") +
    theme_bw() +
    theme(legend.position = "bottom", axis.title.y = element_blank()) +
    scale_y_continuous("Cumulative county population (%)\npast peak season at state peak", limits = c(0,1.1), breaks = c(0,.25,.5,.75,1)) +
    coord_flip()
  ggsave(paste0("graph_outputs/explore_irDt_cumpop/peak_cumpop_state_Sall.png"), plt, width = 2.75, height = 6.5)