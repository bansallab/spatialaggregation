## 8/1/2021
## Draw scatterplot of county pop vs onset or peak timing, colored by whether counties precede or succeed state onset or peak timing
## One panel by season

library(tidyverse)
library(colorspace)

ctyDat <- read_csv("../R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_cty.csv")
stDat <- read_csv("../R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_st.csv")

ctyOnset <- ctyDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  dplyr::group_by(fips, season) %>%
  dplyr::filter(Thu.week == min(Thu.week)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fips_st = str_sub(fips, 1, 2)) %>%
  dplyr::rename(onsetwk_cty = Thu.week) %>%
  dplyr::select(fips, fips_st, season, onsetwk_cty, pop)

stOnset <- stDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  dplyr::group_by(fips_st, state, season) %>%
  dplyr::filter(Thu.week == min(Thu.week)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(onsetwk_st = Thu.week) %>%
  dplyr::select(fips_st, state, season, onsetwk_st)

ctyPeak <- ctyDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  dplyr::group_by(fips, season) %>%
  dplyr::filter(ir.dt == max(ir.dt)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fips_st = str_sub(fips, 1, 2)) %>%
  dplyr::rename(peakwk_cty = Thu.week) %>%
  dplyr::select(fips, fips_st, season, peakwk_cty, pop)

stPeak <- stDat %>%
  dplyr::filter(season > 2 & has.epi & incl.analysis & in.season) %>%
  dplyr::group_by(fips_st, state, season) %>%
  dplyr::filter(ir.dt == max(ir.dt)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(peakwk_st = Thu.week) %>%
  dplyr::select(fips_st, state, season, peakwk_st)

onsetDat <- dplyr::left_join(ctyOnset, stOnset, by = c("fips_st", "season")) %>%
    dplyr::mutate(plt.season = factor(season, levels = 3:9, labels = c("2002-2003", "2003-2004", "2004-2005", "2005-2006", "2006-2007", "2007-2008", "2008-2009"))) %>%
    dplyr::mutate(relative_to_st = ifelse(onsetwk_cty < onsetwk_st, "precedes", ifelse(onsetwk_cty == onsetwk_st, "matches", "succeeds"))) %>%
    dplyr::mutate(relative_to_st = factor(relative_to_st, levels = c("precedes", "matches", "succeeds"))) %>%
    dplyr::filter(!is.na(relative_to_st))
peakDat <- dplyr::left_join(ctyPeak, stPeak, by = c("fips_st", "season")) %>%
    dplyr::mutate(plt.season = factor(season, levels = 3:9, labels = c("2002-2003", "2003-2004", "2004-2005", "2005-2006", "2006-2007", "2007-2008", "2008-2009"))) %>%
    dplyr::mutate(relative_to_st = ifelse(peakwk_cty < peakwk_st, "precedes", ifelse(peakwk_cty == peakwk_st, "matches", "succeeds"))) %>%
    dplyr::mutate(relative_to_st = factor(relative_to_st, levels = c("precedes", "matches", "succeeds"))) %>%
    dplyr::filter(!is.na(relative_to_st))

pal <-sequential_hcl(3, palette = "Hawaii")

## figures by onset
onsetScatter1 <- ggplot(onsetDat, aes(x = onsetwk_cty, y = pop, group = season)) +
    geom_point(aes(colour = relative_to_st), alpha = 0.7) +
    scale_y_log10("County Population") +
    scale_x_date("County Onset Week", date_breaks = "1 month", date_labels = "%b") +
    scale_colour_manual("Relative to state onset week", values = pal) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_blank()) +
    facet_wrap(~state, scales = "free_x")
onsetScatter2 <- ggplot(onsetDat, aes(x = onsetwk_cty, y = pop, group = season)) +
    geom_point(aes(colour = relative_to_st), alpha = 0.7) +
    scale_y_log10("County Population") +
    scale_x_date("County Onset Week", date_breaks = "1 month", date_labels = "%b") +
    scale_colour_manual("Relative to state onset week", values = pal) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~plt.season, scales = "free_x")

ggsave("../graph_outputs/explore_timing_pop_cty/onset_scatter_st.png", onsetScatter1, width = 7, height = 9)
ggsave("../graph_outputs/explore_timing_pop_cty/onset_scatter_seas.png", onsetScatter2, width = 7, height = 9)

## figures by peak
peakScatter1 <- ggplot(peakDat, aes(x = peakwk_cty, y = pop, group = season)) +
    geom_point(aes(colour = relative_to_st), alpha = 0.7) +
    scale_y_log10("County Population") +
    scale_x_date("County Peak Week", date_breaks = "1 month", date_labels = "%b") +
    scale_colour_manual("Relative to state peak week", values = pal) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_blank()) +
    facet_wrap(~state, scales = "free_x")
peakScatter2 <- ggplot(peakDat, aes(x = peakwk_cty, y = pop, group = season)) +
    geom_point(aes(colour = relative_to_st), alpha = 0.7) +
    scale_y_log10("County Population") +
    scale_x_date("County Peak Week", date_breaks = "1 month", date_labels = "%b") +
    scale_colour_manual("Relative to state peak week", values = pal) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~plt.season, scales = "free_x")

ggsave("../graph_outputs/explore_timing_pop_cty/peak_scatter_st.png", peakScatter1, width = 7, height = 9)
ggsave("../graph_outputs/explore_timing_pop_cty/peak_scatter_seas.png", peakScatter2, width = 7, height = 9)

#### figures by season
seasons <- 3:9
for(s in seasons){
    dummyOnset <- dplyr::filter(onsetDat, season == s)
    onsetScatter3 <- ggplot(dummyOnset, aes(x = onsetwk_cty, y = pop, group = interaction(season, state))) +
        geom_point(aes(colour = relative_to_st), alpha = 0.7) +
        geom_smooth(colour = "grey30", method = lm) +
        scale_y_log10("County Population") +
        scale_x_date("County Onset Week", date_breaks = "1 month", date_labels = "%b") +
        scale_colour_manual("Relative to state onset week", values = pal) +
        theme_bw() +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~state, scales = "free_y")
    ggsave(glue::glue("../graph_outputs/explore_timing_pop_cty/onset_scatter_st_S{s}.png"), onsetScatter3, width = 9, height = 8)

    dummyPeak <- dplyr::filter(peakDat, season == s)
    peakScatter3 <- ggplot(dummyPeak, aes(x = peakwk_cty, y = pop, group = interaction(season, state))) +
        geom_point(aes(colour = relative_to_st), alpha = 0.7) +
        geom_smooth(colour = "grey30", method = lm) +
        scale_y_log10("County Population") +
        scale_x_date("County Peak Week", date_breaks = "1 month", date_labels = "%b") +
        scale_colour_manual("Relative to state peak week", values = pal) +
        theme_bw() +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~state, scales = "free_y")
    ggsave(glue::glue("../graph_outputs/explore_timing_pop_cty/peak_scatter_st_S{s}.png"), peakScatter3, width = 9, height = 8)
}