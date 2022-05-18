# Spatial aggregation choice in the era of digital and administrative surveillance data

This repository provides the data and source code for the following paper: Lee, Elizabeth C., Ali Arab, Vittoria Colizza, and Shweta Bansal. (in press) "Spatial aggregation choice in the era of digital and administrative surveillance data." PLOS Digital Health.

Note: This repository uses git lfs due to a handful of large files. You may need to install this extension on your computer to pull these large files when cloning the repository. See more details [here](https://git-lfs.github.com/).

## Data (`R_export/`)
Processed weekly time series of influenza-like illness (ILI) intensity from 2001 through 2009. A codebook for the key column names can be found in this folder.
* county scale: `R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_cty.csv`
* state scale: `R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_st`
* HHS region scale: `R_export/fullIndicAll_periodicReg_irDt_Octfit_span0.4_degree2_analyzeDB_reg`

Other data in this folder may be required data input files to run the code or analysis outputs generated by the code.

## Code (`programs/`)
Scripts in this folder were used to perform analyses described in this paper or produce figures in the main text or supplementary material (See `SM_figures/` folder).

## Other Input Files (`reference_data/`)
Files in this folder may be required input files to run the code. They are county, state, and HHS region coordinates, abbreviations, and shapefiles derived from public sources of information.

