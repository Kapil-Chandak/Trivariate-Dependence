# Trivariate-Dependence: Modeling Compound Flooding Risks

This repository contains the R code and datasets used to investigate compound flooding risks in **Bridgeport, Connecticut** by using trivariate copula models, incorporating non-tidal residual sea level, precipitation, and wind data to improve flood risk assessment and understand the causes and relationships between variables responsible for compound flood risks.

## Structure

* **code**: Contains R scripts and analysis files.
    * `data_processing.R`: Pre-processing of datasets and plots.
    * `TrivariateCopula.Rmd`: Implementation of copula models for tail dependence.
* **data**: Source datasets including wind speed, direction, and sea levels.

## Requirements

The analysis is performed using **R (version 4.0+)**. To run the scripts, install the following dependencies:

```r
install.packages(c("copula", "evd", "MASS", "tidyverse", "dplyr"))
