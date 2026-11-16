# Trivariate-Dependence: Modeling Compound Flooding Risks

This repository contains the R code and datasets used to investigate compound flooding risks in **Bridgeport, Connecticut** by using trivariate copula models, incorporating non-tidal residual sea level, precipitation, and directional wind components (u_t, v_t) to improve flood risk assessment and understand the causes and relationships between variables responsible for compound flood risks.

## Structure

* **code**: Contains R scripts and analysis files.
    * `data_processing.R`: Pre-processing of datasets and plots.
    * `Rubins.R`: Fits five copula families to each triad, with dry-day precipitation zeros jittered 50 times and results pooled by Rubin's rules.
* **additional_analysis**: Adequacy QQ plots, declustering, and missing-data imputation checks.
* **data**: Source datasets including wind speed, direction, precipitation, and sea levels.

## Requirements

The analysis is performed using **R (version 4.0+)**. To run the scripts, install the following dependencies:

```r
install.packages(c("copula", "evd", "MASS", "tidyverse", "dplyr"))
```

## Citation

Chandak, K., Ren, Z., de Vos, M., Ravishanker, N., & O'Donnell, J. (2026). Co-Occurrence of Extreme Precipitation, Sea Level and Wind: Statistical Perspective. *Journal of Marine Science and Engineering*, 14, 2056. https://doi.org/10.3390/jmse14212056
