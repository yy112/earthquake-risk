# earthquake-risk
## Introduction 

This repository contains the data and code used for the empirical analysis of the paper: __Earthquake risk embedded in property prices: Evidence from five Japanese cities__. 
Detailed instructions are provided to replicate all tables and figures in this paper. 

## Data

1. `individual_data.csv`

Data from various sources have been cleaned and compiled into the data file `individual_data.csv` (size: 214.7 MB). This file is compressed in the zip file: `individual_data.zip` (size: 11.4 MB). Each row in this file represents one property transaction record, with information on the characteristics of the transaction, characteristics of the property, earthquake probabilities associated with the area that the property is located in, macroeconomic variables observed at the transaction period, and demographic characteristics of the ward that the property is located in. There are also interaction terms and derived variables. 

2. `Xpsi-1.csv`
  
Short run earthquake risk data obtained through simulation is stored in `Xpsi-1.csv` (size: 2KB). This file contains five time series (quarterly frequency) of the short run probabilities, each column corresponding to one of the five cities included in the scope of our analysis.  

3. `city_range.csv`

The file `city_range.csv` (size: 1KB) contains the chosen space window (for the estimation of ETAS models) corresponding to each city. Each row represents one city and the variables used in our analysis are "latMin", "latMax", "lngMin", "lngMax", corresponding to the minimum and maximum values of the latitude and longitude of the chosen space window (we have chosen to use rectangular space windows).


4. `JMA_records.csv`

Earthquake records collected by JMA are stored in `JMA_records.csv` (size: 24.15 MB). Each row of this file contains one earthquake record, with information on the time and location of the record. Detailed description of the records can be found in http://www.data.jma.go.jp/svd/eqev/data/bulletin/data/shindo/format_e.txt.


  For a more detailed description of the data used in our analysis as well as the data collection process, please see the accompanying data document (_Earthquake Risk Embedded in Property Prices: Evidence from Five Japanese Cities - Data Documentation_). Code used for data collection and data cleaning can be provided upon request.

## Code

1. Prerequisites

  `R` version later than 3.6.0 is needed.
  Additionally the packages `dplyr` (version 1.0.0), `PtProcess` (version 3.3-13), `readr`, `R.utils` and `zoo` are needed.

2. `etas_funcs.R`

This file contains the main functions used in the estimation and simulation of ETAS models. We use the `R` package `PtProcess` for the estimation and simulation.
  
3. `R` package `mvecr`

  We have written an R package `mvecr` ("multivariate error components estimation in R") for the main estimation functions to replicate the analysis in the paper. To install this package, download the source file for R package `mvecr` (`mvecr_0.3.0.tar.gz`) from https://github.com/yy112/earthquake-risk and choose "Install packages from Package Archive File (.zip; .tar.gz)" in `R`. A manual of this package (`mvecr_manual.pdf`) is also included in this repository, which is more detailed about the inputs and outputs of each function in this package. A subsample of the `individual_data` dataset is included in this package to be used in examples.

  Detailed installation instructions as well as code used to replicate tables and figures of the paper, can be found in `replication_instructions.Rmd` or `replication_instructions.pdf`.

4. `replication_instructions.Rmd`

  This file contains all the code used for generating `replication_instructions.pdf`. The structure of the code is as follows:

A. Estimation and simulation of the ETAS model

  - Estimation of the ETAS model (replication of the summary statistics and estimation results in Table 48 and Table 50 of the data documentation). Computation time needed: 5 minutes.
  - Simulation of short-run earthquake probabilities, using the estimated ETAS parameters and historical earthquake catalogue. Computation time needed: around 24 hours for 30000 Monte Carlo runs for each city. The output of the simulation is contained in `Xpsi-1.csv`. Replication of Figures 1 and 2 of the paper can be done using the provided data `Xpsi-1.csv` without having to perform the simulation first.

and 

B. Estimation of the multivariate error components regression model

  - Characteristics of the housing dataset. Replication of Table 1 of the paper.
  - Characteristics of the JSHIS long run probabilities. Replication of Table 2 of the paper.
  - Main estimation results. Replication of Table 3 and Figure 3 of the paper. Computation time needed: 1 - 1.5 hours for each model.
  - Sensitivity analysis regarding probability weighting functions. Replication of Table 4 and Figure 4 of the paper. Computation time needed: 1 - 1.5 hours for each model.
  - Sensitivity analysis regarding other model specifications. Replication of Tables B1 - B5 of the supplementary material. Computation time needed: 1 - 1.5 hours for each model with 2 error components, 3 - 4 hours for the 3 error components model, 40 - 60 minutes for the model where individual data is grouped by nearest station instead of by district.
  - Importance ordering and decomposition of risk premia. Replication of Tables C6 - C7 of the supplementary material.

These two parts can be executed independently of each other. The functions related to the estimation and simulation of the ETAS model are contained in the file `etas_funcs.R`. The functions related to the estimation of the multivariate error components regression model are contained in the `R` package `mvecr`.

## Results

The following output files can be found in the folder `output`:

1. results_long_run_only_model.csv
2. results_objective_short_run_model.csv
3. results_base_model_psi3.74.csv
4. results_SR_TK_psi1.40.csv
5. results_LR_prelec_psi3.78_gamma0.17.csv
6. results_LR_TK_psi3.77_gamma0.32.csv
7. results_attr_psi3.75.csv
8. results_noGDP_psi2.63.csv
9. results_UC_psi3.72.csv
10. results_BS_psi3.89.csv
11. results_LandUse_psi3.76.csv
12. results_noTokyo_psi1.9.csv
13. results_noNagoya_psi4.11.csv
14. results_noOsaka_psi4.04.csv
15. results_Q123_psi4.56.csv
16. results_Q4_psi3.89.csv
17. results_Tohoku_psi3.27.csv
18. results_3error_psi3.52.csv
19. results_station_psi3.41.csv

These are the results obtained from running the code contained in `replication_instructions.Rmd`/`replication_instructions.pdf`. 



## References

Ikefuji, M., Laeven, R. J., Magnus, J. R., & Yue, Y. (2020). Earthquake risk embedded in property prices: Evidence from five Japanese cities.

Ikefuji, M., Laeven, R. J., Magnus, J. R., & Yue, Y. (2020). Earthquake risk embedded in property prices: Evidence from five Japanese cities - Supplementary material.

Ikefuji, M., Laeven, R. J., Magnus, J. R., & Yue, Y. (2020). Earthquake risk embedded in property prices: Evidence from five Japanese cities - Data documentation.