# earthquake-risk
## This repository contains the code used for the empirical analysis of the paper: Earthquake risk embedded in property prices: evidence from five Japanese cities

## data

  Data from various sources have been cleaned and compiled into the data file `individual_data.csv`. This file is too large for the repository and will be stored at another place.
  
  Short run risk data obtained through simulation is stored in `Xpsi-1.csv`. 

## code

  The code is written in `R version 3.4.0`.
  Additionally the package `dplyr` (version 0.7.6) is needed.
  
  The script `main.R` executes the main procedures to generate the results of the three base models in the paper.
  
  The script `supp_funcs.R` contains all the supplementary functions needed to run main.R, for example the gradient and loglikelihood functions.

