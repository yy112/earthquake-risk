# earthquake-risk
## This repository contains the code used for the empirical analysis of the paper: Earthquake risk embedded in property prices: evidence from five Japanese cities

## data

  Data from various sources have been cleaned and compiled into the data file `individual_data.csv`. This file is compressed in the zip file: `individual_data.zip`.
  
  Short run risk data obtained through simulation is stored in `Xpsi-1.csv`. 

## code

  The code is written in `R version 3.4.0`.
  Additionally the package `dplyr` (version 0.7.6), `PtProcess` are needed.
  
  We have written an R package `mvecr` (multivariate error components estimation in R) for the main estimation steps to replicate the analysis in the paper. Detailed installation instructions as well as code used to replicate tables and figures of the paper, can be found in `instructions.Rmd` or `instructions.pdf`.
