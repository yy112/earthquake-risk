# main.R
# This is the main program used to generate the results in paper:
# Earthquake risk embedded in property prices: evidence from five Japanese cities

rm(list = ls())
library(dplyr)
source('supp_funcs.R')

# read data 
data <- read.table("individual_data.csv", 
                   sep = ',', header=TRUE, stringsAsFactors = F)

colName.i <- "Area.Ward.City"
colName.t <- "t"
colName.p <- "Type"
# names of the columns of regressors X 
Xnames_all <- c("constant_LandBldg", "constant_LandOnly", "constant_Condo",
                 "distance.num", "area.m2.num", "total.floor.area.m2.num",
                 "building.age", 
                 "LandBldg_RC", "LandBldg_S", "LandBldg_W",
                 "built.1981_2000", "built.after2000" ,
                 "Urban_Control", 
                "RC", "SRC", "RC_SRC", "S", "W", 
                "LU_Resid", "LU_Comm", "LU_Industr",
                "Region_Residential", "Region_Commercial", 
                "Region_Industrial","Region_PotResidential",
                 "max.building.coverage.ratio", "max.floor.area.ratio",
                 "City_Fukuoka", "City_Nagoya", "City_Osaka", "City_Sapporo",
                 "log.nGDP", "log.CPI",  "Int_rate", "log.TOPIX",
                 "PctImmi", "Ncrime", "PctUnemploy", "PctExec",
                "PctForeign", "Ndaycare", "Nkindergtn", "Nagedhome","Nhosp",                  
                "Nlargeretail", "Ndepstore",      
                 "JSHIS_I45", "JSHIS_I55", "JSHIS_I45_55",
                "JSHIS_I45_station", "JSHIS_I55_station", "JSHIS_I45_55_station",
                "Xpsi_obj",
                "Q1", "Q2", "Q3", "Q4", "Q123", 
                "Q_after_Fukushima", "age_W")

Xnames_base <- c("constant_LandBldg", "constant_LandOnly", "constant_Condo",
            "distance.num", "area.m2.num", "total.floor.area.m2.num",
            "building.age", 
            "LandBldg_RC", "LandBldg_S", "LandBldg_W",
            "built.1981_2000", "built.after2000" ,
            "Urban_Control",  
            "max.building.coverage.ratio", "max.floor.area.ratio",
            "City_Fukuoka", "City_Nagoya", "City_Osaka", "City_Sapporo",
            "log.nGDP", "log.CPI",  
            "PctImmi", "Ncrime", "PctUnemploy", "PctExec",
            "JSHIS_I45_55", "JSHIS_I55", "Xpsi")



# generate averages of y, X and cell counts H
data_vec <- vectorize(data = data, colName.i = colName.i, 
                 colName.t = colName.t, colName.p = colName.p,
                 colName.y = "log.price",
                 colName.X = Xnames_all)
# retrieve each component of the results
H <- data_vec$H
y <- data_vec$y
X <- data_vec$X
district <- data_vec$district
time <- data_vec$time
type <- data_vec$type

# choose the theta parameters to include
# in the 2 error component case: 6 parameters for the zeta matrix and 6 for the epsilon matrix
include_2error <- c(rep(1, 6), rep(0,6), rep(1,6))
# choose initial parameters for optimization
initpar <- c(-1.5, 0.1, -0.1, -2.1, -0.01, -1.1,
             -2.5, 0.02, -0.1, -3.5,  0.1, -3.5,
             -0.9, 0.008, 0.005, -0.9, 0.01, -0.9)
# results for the model with only long run risk variables (model 1)
results_LRonly <- ec_reg(data_X = X, data_y = y, data_H = H,
                  var = setdiff(Xnames_base, "Xpsi"),
                  district = district, time = time, type = type,
                  par.include = include_2error,
                  par.init = initpar)
# results for the model when the short run risk variable is objective (model 2)
results_LR_objSR <- reg_psi(data_X=H, data_y=y, data_H=H,
                        psi = 1, method = "p",
                        district = district, time = time, type = type,
                        var = Xnames_base,
                        par.include = include_2error,
                        par.init = initpar)
# 
# results for the base model (model 3)
results_base <- reg_psi(data_X=X, data_y=y, data_H=H,
                    psi = 3.74, method = "p",
                    district = district, time = time, type = type,
                    var = Xnames_base,
                    par.include = include_2error,
                    par.init = initpar)

#################################################################################
# To find the value of psi that maximizes likelihood, 
# perform a grid search over possible values between 0.1 to 10.
# Refine the grid by each iteration to find the final value (precise up to 2 digits) 

list0 <- seq(from = 0.5, to = 10, by = 0.5)
list1 <- seq(from = 3.65, to = 3.85, by = 0.01)
psi_optim <- opt_psi(data_X=X, data_y=y, data_H=H, 
        list_psi = list1, method = "p",
        district = district, time = time, type = type,
        var = Xnames_base,
        par.include = include_2error,
        par.init = initpar,
        plot = TRUE)
