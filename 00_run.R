#### Please use this code to run the model fitting
#### Step test run of 400 is a test code that can be used before 04_aggregate_chains.R or to test the model fitting code runs without running all 2 million seen in step 2.1
#### The code in this file runs the files: MCMC_functions.R, parameter_functions.R, set_up_data_for_modelling.R, MODEL_C_odeintr.R and 01_MCMC.R
#### Please do not run the above listed files directly as they cannot function alone.

#### To run the model the package odeintr needs to be installed.
#### This package is not available in CRAN so please download below.
#### To download this package you need rtools and BH.
#### If you already have this package please skip this step.

require(BH)
install.packages("packages/odeintr-master.tar.xz",type="src",repo=NULL)

##also needs Rcpp to be installed
require(Rcpp)
require(tidyverse)
require(data.table)
require(ggplot2)
source("parameter_functions.R")
source("MCMC_functions.R")

###############################################################################
###############################################################################

####step test run of 400
TEMP<- "min"
RUN_MCMC<- TRUE
STARTLOAD_MCMC<- FALSE ##start MCMC on a previous state
SMOOTH_CLIMATE<- TRUE
SMOOTH_DAYS<- 7
MCMCstepsN<- 400 ##number of chain steps
MCMC_STARTLOAD_FILE_NAME<- ""
MCMC_OUTPUT_FILE_NAME<- "Outputs/0.csv"
DATA_FILE_NAME<- 'Data/Data for interpolation.csv'

####step 2.1
TEMP<- "min"
RUN_MCMC<- TRUE
STARTLOAD_MCMC<- FALSE ##start MCMC on a previous state
SMOOTH_CLIMATE<- TRUE
SMOOTH_DAYS<- 7
MCMCstepsN<- 2000000 ##number of chain steps
MCMC_STARTLOAD_FILE_NAME<- ""
MCMC_OUTPUT_FILE_NAME<- "Outputs/2millRESULT.csv"
DATA_FILE_NAME<- 'Data/Data for interpolation.csv'


###############################################################################
###############################################################################
###############################################################################
###############################################################################

date_start<- as.Date("2012-01-01","%Y-%m-%d") -1
time_start <- 0
time_stop <- 485
delta_T<- 0.0625 
tps <- seq(time_start , time_stop , by = delta_T)
length(tps)

source("set_up_data_for_modelling.R")
source('MODEL_C_odeintr.R')  ##this sets up the model in C

##fixed parameters

NH<- 270000
prior_f<- 0.5

##define prior knowledge on the ranges of the parameters that are biologically acceptable
prior_range_K<- c(1,NH*30)
prior_range_eta<- c(0.01, 10)
prior_range_alpha<- c(0.01, 10)
prior_range_gamH<- 1/c(1, 15)
prior_range_sigH<- 1/c(1, 15)
prior_range_T<- c(30*5*1/delta_T, 30*10*1/delta_T)
prior_range_psi<- c(0,1)
prior_range_a<- c(0,1)
prior_range_i<- c(0,30)
prior_range_delta<- c(-30,0)


source("01_MCMC.R")
