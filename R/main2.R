library(Rcpp)
library(RcppArmadillo)
library(coda)
library(deSolve)
library(gridExtra)
library(ggplot2)
library(doParallel)
library(truncnorm)
library(RColorBrewer)
library(grid)
library(reshape2)
library(pryr)
source("~/Documents/R_packages/mcmcJH/R/checking_functions.R")
source("~/Documents/R_packages/mcmcJH/R/cost_functions.R")
source("~/Documents/R_packages/mcmcJH/R/helper_functions.R")
source("~/Documents/R_packages/mcmcJH/R/maths_functions.R")
source("~/Documents/R_packages/mcmcJH/R/mcmc_functions.R")
source("~/Documents/R_packages/mcmcJH/R/mcmc_wrapper_script.R")
source("~/Documents/R_packages/mcmcJH/R/model_functions.R")
source("~/Documents/R_packages/mcmcJH/R/optifix.R")
source("~/Documents/R_packages/mcmcJH/R/plotting_functions.R")
source("~/Documents/R_packages/mcmcJH/R/poisson_error_functions.R")
source("~/Documents/R_packages/mcmcJH/R/priors.R")
source("~/Documents/R_packages/mcmcJH/R/mcmc_wrapper_2.R")
sourceCpp("~/Documents/R_packages/mcmcJH/R/rcpp_functions.cpp")

setwd("~/Documents/R_packages/mcmcJH/R/tmp6")

data_file <- "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/ferret_panama_split.csv"
mcmc_param_file <- "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/mcmc_parameters.csv"
infection_file <- "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/infection_times.csv"
infection_times <- read.csv(infection_file)

param_files <- list(
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 1 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 2 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 3 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 4 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 5 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 1 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 2 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 3 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 4 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 5 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 1 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 2 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 3 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 4 (H1N1)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 5 (H1N1)_params.csv"
    
    )

consolidate <- list(
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 1 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 1 (H3N2)_params.csv",
   #' "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 1 (H3N2)_params.csv",
 "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 2 (H3N2)_params.csv",
     "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 2 (H3N2)_params.csv",
  #'  "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 2 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 3 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 3 (H3N2)_params.csv",
  #'  "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 3 (H3N2)_params.csv",
"~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 4 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 4 (H3N2)_params.csv",
  #'  "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 4 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 5 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 5 (H3N2)_params.csv"
  #'  "~/Documents/R_packages/mcmcJH/R/tmp5/inputs/Grp 5 (H3N2)_params.csv"
    )
    
 main_dir <- paste(getwd(),"/outputs",sep="")

greb <- MCMC_main(data_file, param_files, mcmc_param_file, main_dir, infection_times, consolidate)
