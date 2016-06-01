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

setwd("~/Documents/R_packages/mcmcJH/R/simulation")

infection_times <- data.frame(group=c("Grp 1","Grp 1","Grp 2","Grp 2","Grp 2","Grp 2","Grp 3","Grp 3","Grp 3","Grp 4","Grp 4","Grp 4","Grp 4","Grp 5","Grp 5","Grp 5"),infection=c("A","D","A","B","C","D","B","C","D","A","B","C","D","B","C","D"),time=c(0,56,0,28,42,56,28,42,56,0,28,42,56,28,42,56))

groups <- list("Grp 1","Grp 2", "Grp 3", "Grp 4", "Grp 5")
strains <- list("A")
times <- seq(0,100,by=20)
parameters <- list(
    c(0.79,0.2,-1000,0,8,0.5,12,5,0.03,56,10,0.4,12,7,0.01),
    c(0.79,0.2,-1000, 0, 10,0.9,12,9,0.05,28,4,0.6,12,5,0.03,42,4,0.7,12,5,0.03,56,3,0.3,12,8,0.04),
    c(0.79,0.2,-1000,28,5,0.7,12,10,0.04,42,2,0.5,12,10,0.03,8,0.8,12,6,0.03),
    c(0.79,0.2,-1000, 0, 10,0.9,12,9,0.05,28,4,0.6,12,5,0.03,42,4,0.7,12,5,0.03,56,3,0.3,12,8,0.04),
    c(0.79,0.2,-1000,28,5,0.7,12,10,0.04,42,2,0.5,12,10,0.03,8,0.8,12,6,0.03)
    )
    
    

param_files <- list(
    "~/Documents/R_packages/mcmcJH/R/tmp3/inputs/Grp 1 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp3/inputs/Grp 2 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp3/inputs/Grp 3 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp3/inputs/Grp 4 (H3N2)_params.csv",
    "~/Documents/R_packages/mcmcJH/R/tmp3/inputs/Grp 5 (H3N2)_params.csv"
    )

raw_data <- generate_artificial_data_big(groups,strains,times,parameters,5)
melted.data <- melt(raw_data)
melted.data$variable <- as.integer(as.character(melted.data$variable))

file_names <- c("A")
names(file_names) <- c("A")


greb <- MCMC_fit_1.1("~/Documents/R_packages/mcmcJH/R/simulation/",melted.data,file_names,param_files,100000,1000,1,100000,100000,1,0.234,0,20,predict_titres,infection_times,FALSE,FALSE)
