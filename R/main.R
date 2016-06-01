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
setwd("~/Documents/ferret_results")

data <- matrix(ncol=3,nrow=15)
data[,1] <- c(0,0,0,21,21,21,36,36,36,49,49,49,70,70,70)
data[,2] <- c(0,0,0,0,0,0,6,6,7,13,13,12,12,11,11)
data[,2] <- c(0,0,0,0,0,0,4,5,4,4,4,5,9,10,10)
data[,2] <- c(0,0,0,0,0,0,0,0,0,0,0,0,5,5,4)
data[,2] <- c(0,0,0,8,8,8,3,6,4,5,5,4,4,4,2)
data[,2] <- c(0,0,0,8,8,6,8,8,7,11,11,6,7,8,7)
data[,2] <- c(0,0,0,8,8,8,8,8,8,8,8,7,6,5,5)
data <- predict_titres(c(-1000,10,5,0.3,12,25,0.01,80,8,0.7,12,30,0.02),seq(0,200,by=5))
data <- predict_titres(c(-1000,0,8,0.5,12,25,0.03,70,10,0.6,12,25,0.03),seq(0,200,by=20))
                       c(0,28,37,49,70))

parameters <- c(-1000,10,8,0.5,15,25,0.02)
data <- predict_titres(parameters,seq(0,200,by=25))
data[,2] <- floor(data[,2])
plot(data)
param_table <- param_table[1:7,]


data[,2] <- as.integer(data[,2])
param_table <- as.matrix(read.csv("~/Documents/R_packages/mcmcJH/R/param_table.csv",header=0))
start <- c(0.9,0.05,-1000,28,6,0.5,12,25,0.02,42,10,0.65,12,55,0.03,56,10,0.65,12,55,0.03)
start <- c(0.8,0.05,-1000,0,4,0.8,12,30,0.03,56,12,0.65,12,25,0.01)
start <- c(0.9,0.05,-1000,0,4,0.8,12,30,0.03,28,6,0.5,12,25,0.02,42,10,0.65,12,55,0.03,56,10,0.65,12,55,0.03)
start <- grebbons
start <- c(0.9,-1000,0,5,0.9,15,31,0.05)
win <- run_MCMC(start,data,param_table,500000,0.44,10000,1,200000,200000,"test",500)

greb1 <- read.csv(win)
greb1 <- greb1[400000:nrow(greb1),]
#'plot(as.mcmc(greb1))
c1 <- cov(greb1[,2:(ncol(greb1)-1)])

win1 <- run_MCMC_test(start,data,param_table,1000000,0.234,10000,5,1000000,1000000,"test1",1000,c1,c(1,2.38,2.38,2.38,2.38),0.8)
greb <- read.csv(win1)
greb <- greb[seq(400000,nrow(greb),by=1),]
plot(as.mcmc(greb))


pars <- as.numeric(greb1[which.max(greb1[,ncol(greb1)]),4:(ncol(greb1))])
pars <- NULL
for(i in 4:(ncol(greb1))){
    pars[i-3] <- mean(greb1[,i])
}


grebbons <- c(-1000,10,5,0.3,12,25,0.01,80,8,0.7,12,30,0.02)
y <- generate_prediction_intervals_1.2(list(greb1[,4:(ncol(greb1)-1)]),0,data[,1],100000,0.95,0.01,predict_titres)
x <- predict_titres(pars,seq(0,70,by=1))
plot(data,ylim=c(0,13),xlim=c(0,70))
lines(x)
lines(y$upper$y~y$upper$x)
lines(y$lower$y~y$lower$x)
curwd <- getwd()

discrete_norm <- function(a, b, sd){
    c <- pnorm(a+0.5,b,sd) - pnorm(a-0.5,b,sd)
    return(c)
}

cppFunction('double discrete_norm2(double a, double b, double sd){
double c;
c = R::pnorm((a+0.5),b,sd,1,0) - R::pnorm((a-0.5),b,sd,1,0);
return(c);
}')

#' Run 1
data_file <- "/inputs/inputs_obs/ferret_panama_split.csv"
mcmc_param_file <- "/inputs/inputs_obs/mcmc_parameters.csv"
infection_file <- "/inputs/inputs_obs/infection_times.csv"
output_file <- NULL
#'main_dir <- paste(getwd(),"/outputs/outputs_tp",sep="")
 main_dir <- "~/Documents/ferret_results/mcmc_outputs/outputs_obs"
VERBOSE <- TRUE

#' Data to be fitted
raw_data <- read.csv(paste(getwd(),data_file,sep=""),stringsAsFactors=FALSE,na.strings=c("NA","-","?"),header=0)
colnames(raw_data) <- raw_data[1,]
raw_data <- raw_data[2:nrow(raw_data),]

#' This project specific - getting the experimental groups such that model fitting can be done for each group sequentially
melted.data <- melt(raw_data)
melted.data$variable <- as.integer(as.character(melted.data$variable))

#' Infection times
infection_table <- read.csv(paste(getwd(),infection_file,sep=""),header=1)
groups <- unique(melted.data$group)
groups <- c("Grp 1 (H3N2)", "Grp 2 (H3N2)", "Grp 3 (H3N2)", "Grp 4 (H3N2)","Grp 5 (H3N2)")
groups <- c("Grp 1 (H3N2)")
ALL_RESULTS <- NULL
index <- 1

#' For each group
for(group in groups){
    param_file <- paste("/inputs/inputs_final/",group,"_params.csv",sep="")
    ALL_RESULTS[[index]]<- MCMC_main(data_file,param_file,mcmc_param_file,output_file,main_dir,group,VERBOSE,infection_table)
    index <- index + 1
}

cppFunction('NumericVector subset_test(NumericVector a, int low, int high){
  return(a[a >= low & a < high]);
}
')
