#' Run 1
data_file <- "/inputs/inputs_tp/ferret_panama_split.csv"
mcmc_param_file <- "/inputs/inputs_tp/mcmc_parameters.csv"
infection_file <- "/inputs/inputs_tp/infection_times.csv"
output_file <- NULL
#'main_dir <- paste(getwd(),"/outputs/outputs_tp",sep="")
 main_dir <- "~/Documents/ferret_results/mcmc_outputs/outputs_PRIORTEST2"
VERBOSE <- FALSE

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
groups <- c("Grp 3 (H3N2)")
#' For each group
for(group in groups){
    param_file <- paste("/inputs/inputs_tp/",group,"_params.csv",sep="")
    end <- MCMC_main(data_file,param_file,mcmc_param_file,output_file,main_dir,group,VERBOSE,infection_table)
    gc()  
}
