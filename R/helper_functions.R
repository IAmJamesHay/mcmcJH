#' Parameter table loading function
#'
#' Reads in a .csv file to generate a table of parameter data for use in MCMC. The table should be have the following columns:
#' \itemize{
#' \item{name: }{Name of the parameter (string/character)}
#' \item{value: }{arbritary numeric value that represents a typical value (this is not that important - it is just used for some upkeep functions)}
#' \item{lower_bound: }{numeric value for the lower allowable bound of the parameter (this can be -Inf)}
#' \item{upper_bound: }{numeric value for the upper allowable bound of the parameter (this can be Inf)}
#' \item{use_logistic: }{boolean value (1/0) indicating whether or not the parameter should be logistic transformed during optimisation, ensuring that the parameter search space is between the specified bounds}
#' \item{use_log: }{boolean value (1/0) indicating whether parameter optimisation should be on a log scale (note that this and use_logistic should NOT both be specified}
#' \item{step: }{initial step size for the MCMC random walk}
#' \item{fixed: }{boolean value indicating whether parameter should be fixed during MCMC/optimisation or not (1/0)}
#' \item{nuisance: }{boolean value indicating if parameter is a nuisance paramter (1/0)}
#' \item{lower_seed: }{numeric value for the lower bound of the random seed generation}
#' \item{upper_seed: }{numeric value for the upper bound of the random seed generation}
#' \item{prior_func: }{this should be a string that matches the name of a function that will be used as the prior function for MCMC. This function should be available in the current R environment, or an error will be thrown. See below.}
#' \item{prior_args: }{a vector of default arguments for the prior function. For example, this might contain the standard deviation or mean of a normal distribution. Most prior functions have a flag to indicate log probabilities, which should be set to true by this argument}
#' \item{log_proposal: }{boolean flag to indicate whether proposals in the MCMC function should be made on a log or linear scale (1 for log, 0 for linear). If in doubt, set this to 0.}
#' }
#' @param param_file the full file location of the .csv file to be read
#' @return a data frame containing the formatted data imported from the given .csv file
#' @seealso \code{\link{param_check}}
#' @export
load_param_table <- function(param_file){
       # Parameters for model
    param_table <- read.csv(paste(getwd(),param_file,sep=""),header=1)
    # Make some slight formatting adjustments to data and parameter tables
    # Parameters
    colnames(param_table) <- c("names","value","lower_bound","upper_bound","use_logistic","use_log","step","fixed","nuisance","lower_seed","upper_seed","prior_func","prior_args","log_proposal")
    param_table$names <- as.character(param_table$names)
    
    # Priors - match up prior functions and arguments
    param_table$prior_func<- as.character(param_table$prior_func)
    param_table$prior_args <- as.character(param_table$prior_args)
    param_table$prior_func <- lapply(param_table$prior_func, function(x) match.fun(x))
    param_table$prior_args <- lapply(param_table$prior_args, function(x) as.numeric(unlist(strsplit(x,","))))
    return(param_table)
}

#' Parameter logit/log transform function
#' 
#' Takes a vector of parameters and a table of parameter bounds, transformation instructions etc. Idea is to convert all parameters to log/specified logit scale.
#' Table should have four columns: use_logistic, use_log, lower_bound, and upper_bound.
#' @param params vector of parameters to be transformed
#' @param param_table a matrix of boolean flags and bounds
#' @return a vector of the transformed parameters
#' @seealso \code{\link{param_check}}, \code{\link{transform_params_logistic}}
#' @export
transform_params_logit<- function(params, param_table){
    tmp <- params
    use_logistic <- param_table[,1]
    use_logs <- param_table[,2]
    lower_bound <- param_table[,3]
    upper_bound <- param_table[,4]
    for(i in 1:nrow(param_table)){
        tmp[i] <- ifelse(use_logistic[i],transform_logit(params[i],lower_bound[i],upper_bound[i]),
                         ifelse(use_logs[i], log(params[i]), params[i]))
    }
    return(tmp)
}

#' Parameter logistic/exp transform function
#' 
#' Takes a list of parameters and table of parameter bounds. Idea is to convert all parameters from log/specified logit scale to normal scale (-Inf to +Inf).
#' Table should have four columns: use_logistic, use_log, lower_bound, and upper_bound.
#' @param params vector of parameters to be transformed
#' @param param_table a matrix of boolean flags and bounds
#' @return a vector of the transformed parameters
#' @seealso \code{\link{param_check}}, \code{\link{transform_params_logit}}
#' @export
transform_params_logistic<- function(params, param_table){
    tmp <- params
    use_logistic <- param_table[,1]
    use_logs <- param_table[,2]
    lower_bound <- param_table[,3]
    upper_bound <- param_table[,4]
    for(i in 1:nrow(param_table)){
        tmp[i] <- ifelse(use_logistic[i],transform_logistic(params[i],lower_bound[i],upper_bound[i]),
                         ifelse(use_logs[i], exp(params[i]), params[i]))
    }
    return(tmp)
}

#' Random parameter generator
#' 
#' From a table of parameters, generates random values within the given bounds. Idea is to use these as seeds for the MCMC chains
#' Table should be as imported by \code{\link{load_param_table}}.
#' @param param_table a data frame of parameters containing columns with the lower and upper seed bounds (titled lower_seed and upper_seed)
#' @return a vector of random parameters
#' @seealso \code{\link{load_param_table}}
#' @export
rand.params <- function(param_table){
    params <- NULL
    for(i in 1:nrow(param_table)){
        if(param_table$fixed[i]==0){
            params[i] <- runif(1,param_table$lower_seed[i],param_table$upper_seed[i])
        }
        else{
            params[i] <- param_table$value[i]
        }
    }
    return(params)
}

#' Condenses a table of data to provide mean and standard deviations for each time point
#' 
#' Takes a data frame of uncondensed data (ie. individuals rather than mean) and returns the mean and standard deviation
#' of the data. Given data should have 3 levels of grouping (eg. group, strain and individual), with "individual" being essential.
#' Should also have columns for "variable" and "value", as returned when melting data.
#' @param data a data frame containing the data to be condensed
#' @return a matrix of the condensed data
#' @export
condense.data <- function(data){
    groupings <- colnames(data)
    groupings <- groupings[!(groupings %in% c("individual","value"))]
    temp <- unique(data[,groupings])
    for(i in unique(data[,groupings[1]])){
        for(j in unique(data[,groupings[2]])){
            for(q in unique(data[,groupings[3]])){
                temp$value[temp[,groupings[1]] == i & temp[,groupings[2]] == j & temp[,groupings[3]] == q] <- mean(data[data[,groupings[1]] == i & data[,groupings[2]] == j & data[,groupings[3]] == q,"value"])
                temp$sd[temp[,groupings[1]] == i & temp[,groupings[2]] == j & temp[,groupings[3]] == q] <- sd(data[data[,groupings[1]] == i & data[,groupings[2]] == j & data[,groupings[3]] == q,"value"])
            }
        }
    }
    return(temp)
}

#' Condenses data from the same times points into mean values
#'
#' Does what it says on the tin
#' @param data a data frame of the data to be condensed. The first column contains all of the time points (which may/should be duplicates)
#' @return a matrix of condensed data
#' @export
condense.data.simple <- function(data){
    times <- unique(data[,1])
    tmp <- data.frame(time=times,y=NA)
    for(i in times){
        tmp[tmp$time==i,2] <- mean(data[data[,1] == i,2])
    }
    return(tmp)
}

#' Second tier function for generating prediction intervals from the MCMC chain
#'
#' This function is called by a wrapper function, \code{\link{generate_prediction_intervals_1.1}} to generate plottable prediction intervals for a model plot.
#' @param mcmc_chains a matrix containing the whole MCMC chain (ie. parameter values at each iteration)
#' @param burnin optional parameter specifying the number of iterations from the MCMC chain that should be discarded as burn in
#' @param times a vector of time points to calculate the model trajectory over
#' @param runs number of random draws from the multivariate posterior with which to calculate the prediction intervals. Defaults to 100
#' @param level desired prediction interval size. Defaults to 0.95 ie. 95% prediction intervals
#' @param smoothing uses smoothing splines to give smoother prediction intervals where these might be aethsetically "spikey". Defaults to 0
#' @param MODEL_FUNCTION a pointer to an existing R function that is used to calculate the model trajectories
#' @return a list of the lower and upper prediction lines
#' @seealso \code{\link{generate_prediction_intervals_1.1}}
#' @export
generate_prediction_intervals_1.2 <- function(mcmc_chains,burnin=NULL,times,runs=100,level=0.95,smoothing=0,MODEL_FUNCTION){
    upper_level <- level+((1-level)/2)
    lower_level <- (1-level)/2

    if(!is.null(burnin)) burnin_start<- burnin
    else burnin_start <- 1

    # Combine chains so can sample from all
    mcmc <- NULL
    for(i in 1:length(mcmc_chains)){
        mcmc <- rbind(mcmc, mcmc_chains[[i]][burnin_start:nrow(mcmc_chains[[i]]),])
    }

    # Get sample indices
    samples <- sample(1:nrow(mcmc),runs,replace=T)
    
    times <- seq(times[1],max(times),by=1)

    max_time <- length(times)
    temp_matrix <- matrix(nrow=max_time,ncol=runs)

    # For each sample from the multivariate posterior, get the trajectory
    for(i in 1:runs){
        params <- as.numeric(mcmc[samples[i],])
        temp_matrix[,i] <- MODEL_FUNCTION(params,times)[,2]
    }
    
    upper <- lower <- matrix(nrow=max_time,ncol=2)
    upper[,1] <- lower[,1] <- MODEL_FUNCTION(params,times)[,1]
    # Go through the model trajectories and for each time index, get the upper and lower quantiles
    for(i in 1:max_time){
        upper[i,2] <- quantile(temp_matrix[i,],upper_level)
        lower[i,2] <- quantile(temp_matrix[i,],lower_level)
    }

    # Smoothing of upper and lower bounds for aesthetic reasons.
    upper_smooth <- predict(smooth.spline(upper,spar=smoothing))
    lower_smooth <- predict(smooth.spline(lower,spar=smoothing))
   
    return(list(lower=lower_smooth,upper=upper_smooth))
}

#' Prediction interval generating function
#'
#' Generates data and prediction intervals that are used to plot model fits with prediction intervals and original data from provided MCMC chain file and parameters. Calls a second helper function, \code{\link{generate_prediction_intervals_1.2}}
#' @param results_MCMC a full, valid file path for a .csv file containing the full MCMC chain
#' @param burnin numeric value for the number of MCMC chain iterations to be discarded as burn in
#' @param param_table a correctly formatted parameter table as loaded by \code{\link{load_param_table}}
#' @param all_data a data frame containing all of the original data
#' @param MODEL_FUNCTION a pointer to an existing R function that is used to calculate the model trajectories
#' @return a list containing the following items:
#' \itemize{
#' \item{model_data}{a data frame with the original data (ie. fitted) for use by a plotting function}
#' \item{prediction_data}{a data frame with the best fitting model's data for use by a plotting function}
#' \item{lower_prediction_bounds}{a data frame with the necessary information required to plot the lower prediction bound}
#' \item{upper_prediction_bounds}{a data frame with the necessary information required to plot the upper prediction bound}
#' }
#' @seealso \code{\link{generate_prediction_intervals_1.2}}
#' @export
generate_prediction_intervals_1.1 <- function(results_MCMC, burnin, param_table, all_data,MODEL_FUNCTION){
    upper_bounds <- NULL
        lower_bounds <- NULL
        predictions <- NULL
        index <- 1
        # Generate prediction intervals for each fitted strain
        for(j in 1:length(results_MCMC)){
            if(is.null(results_MCMC[[j]]) | length(results_MCMC[[j]]) < 1){
                next
            }
            tmp_chains <- NULL
            for(z in 1:length(results_MCMC[[j]]$files)){
                tmp_chains[[z]] <- read.csv(results_MCMC[[j]]$files[[z]], header=1)
                tmp_chains[[z]] <- tmp_chains[[z]][,2:(nrow(param_table)+1)]
            }
                                        # Generate model predictions for all fits
                                        # Generate prediction intervals
            bounds <- generate_prediction_intervals_1.2(tmp_chains,
                                                        burnin,
                                                        results_MCMC[[j]]$times,
                                                        runs=10000,
                                                        level=0.95,
                                                        smoothing=0.01,
                                                        MODEL_FUNCTION
                                                    )
            upper_bounds[[index]] <- bounds$upper
            lower_bounds[[index]] <- bounds$lower
            predictions[[index]] <- as.data.frame(MODEL_FUNCTION(results_MCMC[[j]]$params,
                                                                     seq(results_MCMC[[j]]$times[1],max(results_MCMC[[j]]$times),by=0.1)))
            colnames(predictions[[index]]) <- c("variable","value")
            lower_bounds[[index]]$strain <- upper_bounds[[index]]$strain <- predictions[[index]]$strain <- unique(all_data$strain)[j]
            lower_bounds[[index]]$group <- upper_bounds[[index]]$group <- predictions[[index]]$group <- unique(all_data$group)[j]
            index <- index + 1
        }
        all_model_dat <- NULL
        for(j in 1:length(predictions)){
            if(length(predictions) >0 & nrow(predictions[[j]]) > 0){
                all_model_dat <- rbind(all_model_dat, predictions[[j]])
            }
        }

        # Plot model fits
    temp_all <- all_data[all_data$group == group,]
    
    return(list("model_data"=all_model_dat,"prediction_data"=temp_all,"lower_prediction_bounds"=lower_bounds,"upper_prediction_bounds"=upper_bounds))
}
    


#' Pads lists of different sizes
#'
#' When passed a list of lists of different lengths, returns the same list such that each list is the same size, with NAs used to pad the smaller lists. ie. all lists same size as biggest list
#' @param my_list list of lists to be padded
#' @return the padded list of lists
#' @export
pad_lists <- function(my_list){
    max_length <- max(sapply(1:length(my_list),function(x) length(my_list[[x]])))
    for(i in 1:length(my_list)){
        my_list[[i]] <- c(my_list[[i]], rep(NA, max_length-length(my_list[[i]])))
    }
    return(my_list)
}

#' Remove zero value groups
#'
#' Removes those sub-groups from the data that have only zero values for all time points, and the corresponding strain name from the strain name list
#' @param all_data a data frame to be processed, with time points, values, and a column "strain" that represent the sub-groups
#' @param file_strains a vector of the strain names for the all_data "strain" column
#' @return the same data table and vector of strain names, less those that have all zero values
#' @export
remove_zero_data <- function(all_data,file_strains){
    tmp_all_data <- all_data
    for(i in 1:length(unique(all_data$strain))){
        temp_dat <- all_data[all_data$strain==unique(all_data$strain)[i],]
        if(all(temp_dat$value == 0)){
            tmp_all_data <- tmp_all_data[tmp_all_data$strain != unique(tmp_all_data$strain)[i],]
            file_strains <- file_strains[-i]
        }
    }
    all_data <- tmp_all_data
    return(list(all_data,file_strains))
}

    
#' OLD transform parameters to logit
#' 
#' Takes a list of parameters and a table of parameter bounds, transformation instructions etc. Idea is to convert all parameters to log/specified logit scale. 
#' @seealso \code{\link{transform_params_logit}}
transform_params_logitOLD<- function(params, param_table){
    tmp <- params

    for(i in 1:nrow(param_table)){
        if(param_table$use_log[i]) tmp[i] <- log(params[i])
        else if(param_table$use_logistic[i]) tmp[i] <- transform_logit(params[i], param_table$lower_bound[i], param_table$upper_bound[i])
        else tmp[i] <- params[i]
    }
    
    return(tmp)
}

#' OLD transform parameters to logistic
#' 
#'  Takes a list of parameters and table of parameter bounds. Idea is to convert all parameters from log/specified logit scale to normal scale (-Inf to +Inf).
#'@seealso \code{\link{transform_params_logistic}}
transform_params_logisticOLD<- function(params, param_table){
    tmp <- params

    for(i in 1:nrow(param_table)){
        if(param_table$use_log[i]) tmp[i] <- exp(params[i])
        else if(param_table$use_logistic[i]) tmp[i] <- transform_logistic(params[i], param_table$lower_bound[i], param_table$upper_bound[i])
        else tmp[i] <- params[i]
    }
    
    return(tmp)
}

