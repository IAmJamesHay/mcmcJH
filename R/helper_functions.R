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


# Takes a list of parameters and a table of parameter bounds, transformation instructions etc. Idea is to convert all parameters to log/specified logit scale. 
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

# Takes a list of parameters and table of parameter bounds. Idea is to convert all parameters from log/specified logit scale to normal scale (-Inf to +Inf).
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

# From a table of parameters, generates random values within the given bounds
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
    
# Takes a data frame of uncondensed data (ie. individuals rather than mean) and returns the mean and standard deviation
# of the data. Given data should have 3 levels of grouping (eg. group, strain and individual), with "individual" being essential.
# Should also have columns for "variable" and "value", as returned when melting data.
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

# Condenses data from the same times points into mean values
condense.data.simple <- function(data){
    times <- unique(data[,1])
    tmp <- data.frame(time=times,y=NA)
    for(i in times){
        tmp[tmp$time==i,2] <- mean(data[data[,1] == i,2])
    }
    return(tmp)
}

# Used to generate prediction intervals for a given MCMC chain.
# Takes an MCMC chain (ie. the posterior distribution), a list of infection times,
# burnin period (optional), the corresponding data, number of samples from the
# posterior, the desired prediction level (default 95%), and level of smoothing of bounds
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
    




# Idea here is that when passed a list of lists of different lengths, returns the same list such that
# each list is the same size, with NAs used to pad the smaller lists. ie. all lists same size as
# biggest list
pad_lists <- function(my_list){
    max_length <- max(sapply(1:length(my_list),function(x) length(my_list[[x]])))
    for(i in 1:length(my_list)){
        my_list[[i]] <- c(my_list[[i]], rep(NA, max_length-length(my_list[[i]])))
    }
    return(my_list)
}

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

    

# Takes a list of parameters and a table of parameter bounds, transformation instructions etc. Idea is to convert all parameters to log/specified logit scale. 
transform_params_logitOLD<- function(params, param_table){
    tmp <- params

    for(i in 1:nrow(param_table)){
        if(param_table$use_log[i]) tmp[i] <- log(params[i])
        else if(param_table$use_logistic[i]) tmp[i] <- transform_logit(params[i], param_table$lower_bound[i], param_table$upper_bound[i])
        else tmp[i] <- params[i]
    }
    
    return(tmp)
}



# Takes a list of parameters and table of parameter bounds. Idea is to convert all parameters from log/specified logit scale to normal scale (-Inf to +Inf).
transform_params_logisticOLD<- function(params, param_table){
    tmp <- params

    for(i in 1:nrow(param_table)){
        if(param_table$use_log[i]) tmp[i] <- exp(params[i])
        else if(param_table$use_logistic[i]) tmp[i] <- transform_logistic(params[i], param_table$lower_bound[i], param_table$upper_bound[i])
        else tmp[i] <- params[i]
    }
    
    return(tmp)
}

