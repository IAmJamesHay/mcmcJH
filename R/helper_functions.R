combine_numbers <- function(middle, lower, upper,sample_size){
  tmp <- NULL
  for(i in 1:length(middle)){
    tmp[i] <- paste(signif(middle[i],3), " (", signif(lower[i],3), "-",signif(upper[i],3),")",sep="")
    if(sample_size[i] != 0 & sample_size[i] < 200){
        tmp[i] <- paste(tmp[i], "*",sep="")
    }
  }
  return(tmp)
}

generate_artificial_data_big <- function(groups, strains, times, parameters, N_indivs){
    index <- 1
    all_data <- NULL
    for(group in groups){
        for(strain in strains){
            tmp <- generate_artificial_data(parameters[[index]],times,N_indivs)
            tmp$strain <- strain
            tmp$group <- group
            all_data <- rbind(all_data,tmp)
            index <- index + 1
        }
    }
    all_data <- all_data[,c(ncol(all_data),(ncol(all_data)-1),seq(1,(ncol(all_data)-2),by=1))]
    return(all_data)
}
generate_artificial_data <- function(parameters,times,individuals){
    pars <- parameters[3:length(parameters)]
    S <- parameters[1]
    EA <- parameters[2]


    y <- predict_titres(pars,times)
    y[y[,2] < 0,2] <- 0
    y[y[,2] > 13,2] <- 13
    
    data <- matrix(ncol=length(times)+1,nrow=individuals)
    data[,1] <- as.character(seq(1,individuals,by=1))
    
    for(i in 1:length(times)){
        tmp <- generate_sample(individuals,y[i,2],S,EA)
       
        data[,i+1] <- tmp
     }

    data <- as.data.frame(data)
    for(i in 2:ncol(data)){
        data[,i] <- as.numeric(as.character(data[,i]))
    }
    colnames(data) <- c("indiv",as.character(times))
    
    return(data)
}

generate_sample <- function(N,titre,S,EA){
    dist <- create_sample_distribution(titre,S,EA)
    results <- numeric(N)
    for(j in 1:N){
        a_number <- runif(1,0,1)
        p <- dist[1]
        i <- 1

        while(p < a_number){
            i <- i + 1
            p <- dist[i]
        }
        results[j] <- i-1
    }
    return(as.numeric(results))
}
        
        

create_sample_distribution <- function(titre, S, EA){
    test <- NULL
    titre <- floor(titre)
    test[1] <- obs_error(titre,0,S,EA)
    for(i in 1:13){
        test[i+1] <- obs_error(titre,i,S,EA) + test[i]
    }
    return(test)
}
    

consolidate_data <- function(param_files, start_dir){
    setwd(start_dir)
    cur_dir <- getwd()
    all_dirs <- list.files()
    not_dirs1 <- Sys.glob("*.pdf")
    not_dirs2 <- Sys.glob("*.png")
    not_dirs3 <- Sys.glob("*.csv")

    all_dirs <- setdiff(all_dirs, not_dirs1)
    all_dirs <- setdiff(all_dirs, not_dirs2)
    all_dirs <- setdiff(all_dirs, not_dirs3)
    
    print(all_dirs)
    column_names <- c("name","mean","sd","se","se2","X2.5","X25","X50","X75","X97.5","samp","mode","mle")

    tmp_table <- NULL
    tmp_table2 <- NULL

    largest_table <- NULL
    largest_nrows <- 0
    for(i in 1:length(param_files)){
        param_table <- load_param_table(param_files[[i]])
        if(nrow(param_table) > largest_nrows){
            largest_nrows <- nrow(param_table)
            largest_table <- param_table
        }
    }

    overall_names <- NULL
    for(i in 1:nrow(largest_table)){
        overall_names[i] <- paste(largest_table$name[i]," (",largest_table$lower_bound[i],"-",largest_table$upper_bound[i],")",sep="")
    }


    for(i in 1:length(all_dirs)){
        setwd(cur_dir)
        setwd(all_dirs[i])
        print(getwd())
        all_files <- Sys.glob("*.csv")
        #'print(param_files[[i]])
        param_table <- load_param_table(param_files[[i]])
        group <- all_dirs[i]
        
        for(j in 1:length(all_files)){
            strain <- all_files[j]
            #'print(all_files[j])
            tmp_summary <- read.csv(all_files[j])
            colnames(tmp_summary) <- column_names
            #'print(colnames(tmp_summary))
            names <- NULL
            values <- NULL
#'            print(tmp_summary$name)
            for(q in 1:(nrow(tmp_summary)-1)){
                name <- as.character(tmp_summary$name[q])

                mle <- tmp_summary$mle[q]
                mode <- tmp_summary$mode[q]
                mean <- tmp_summary$mean[q]
                lower_quantile <- tmp_summary$X2.5[q]
                upper_quantile <- tmp_summary$X97.5[q]

                lower_bound <- param_table$lower_bound[which(param_table$name == name)]
                upper_bound <- param_table$upper_bound[which(param_table$name == name)]
#'                print(upper_bound)
                samp <- tmp_summary$samp[q]

                combined <- as.character(combine_numbers(mean, lower_quantile, upper_quantile, samp))
                values[q] <- combined
                long_name <- paste(name, " (",as.character(lower_bound),"-",as.character(upper_bound),")",sep="")
                names[q] <- long_name
                
                row <- c(group, strain, name,mle,mode,mean,lower_quantile,upper_quantile,samp,lower_bound,upper_bound)
                tmp_table <- rbind(tmp_table,row)
            }
            names(values) <- names
            tmp_names <- rep(NA, length(overall_names))
            index <- 1
            for(f in 1:length(tmp_names)){
                if(overall_names[f] == names(values)[index]){
                    tmp_names[f] <- as.character(values)[index]
                    index <- index + 1
                }
            }
                    
            #'            tmp_names[which(overall_names==names(values))] <- as.character(values)
            #'print(tmp_names)
            tmp_table2 <- rbind(tmp_table2,c(group, strain,tmp_names))
        }
    }
    tmp_table <- as.data.frame(tmp_table)
    row.names(tmp_table) <- NULL
    colnames(tmp_table) <- c("group","strain","param","mle","mode","mean","lower_quantile","upper_quantile","samp_size","lower_bound","upper_bound")

    for(i in 1:ncol(tmp_table)){
        tmp_table[,i] <- as.character(tmp_table[,i])
    }


    tmp_table2 <- as.data.frame(tmp_table2)
    colnames(tmp_table2) <- c("Group","Strain",overall_names)
    row.names(tmp_table2) <- NULL


    for(i in 1:ncol(tmp_table2)){
        tmp_table2[,i] <- as.character(tmp_table2[,i])
    }
    tmp_table[tmp_table$group == "Grp 1 (H3N2)","group"] <- "Group 1"
    tmp_table[tmp_table$group == "Grp 2 (H3N2)","group"] <- "Group 2"
    tmp_table[tmp_table$group == "Grp 3 (H3N2)","group"] <- "Group 3"
    tmp_table[tmp_table$group == "Grp 4 (H3N2)","group"] <- "Group 4"
    tmp_table[tmp_table$group == "Grp 5 (H3N2)","group"] <- "Group 5"
    tmp_table[tmp_table$group == "Grp 1 (H1N1)","group"] <- "Group 1"
    tmp_table[tmp_table$group == "Grp 2 (H1N1)","group"] <- "Group 2"
    tmp_table[tmp_table$group == "Grp 3 (H1N1)","group"] <- "Group 3"
    tmp_table[tmp_table$group == "Grp 4 (H1N1)","group"] <- "Group 4"
    tmp_table[tmp_table$group == "Grp 5 (H1N1)","group"] <- "Group 5"
    tmp_table[tmp_table$group == "Grp 1 (B)","group"] <- "Group 1"
    tmp_table[tmp_table$group == "Grp 2 (B)","group"] <- "Group 2"
    tmp_table[tmp_table$group == "Grp 3 (B)","group"] <- "Group 3"
    tmp_table[tmp_table$group == "Grp 4 (B)","group"] <- "Group 4"
    tmp_table[tmp_table$group == "Grp 5 (B)","group"] <- "Group 5"

    a <- intersect(grep("H3N2",tmp_table$strain),grep("panama",tmp_table$strain))
    tmp_table[a,"strain"] <- "A/Panama/2007/1999 (H3N2)"

    a <- intersect(grep("H3N2",tmp_table$strain),grep("brisbane",tmp_table$strain))
    tmp_table[a,"strain"] <- "A/Brisbane/10/2007 (H3N2)"

    a <-  intersect(grep("H3N2",tmp_table$strain),grep("wisconsin",tmp_table$strain))
    tmp_table[a,"strain"] <- "A/Wisconsin/67/2005 (H3N2)"

    a <- intersect(grep("H1N1",tmp_table$strain),grep("fukushima",tmp_table$strain))
    tmp_table[a,"strain"] <- "A/Fukushima/141/06 (H1N1)"

    a <- intersect(grep("H1N1",tmp_table$strain),grep("solomon",tmp_table$strain))
    tmp_table[a,"strain"] <- "A/Solomon Islands/3/2006 (H1N1)"
    
    a <- intersect(grep("B",tmp_table$strain),grep("brisbane",tmp_table$strain))
    tmp_table[a,"strain"] <- "B/Brisbane/3/2007 (B)"
    
    a <- intersect(grep("B",tmp_table$strain),grep("malaysia",tmp_table$strain))
    tmp_table[a,"strain"] <- "B/Malaysia/2506/2004 (B)"


    tmp_table2[tmp_table2$Group == "Grp 1 (H3N2)","Group"] <- "Group 1"
    tmp_table2[tmp_table2$Group == "Grp 2 (H3N2)","Group"] <- "Group 2"
    tmp_table2[tmp_table2$Group == "Grp 3 (H3N2)","Group"] <- "Group 3"
    tmp_table2[tmp_table2$Group == "Grp 4 (H3N2)","Group"] <- "Group 4"
    tmp_table2[tmp_table2$Group == "Grp 5 (H3N2)","Group"] <- "Group 5"
    tmp_table2[tmp_table2$Group == "Grp 1 (H1N1)","Group"] <- "Group 1"
    tmp_table2[tmp_table2$Group == "Grp 2 (H1N1)","Group"] <- "Group 2"
    tmp_table2[tmp_table2$Group == "Grp 3 (H1N1)","Group"] <- "Group 3"
    tmp_table2[tmp_table2$Group == "Grp 4 (H1N1)","Group"] <- "Group 4"
    tmp_table2[tmp_table2$Group == "Grp 5 (H1N1)","Group"] <- "Group 5"
    tmp_table2[tmp_table2$Group == "Grp 1 (B)","Group"] <- "Group 1"
    tmp_table2[tmp_table2$Group == "Grp 2 (B)","Group"] <- "Group 2"
    tmp_table2[tmp_table2$Group == "Grp 3 (B)","Group"] <- "Group 3"
    tmp_table2[tmp_table2$Group == "Grp 4 (B)","Group"] <- "Group 4"
    tmp_table2[tmp_table2$Group == "Grp 5 (B)","Group"] <- "Group 5"

    a <- intersect(grep("H3N2",tmp_table2$Strain),grep("panama",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "A/Panama/2007/1999 (H3N2)"
    
    a <- intersect(grep("H3N2",tmp_table2$Strain),grep("brisbane",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "A/Brisbane/10/2007 (H3N2)"

    a <- intersect(grep("H3N2",tmp_table2$Strain),grep("wisconsin",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "A/Wisconsin/67/2005 (H3N2)"

    a <- intersect(grep("H1N1",tmp_table2$Strain),grep("fukushima",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "A/Fukushima/141/06 (H1N1)"

    a <- intersect(grep("H1N1",tmp_table2$Strain),grep("solomon",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "A/Solomon Islands/3/2006 (H1N1)"
     
    a <- intersect(grep("B",tmp_table2$Strain),grep("brisbane",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "B/Brisbane/3/2007 (B)"
    
    a <- intersect(grep("B",tmp_table2$Strain),grep("malaysia",tmp_table2$Strain))
    tmp_table2[a,"Strain"] <- "B/Malaysia/2506/2004 (B)"

    #'a <- grep("B",tmp_table2$Strain)[grep("brisbane",tmp_table2$Strain)]
    #'a <- a[!is.na(a)]
    #'tmp_table2[a,"Strain"] <- "Brisbane (B)"

    #'a <- grep("B",tmp_table2$Strain)[grep("malaysia",tmp_table2$Strain)]
    #'a <- a[!is.na(a)]
   #' tmp_table2[a,"Strain"] <- "Malaysia (B)"


    setwd(cur_dir)
    write.table(tmp_table,"all_results_shiny.csv",sep=",",row.names=FALSE,col.names=TRUE,append=FALSE)
    write.table(tmp_table2,"all_results_shiny_table.csv",sep=",",row.names=FALSE,col.names=TRUE,append=FALSE)
}



#' Takes a list of filenames, and reads these in, returning an MCMC chain
read_chain_list <- function(filenames, col_names, disregard, thin=1, indices= NULL){
    tmp_chains <- NULL
    for(i in 1:length(filenames)){
        tmp <- read.csv(filenames[[i]],header=0)
        colnames(tmp) <- col_names
        tmp <- tmp[tmp[,1] >= disregard,]
         if(!is.null(indices)){
            tmp <- tmp[,indices]
        }
        tmp <- tmp[seq(1,nrow(tmp),by=thin),]
        tmp_chains[[i]] <- as.mcmc(tmp)
    }
    return(tmp_chains)
}



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
    param_table <- read.csv(param_file,header=1)
    # Make some slight formatting adjustments to data and parameter tables
    # Parameters
    colnames(param_table) <- c("names","value","lower_bound","upper_bound","use_logistic","use_log","step","fixed","nuisance","lower_seed","upper_seed","prior_func","prior_args","log_proposal","block")
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
#' @param level desired prediction interval size. Defaults to 0.95 ie. 95 percent prediction intervals
#' @param smoothing uses smoothing splines to give smoother prediction intervals where these might be aethsetically spikey. Defaults to 0
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
        mcmc <- rbind(mcmc, mcmc_chains[[i]])
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
#' \item{model_data: }{a data frame with the original data (ie. fitted) for use by a plotting function}
#' \item{prediction_data: }{a data frame with the best fitting model's data for use by a plotting function}
#' \item{lower_prediction_bounds: }{a data frame with the necessary information required to plot the lower prediction bound}
#' \item{upper_prediction_bounds: }{a data frame with the necessary information required to plot the upper prediction bound}
#' }
#' @seealso \code{\link{generate_prediction_intervals_1.2}}
#' @export
generate_prediction_intervals_1.1 <- function(results_MCMC, burnin, param_table, all_data,MODEL_FUNCTION,group){
    upper_bounds <- NULL
    lower_bounds <- NULL
    predictions <- NULL
    index <- 1
    #' Generate prediction intervals for each fitted strain
    for(j in 1:length(results_MCMC)){
        if(is.null(results_MCMC[[j]]) | length(results_MCMC[[j]]) < 1){
            next
        }
        
        tmp_chains <- read_chain_list(
            results_MCMC[[j]]$files,
            c("sampno", param_table$names,"lnlike"),
            burnin,
            1,
            (which(param_table$nuisance==0)+1)
            )
        
        #' Generate model predictions for all fits
        #' Generate prediction intervals
        bounds <- generate_prediction_intervals_1.2(tmp_chains,
                                                    NULL,
                                                    results_MCMC[[j]]$times,
                                                    runs=10000,
                                                    level=0.95,
                                                    smoothing=0.01,
                                                    MODEL_FUNCTION
                                                    )
        upper_bounds[[index]] <- bounds$upper
        lower_bounds[[index]] <- bounds$lower
        predictions[[index]] <- as.data.frame(
            MODEL_FUNCTION(
                results_MCMC[[j]]$params[which(param_table$nuisance ==0)],
                seq(results_MCMC[[j]]$times[1],max(results_MCMC[[j]]$times),by=0.1)
                )
            )
        colnames(predictions[[index]]) <- c("variable","value")
        lower_bounds[[index]]$strain <- upper_bounds[[index]]$strain <- predictions[[index]]$strain <- unique(all_data$strain)[j]
        lower_bounds[[index]]$group <- upper_bounds[[index]]$group <- predictions[[index]]$group <- group
        index <- index + 1
    }
    all_model_dat <- NULL
    for(j in 1:length(predictions)){
        if(length(predictions) >0 & nrow(predictions[[j]]) > 0){
            all_model_dat <- rbind(all_model_dat, predictions[[j]])
        }
    }
    return(list("model_data"=all_model_dat,"prediction_data"=all_data,"lower_prediction_bounds"=lower_bounds,"upper_prediction_bounds"=upper_bounds))
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

