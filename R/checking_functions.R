#' Model parameter table check
#' 
#' Takes the table of parameters to be used for model calculation and MCMC fitting, and returns a string error if any of the inputs are invalid. The table should have the following column names:
#' \itemize{
#' \item{"name"}{Name of the parameter, string or character}
#' \item{"value"}{arbritary numeric value that represents a typical value. This is not that important, it is just used for some upkeep functions}
#' \item{"lower_bound"}{numeric value for the lower allowable bound of the parameter. This can be -Inf}
#' \item{"upper_bound"}{numeric value for the upper allowable bound of the parameter. This can be Inf}
#' \item{"use_logistic"}{boolean value (1/0) indicating whether or not the parameter should be logistic transformed during optimisation, ensuring that the parameter search space is between the specified bounds}
#' \item{"use_log"}{boolean value (1/0) indicating whether parameter optimisation should be on a log scale. Note that this and use logistic should NOT both be specified}
#' \item{"step"}{initial step size for the MCMC random walk}
#' \item{"fixed"}{boolean value (1/0) indicating whether parameter should be fixed during MCMC/optimisation or not}
#' \item{"nuisance"}{boolean value (1/0) indicating if parameter is a nuisance parameter}
#' \item{"lower_seed"}{numeric value for the lower bound of the random seed generation}
#' \item{"upper_seed"}{numeric value for the upper bound of the random seed generation}
#' \item{"prior_func"}{this should be a string that matches the name of a function that will be used as the prior function for MCMC. This function should be available in the current R environment, or an error will be thrown. See below.}
#' \item{"prior_args"}{a vector of default arguments for the prior function. For example, this might contain the standard deviation or mean of a normal distribution. Most prior functions have a flag to indicate log probabilities, which should be set to true by this argument}
#' \item{"log_proposal"}{boolean flag to indicate whether proposals in the MCMC function should be made on a log or linear scale. 1 for log, 0 for linear. If in doubt, set this to 0.}
#' }
#' 
#' @param param_table a data frame or matrix containing the parameters and corresponding meta-data.
#' @return returns a string containing an error message, or NULL if no error
#' @seealso \code{\link{mcmc_param_check}}, \code{\link{prior}}, \code{\link{prior_wrapper}}
#' @export
param_check <- function(param_table){
    if(ncol(param_table) != 14){
        return("Error in parameters - incorrect number of parameter table columns. Should be value, lower_bound, upper_bound, use_logistic, use_log, step")
    }
    correct_length <- max(sapply(1:ncol(param_table), function(x) length(param_table[,x])))
    for(i in 1:ncol(param_table)){
        if(length(param_table[,i]) < correct_length){
            print(length(param_table[,i]))
            return("Error in parameters - inputs not provided for all columns!")
        }
    }
    correct_height <- max(sapply(1:nrow(param_table),function(x) length(param_table[x,])))
    for(i in 1:nrow(param_table)){
        if(length(param_table[i,]) < correct_height){
            return("Error in parameters - inputs not provided for all rows!")
        }
    }
    for(i in 1:nrow(param_table)){
        if(!is.numeric(param_table[i,2])){
            return("Error in parameters - invalid value provided")
        }
        if(!is.numeric(param_table[i,3]) | !is.numeric(param_table[i,4]) | (param_table[i,3] > param_table[i,4])){
            return("Error in parameters - invalid bounds provided")
        }
        if((param_table[i,5] != 0 & param_table[i,5] != 1) | (param_table[i,6] != 0 & param_table[i,6] != 1) | (param_table[i,14] !=0 & param_table[i,14] != 1)){
            return("Error in parameters - invalid log/logistic flags provided")
        }
        if(!is.numeric(param_table[i,7]) | param_table[i,7] > 1){
            return("Error in parameters - invalid step size provided")
        }
        if((param_table[i,8] != 0 & param_table[i,8] != 1) | (param_table[i,9] != 0 & param_table[i,9] != 1)){
            return("Error in parameters - have not specified which parameters to fix in optimisation")
        }
        if(!is.numeric(param_table[i,10]) | !is.numeric(param_table[i,11]) | (param_table[i,10] > param_table[i,11])){
            return("Error in parameters - invalid seeding bounds provided")
        }
    }
    return(NULL)
}
#' MCMC parameter table check
#'
#' Takes the table of parameters for MCMC and checks for any invalid inputs. Table should have the following columns:
#' \itemize{
#' \item{iterations}{numeric value for the number of iterations in the MCMC algorithm}
#' \item{opt_freq}{numeric value, specifying that the step sizes should be adapted every opt_freq iterations}
#' \item{thin}{numeric thinning value for the MCMC chain}
#' \item{burnin}{numeric value specifying the number of iterations to be discarded as burnin. Note that this is in addition to the adaptive period}
#' \item{adaptive_period}{numeric value specifying the number of iterations to be discarded as the adpative period (ie. adapting step size)}
#' \item{nchain}{number of chains to be run}
#' \item{save_results}{currently does not have a use, but would take a value 1 or 0 to represent boolean value}
#' \item{likelihood_function}{this should be a string that matches the name of a function available in the R environment. See the cost_functions documentation for details. Basically, this should be a function that takes the vector of parameters to be evaluated, the data to be evaluated against, the parameter table, the optimisation direction (ie. negative or positive log likelihood), and a function pointer for the model function to be evaluated}
#' \item{model_function}{this should be a strnig that matches the name of a function in the R environment. This function will take a vector of parameters and a vector of time points, and will return the evaluated model results}
#' \item{lower_plot_bound}{numeric value specifying the lower y bound for the model plot}
#' \item{upper_plot_bound}{numeric value specifying the upper y bound for the model plot}
#'}
#' @param mcmc_params a data frame or matrix containing the parameters and corresponding meta-data
#' @return returns a string containing an error message, or NULL if no error
#' @seealso \code{\link{param_check}}
#' @export
mcmc_param_check <- function(mcmc_params){
    if(ncol(mcmc_params) != 12){
        return("Error in MCMC parameters - incorrect number of mcmc inputs provided")
    }
   # if(mcmc_params$iterations < mcmc_params$burnin){
     #   return("Error in MCMC parameters - burnin period longer than chain length")
    #}
    if(mcmc_params$nchain < 1){
        return("Error in MCMC parameters - 0 chains specified")
    }
    if(mcmc_params$thin < 1){
        return("Error in MCMC parameters - thinning less than 1")
    }
    if(mcmc_params$popt < 0 | mcmc_params$popt > 1){
        return("Error in MCMC parameters - desired acceptance rate invalid")
    }
    if(!is.numeric(mcmc_params$lower_plot_bound) | !is.numeric(mcmc_params$upper_plot_bound) | (mcmc_params$upper_plot_bound <= mcmc_params$lower_plot_bound)){
        return("Error in MCMC parameters - plot bounds invalid")
    }
    
    final <- NULL
    error <- tryCatch({
        func <- match.fun(as.character(mcmc_params$likelihood_function))
        final <- NULL
    }, warning=function(war){
                final <- (paste("A warning occured in MCMC loading: ", war))
    },error=function(err){
        final <- (paste("Error in MCMC parameters - ",err))
    }, finally = {
    })

    final <- NULL
    error <- tryCatch({
        func <- match.fun(as.character(mcmc_params$model_function))
        final <- NULL
    }, warning=function(war){
        final <- (paste("A warning occured in MCMC loading: ", war))
    },error=function(err){
        final <- (paste("Error in MCMC parameters - ",err))
    }, finally = {
    })
    
    return(final)
}

# Checks that provided data is all valid
data_check <- function(raw_data){
   if(colnames(raw_data)[1] != "group" | colnames(raw_data)[2] != "strain" | colnames(raw_data)[3] != "indiv"){
        return("Error in data - invalid column names!")
    }
    for(i in 4:ncol(raw_data)){
        if(is.na(as.numeric(colnames(raw_data)[i]))){
            return("Error in data - invalid time point!")
        }
    }
    for(i in 1:nrow(raw_data)){
        for(j in 4:ncol(raw_data)){
            if(!is.numeric(raw_data[i,j])){
                return("Error in data - invalid titre value given")
            }
        }
    }
    return(NULL)
}

MCMC_main_checks <- function(mcmc_params, param_table, raw_data,infection_table,prior_table){
    ERROR_MCMC<- mcmc_param_check(mcmc_params)
    ERROR_PARAM<- param_check(param_table)
    ERROR_DATA<- data_check(raw_data)
    
    STOP <- FALSE
    if(!is.null(ERROR_MCMC)){
        print("Error in MCMC parameter file")
        print(ERROR_PARAM)
        STOP <- TRUE
    }
    if(!is.null(ERROR_PARAM)){
        print("Error in parameter file")
        print(ERROR_PARAM)
        STOP <- TRUE
    }
    if(!is.null(ERROR_DATA)){
        print("Error in provided data file")
        print(ERROR_DATA)
        STOP <- TRUE
    }
    
    if(STOP) return(1)
    return(0)
}


