#' Normally distributed errors likelihood function
#' 
#' Calculates the likelihood of a given vector of parameters for the provided data. Assumes error is normally distributed where standard deviation is the final argument of the parameter vector.
#' Restrictions on param_table: must be a data frame or matrix with headers:
#' value, lower_bound, upper_bound, use_logistic and use_log. Note that use_logistic and use_log should not both be TRUE.
#' Optimisation direction should be -1 or 1 depending on if we are maximising negative or positive log likelihood.
#' This should be -1 for optim, and +1 for MCMC (-1 for minimisation and +1 for maximisation).
#' @param params vector of parameters (ie. theta) for which likelihood is to be evaluated
#' @param data matrix or data frame of the data under evaluation. The first column should be a set of time points, whereas the second column is the values at these times
#' @param param_table optional table to accompany the parameters, specifying options such as log or logistic transformation
#' @param optimisation_direction numeric value taking +1 or -1 to specify whether returned value is a negative or positive log likelihood
#' @param MODEL_FUNCTION a function pointer for the model function that will be used to evaluated the model with the given set of parameters
#' @param ... any additional arguments that are required by the model function
#' @return returns a numeric value of the total log likelihood of the give set of parameters for the data
#' @export
#' @examples
#' params <- c(2,2)
#' data <- data.frame(t=seq(1,100,by=1),y=rnorm(100,1,1))
#' parameter_table <- NULL
#' optim_direction <- 1
#' model_function <- rnorm
#' nll <- llikelihood_norm(params,data,parameter_table,optim_direction,rnorm)
likelihood_norm <- function(params,data,param_table=NULL,optimisation_direction=-1, MODEL_FUNCTION, ...){
    # Convert parameters to normal space for prediction ONLY IF upper and lower bounds have been provided
    if(!is.null(param_table)){
        params <- transform_params_logistic(params,param_table)
    }
    out <- MODEL_FUNCTION(params, data[,1], ...)
    sigma <- params[length(params)]
    R <- data[,2] - out[,2]
    nll <- optimisation_direction*sum(dnorm(0,R,sigma,log=TRUE))
}


#' Poisson distributed errors likelihood function
#'
#' Calculates likelihood with given parameters for data. Assuming Poisson error distribution model.
#' Restrictions on param_table: must be a data frame or matrix with headers:
#' value, lower_bound, upper_bound, use_logistic and use_log. Note that use_logistic and use_log should not both be TRUE.
#' Optimisation direction should be -1 or 1 depending on if we are maximising negative or positive log likelihood.
#' This should be -1 for optim, and +1 for MCMC (-1 for minimisation and +1 for maximisation).
#' @param params vector of parameters (ie. theta) for which likelihood is to be evaluated
#' @param data matrix or data frame of the data under evaluation. The first column should be a set of time points, whereas the second column is the values at these times
#' @param param_table optional table to accompany the parameters, specifying options such as log or logistic transformation
#' @param optimisation_direction numeric value taking +1 or -1 to specify whether returned value is a negative or positive log likelihood
#' @param MODEL_FUNCTION a function pointer for the model function that will be used to evaluated the model with the given set of parameters
#' @param ... any additional arguments that are required by the model function
#' @return returns a numeric value of the total log likelihood of the give set of parameters for the data
#' @export
likelihood_poisson <- function(params,data,param_table=NULL,optimisation_direction=-1, MODEL_FUNCTION, ...) {
   # Convert parameters to normal space for prediction ONLY IF upper and lower bounds have been provided
    if(!is.null(param_table)){
        params <- transform_params_logistic(params,param_table)
    }
    out <- MODEL_FUNCTION(params, data[,1], ...)
    out[out<0] <- 0
    nll <- 0
    a <- out[,2]
    b <- data[,2]
    for(i in 1:nrow(out)){
        nll <- nll +  log(likelihood_simple(a[i],b[i]))
    }
    if(is.nan(nll) | is.infinite(nll)) nll <- -99999999999
    nll <- optimisation_direction*nll
}


#' Observation error only function
#'
#'
#' @export
likelihood_observation <- function(params, data, param_table=NULL,optimisation_direction=-1, MODEL_FUNCTION, ...) {
   # Convert parameters to normal space for prediction ONLY IF upper and lower bounds have been provided
    if(!is.null(param_table)){
        params <- transform_params_logistic(params,param_table)
    }
    out <- MODEL_FUNCTION(params, data[,1], ...)
    out[out<0] <- 0
    nll <- 0
    a <- out[,2]
    b <- data[,2]
    for(i in 1:nrow(out)){
        nll <- nll +  log(observation_model_kucharski(as.integer(a[i]),as.integer(b[i])))
    }
    if(is.nan(nll) | is.infinite(nll)) nll <- -99999999999
    nll <- optimisation_direction*nll
}
