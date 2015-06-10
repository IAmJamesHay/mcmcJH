#' Calculates and sums priors for posterior
#'
#' Uses the prior function pointers provided by the parameter table to calculate prior values for the current parameters
#' @param params a vector of current parameter values
#' @param param_table parameter table as loaded by \code{\link{load_param_table}}. The idea is to provide prior function pointers and arguments for these, and pass them to a prior wrapper function
#' @return a single, sum of prior values
#' @export
#' @seealso \code{\link{prior_wrapper}}
prior <- function(params, param_table){
    sum <- sum(sapply(which(param_table$fixed==0), function(x) prior_wrapper(param_table$prior_func[[x]],params,x,param_table$prior_args[[x]])))
    return(sum)
}

#' Prior wrapper function
#' 
#' This is a wrapper function to return the value provided from FUNC when given the inputs "value" and "args"
#' @param FUNC a function pointer
#' @param values the current parameter values
#' @param the index of the parameter currently of interest
#' @param args arguments to be passed to the provided function
#' @return a single value for the prior
#' @export
#' @seealso \code{\link{prior}}
prior_wrapper <- function(FUNC, values, index,args){
    return(FUNC(values, index, args))
}

########################################################################
# THIS IS WHERE YOU SHOULD ADD ANY EXTRA PRIORS
# Functions should take only two arguments:
# 1. The current parameter value
# 2. An array of arguments to be given to the prior function (eg. mean, sd, bounds etc)
# The function MUST return a single numeric value, which will be summed (therefore take log priors)
########################################################################

# Normally distributed prior
mu_prior_norm <- function(values, index, args){
    return(dnorm(values[index],mean=args[1],sd=args[2],log=args[3]))
}

# Normally distributed prior
dRF_prior_norm <- function(values, index, args){
    return(dnorm(values[index],mean=args[1],sd=args[2],log=args[3]))
}

# Truncated normally distributed prior, with upper bound being Tp (ie. next index)
tp_prior_norm <- function(values, index, args){
    return(dnorm(values[index], mean=args[1], sd=args[2], log=args[3]))
}

tp_prior_norm2 <- function(values, index, args){
    return(log(dtruncnorm(values[index],mean=args[1],sd=args[2],b=values[index+1])))
}

# Truncated normally distributed prior, with lower bound being Ts (ie. previous index)
ts_prior_norm <- function(values, index, args){
    return(dnorm(values[index], mean=args[1], sd=args[2], log=args[3]))
}
ts_prior_norm2 <- function(values, index, args){
    return(log(dtruncnorm(values[index],mean=args[1],sd=args[2],a=values[index-1])))
}

# Normally distributed prior
m_prior_norm <- function(values, index, args){
    return(dnorm(values[index],mean=args[1],sd=args[2],log=args[3]))
}

sigma_prior_norm <- function(values, index, args){
    return(dnorm(values[index],mean=args[1],sd=args[2],log=args[3]))
}

no_prior <- function(){
    return(1)
}



#prior_table <- read.csv(prior_file,header=1)
#prior_table<- lapply(prior_table,as.character)
#prior_table$func <- lapply(prior_table$func, function(x) match.fun(x))
#prior_table$args <- lapply(prior_table$args, function(x) as.numeric(unlist(strsplit(x,","))))
