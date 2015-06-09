#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size. Requires \code{\link{toUnitScale}} and \code{\link{fromUnitScale}} from the mathfuncs package. Random walk may be on a linear or log scale
#' @param param a vector of the parameters to be explored
#' @param param_table a matrix of flags and bounds that specify how the proposal should be treated. Crucial elements are the 3rd column (lower bound), 4th column (upper bound), 5th column )step size) and 6th column (flag for log proposal)
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
proposalfunction <- function(param,param_table,index){
    # 4th index is upper bound, 3rd is lower
    # 1st and 2nd index used in other functions
    mn <- param_table[index,3]
    mx <- param_table[index,4]

    rtn <- param
    
    x <- rtn[index]
    
    x <- toUnitScale(param[index],mn,mx,param_table[index,6])

    # 5th index is step size
    stp <- param_table[index,5]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

    # Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

    # Cyclical boundary conditions
    #if (x < 0) x <- 1 + x	
    #if (x > 1) x <- x - 1
    
    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx,param_table[index,6])
    rtn
}

#' Posterior function
#'
#' Combines the likelihood and priors to produce a likelihood value for the given parameters against he provided data. This is used as the overall cost function for the MCMC algorithm.
#' @param params a vector of parameters to be explored
#' @param param_transform_table a matrix of flags and bounds as required by \code{\link{proposalfunction}} for the random walk proposal, for use in the likelihood function
#' @param param_table a data table as imported by \code{\link{load_param_table}} to be used by the prior wrapper function in finding the correct prior function and arguments
#' @param data the data of interest
#' @param LIKELIHOOD_FUNCTION a valid pointer to an R function which returns a single, log likelihood of the data given the current parameters
#' @param MODEL_FUNCTION a valid pointer to an R function which is used to evaluate the model for the current set of parameters
#' @return a single value for the posterior (by Bayes rule)
#' @export
#' @seealso \code{\link{prior}}
posterior <- function(params, param_transform_table, param_table, data, LIKELIHOOD_FUNCTION, MODEL_FUNCTION){
    return(LIKELIHOOD_FUNCTION(params,data,param_transform_table, 1, MODEL_FUNCTION) + prior(params, param_table))
}


#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{proposalfunction}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param startvalue a vector of parameter values used as the starting point for the MCMC chain. MUST be valid parameters for the model function
#' @param iterations number of iterations to run the MCMC chain for. Note that each parameter is moved once independently for each iteration. Defaults to 1000
#' @param data the data against which the likelihood is calculated
#' @param param_table a table of parameter data used for information such as bounds and prior function pointers. See \code{\link{load_param_table}}
#' @param popt the desired acceptance rate. Defaults to 0.44
#' @param opt_freq how frequently the acceptance rate is adapted. Defaults to 50
#' @param thin thinning value for the MCMC chain. Default is 1
#' @param burnin the length of the burn in period. Defaults to 100
#' @param adaptive_period length of the adaptive period. Defaults to 1
#'  @param LIKELIHOOD_FUNCTION a valid pointer to an R function which returns a single, log likelihood of the data given the current parameters
#' @param MODEL_FUNCTION a valid pointer to an R function which is used to evaluate the model for the current set of parameters
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param save_block the number of iterations that R will keep in memory before writing to .csv. Defaults to 500
#' @param VERBOSE boolean flag for additional output. Defaults to FALSE
#' @return the full file path at which the MCMC chain is saved as a .csv file
#' @export
#' @seealso \code{\link{posterior}}, \code{\link{proposalfunction}}
run_metropolis_MCMC <- function(startvalue,
                                iterations=1000,
                                data,
                                param_table,
                                popt=0.44,
                                opt_freq=50,
                                thin=1,
                                burnin=100,
                                adaptive_period=1,
                                LIKELIHOOD_FUNCTION,
                                MODEL_FUNCTION,
                                filename,
                                save_block = 500,
                                VERBOSE=FALSE
                                ){
    TUNING_ERROR<- 0.1

    if(opt_freq ==0 && VERBOSE){ print("Not running adaptive MCMC - opt_freq set to 0")}
    else if(VERBOSE){ print("Adaptive MCMC - will adapt step size during specified burnin period")}

    # Set up quicker tables
    # Turns the data frame table into a matrix that can allow faster indexing
    param_transform_table <- as.matrix(param_table[,c("use_logistic","use_log","lower_bound","upper_bound","step","log_proposal")])

    # Get those parameters which should be optimised
    non_fixed_params <- which(param_table$fixed==0)
    non_fixed_params_length <- length(non_fixed_params)
    all_param_length <- length(startvalue)

    # Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
    chain_colnames <- c("sampno",param_table$name,"lnlike")
  
    # Arrays to store acceptance rates
    tempaccepted <- tempiter <- reset <- integer(all_param_length)
    reset[] <- 0

    # Create empty chain to store "save_block" iterations at a time
    empty_chain <- chain <- matrix(nrow=save_block,ncol=all_param_length+2)
    
    # Set starting value and params
    current_params <- startvalue
    
    # Create array to store values
    empty_values <- values <- sample <- numeric(save_block)

    startvalue <- transform_params_logit(startvalue, param_transform_table)
    probab <- posterior(startvalue, param_transform_table, param_table, data, LIKELIHOOD_FUNCTION, MODEL_FUNCTION)
    startvalue <- transform_params_logistic(startvalue, param_transform_table)

    # Set up initial csv file
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,startvalue,probab)
    colnames(tmp_table) <- chain_colnames
    
    # Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    no_recorded <- 1
    sampno <- 2
    
    # Go through chain
    for (i in 1:(iterations+adaptive_period+burnin)){
        # For each parameter (Gibbs)
        for(j in non_fixed_params){
            # Propose new parameters and calculate posterior
            proposal <- transform_params_logit(proposalfunction(current_params,param_transform_table,j),param_transform_table)
            newprobab <- posterior(proposal, param_transform_table, param_table, data, LIKELIHOOD_FUNCTION, MODEL_FUNCTION)
            proposal <- transform_params_logistic(proposal, param_transform_table)

            # Calculate log difference in posteriors and accept/reject
            difflike <- newprobab - probab
         
            if ((!is.nan(difflike) & !is.infinite(newprobab)) & (runif(1) < exp(difflike) |  difflike > 0)){
                current_params <- proposal
                probab <- newprobab
                tempaccepted[j] <- tempaccepted[j] + 1
            }
            tempiter[j] <- tempiter[j] + 1

           # If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
           # save this block of chain to file and reset chain
            if(sampno %% thin ==0){
                chain[no_recorded,1] <- sampno
                chain[no_recorded,2:(ncol(chain)-1)] <- current_params
                chain[no_recorded,ncol(chain)] <- probab
                no_recorded <- no_recorded + 1
                
                if(no_recorded > save_block){
                    write.table(chain[1:(no_recorded-1),],file=mcmc_chain_file,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
                    chain <- empty_chain
                    no_recorded <- 1
                }
            }
            sampno <- sampno + 1
        }
        
        # Update step sizes based on acceptance rate
        # Note that if opt_freq is 0, then no tuning will take place
        if(opt_freq != 0 & i <= adaptive_period & i%%opt_freq== 0) {
            pcur <- tempaccepted/tempiter
            tempaccepted <- tempiter <- reset
            tmp_transform <- param_transform_table[,5]
            for(x in non_fixed_params){
                if(pcur[x] < popt - (TUNING_ERROR*popt) | pcur[x] > popt + (TUNING_ERROR*popt)){
                    tmp_transform[x] <- scaletuning(tmp_transform[x],popt,pcur[x])
                }
            }
            param_transform_table[,5] <- tmp_transform
        }
    }

    # If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    # that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    # rather than an array
    if(no_recorded > 2){
        write.table(chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }
    
    if(VERBOSE){
        print("Final step sizes:")
        print(param_table$step)
    }
    return(mcmc_chain_file)
}

#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param the current acceptance rate
#' @return the scaled step size
#' @export
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}

#' MCMC diagnostic tests
#'
#' Runs some basic MCMC diagnostics on the given chain and saves a few plots. The diagnostics are:
#' \itemize{
#' \item{Gelman Diagnostics: }{Saves the Gelman plots at the given file location}
#' \item{Auto-correlation: }{Saves autocorrelation plots at the given file location}
#' }
#' @param mcmc_chains the entire MCMC chain to be tested
#' @param filename the full file path at which to save the diagnostics. _gelman.pdf will be appended, for example
#' @param param_table the parameter table as loaded by \code{\link{load_param_table}}
#' @param VERBOSE boolean flag for additional output. Defaulst to FALSE
#' @return returns any error messages raised during the tests
#' @export
mcmc_diagnostics <- function(mcmc_chains, filename, param_table,VERBOSE=FALSE){
    errors <- NULL
    final <- NULL
    if(length(mcmc_chains) > 1){
        gelman.error <- tryCatch({
            if(VERBOSE) print("Saving Gelman diagnostics")
            gelman.filename <- paste(filename,"_gelman.pdf",sep="")
            pdf(gelman.filename)
            gelman.plot(as.mcmc.list(mcmc_chains)[,which(param_table$fixed==0)])
            if(VERBOSE) print("Gelman diagnostics:")
            gelman.diag(as.mcmc.list(mcmc_chains)[,which(param_table$fixed==0)])
            final <- NULL
        }, warnings=function(war1){
            if(VERBOSE) print(paste("A warnings occured in gelman diagnostics: ", war1))
            final <- war1
        }, error=function(err1){
            if(VERBOSE) print(paste("An error occured in gelman diagnostics: ", err1))
            final <- err1
        }, finally = {
            dev.off()
            errors <- c(errors, final)
        })
    }

    for(i in 1:length(mcmc_chains)){
        autocorr.error <- tryCatch({
            if(VERBOSE) print("Saving auto-correlation plot")
            autocorr.filename <- paste(filename,"_autocor_",i,".pdf",sep="")
            pdf(autocorr.filename)
            autocorr.plot(mcmc_chains[[i]][,which(param_table$fixed==0)])
            final <- NULL
        }, warnings=function(war2){
            if(VERBOSE) print(paste("A warnings occured in autocorr plot: ", war2))
            final <- war2
        }, error=function(err2){
            if(VERBOSE) print(paste("An error occured in autocorr plot diagnostics: ", err2))
            final <- err2
        }, finally = {
            dev.off()
            errors <- c(errors, final)            
        })
    }
    
    return(errors)
}
