#' Main MCMC wrapper function
#' 
#' This is the main wrapper function for the ferret model data MCMC fits. This will load any necessary data or parameter files and then call the next level MCMC wrapper function, \code{\link{MCMC_fit_1.1}}. Brief description:
#' Firstly, reads in the parameter files and carries out some checks on these. eg. numeric values where needed, bounds being the right way around etc. Then formats these parameter tables correctly and matches function names with the correct functions. Also melts data to correct format. Creates a set of file names that are used for each sub-group to be fitted (ie. strain). Fairly specific to this problem, also creates an array of infection times
#' @param data_file the correct location of a .csv file containing the data to be fitted. In this case, should have four columns: group, strain, individual, titre and time
#' @param  param_file the correct location of a .csv file containing information about the parameters that describe the model of interest. See \code{\link{load_param_table}} for further details
#' @param mcmc_param_file the correct location of a .csv file containing the parameters for the MCMC fitting. See \code{\link{mcmc_param_check}} for further details
#' @param output_file (OPTIONAL) if specified, all standard output will be parsed to this specified file. NOT WELL TESTED, SO RECOMMEND AVOIDING
#' @param main_dir top level directory in which to work. Defaults to the current directory
#' @param group this is fairly specific to this data. Just specifies which group from the data file should be explored
#' @param VERBOSE (OPTIONAL) - if true, will give printed outputs from the program
#' @param infection_labels a vector of times and labels to represent infection events. These are used by the model plotting function, and are optional                                        
#' @return 1 if an error, 0 otherwise.
#' @seealso \code{\link{MCMC_fit_1.1}}
#' @export               
MCMC_main <- function(
    data_file,
    param_file,
    mcmc_param_file,
    output_file=NULL,
    main_dir="",
    group,
    VERBOSE=FALSE,
    infection_labels=NULL
    ){

    # Read in all data and parameters
    # MCMC parameters
    mcmc_params <- read.csv(paste(getwd(),mcmc_param_file,sep=""),header=1)
    
    # Parameters for model
    param_table <- read.csv(paste(getwd(),param_file,sep=""),header=1)

    # Data to be fitted
    raw_data <- read.csv(paste(getwd(),data_file,sep=""),stringsAsFactors=FALSE,na.strings=c("NA","-","?"),header=0)
    colnames(raw_data) <- raw_data[1,]
    raw_data <- raw_data[2:nrow(raw_data),]
    
    # Check that all imported tables fulfill necessary checks. If not, return early.
    ERROR <- MCMC_main_checks(mcmc_params,param_table,raw_data)

    # End if there was an error
    if(ERROR == 1){
        return(1)
    }

    param_table <- load_param_table(param_file)
    
    # Data
    melted.data <- melt(raw_data)
    melted.data$variable <- as.integer(as.character(melted.data$variable))
    melted.data <- melted.data[melted.data$group == group,]
    
    # Save current WD and then move to specified WD
    oldwd <- getwd()
    dir.create(main_dir)
    setwd(main_dir)

    # Match up correct likelihood function
    mcmc_params$likelihood_function <- match.fun(as.character(mcmc_params$likelihood_function))
    mcmc_params$model_function <- match.fun(as.character(mcmc_params$model_function))
    # Create filenames from the provided data
    file_names <- tolower(unlist(lapply(unique(melted.data$strain),function(x) unlist(strsplit(x,"/"))[2])))
    file_names <- gsub("[[:punct:]]|\\s+","_",file_names)

    # Find infection times from params
    infection_times <- infection_labels[infection_labels$group==group,c("time","infection")]
    
    ##################
    # RUN MAIN MCMC SCRIPT
    ##################
    final <- MCMC_fit_1.1(main_dir,
                          melted.data,
                          output_file,
                          file_names,
                          param_table,
                          mcmc_params$iterations,
                          mcmc_params$opt_freq,
                          mcmc_params$thin,
                          mcmc_params$burnin,
                          mcmc_params$adaptive_period,
                          mcmc_params$nchain,
                          mcmc_params$popt,
                          group,
                          mcmc_params$lower_plot_bound,
                          mcmc_params$upper_plot_bound,
                          mcmc_params$likelihood_function,
                          mcmc_params$model_function,
                          VERBOSE,
                          infection_times
                          )

    # Return to old WD and finish
    setwd(oldwd)
    return(final)
}

#' First level wrapper function for the MCMC algorithm.
#' 
#' Function to provide MCMC (Adaptive Metropolis within Gibbs) fits to all data subsets from a given data frame. No return value - all relevant outputs are saved as .csv files. Brief desription:
#' 
#' Firstly, sets up number of cores to be used. If only one chain, then parallelises the strains. If more than one chain, puts each chain on a different core instead. Checks for any all-zero data, and removes these from further analysis. Goes through each group in the data frame and creates a new directory for each. Gets the number of strains to be fitted and goes through each one (parallel or sequentially depending on number of cores)
#' @param top_dir  location of the top level directory to work in
#' @param all_data  data frame of the data to be fit. Should be of the form "group", "strain", "variable", "value", where variable is the numeric time point and value is the titre level
#' @param output_file (OPTIONAL) - should be a text file to output model fitting printouts
#' @param param_table the data.frame version of the parameter table read in from the .csv file above. See \code{\link{load_param_table}}
#' @param iterations number of iterations for the MCMC algorithm
#' @param opt_freq how frequently the acceptance rate is adapted. Defaults to 50
#' @param thin thinning value for the MCMC chain. Default is 1
#' @param burnin the length of the burn in period. Defaults to 100
#' @param adaptive_period length of the adaptive period. Defaults to 1
#' @param nchain number of chains to run
#' @param popt the desired acceptance rate. Defaults to 0.44
#' @param group the name of the current group under investigation. Mostly just used to create a directory for the group specific results
#' @param lower_plot_bound lowest value to be plotted ont he y axis for model trajectory. Defaults to 0
#' @param upper_plot_bound highest value to be plotted ont he y axis for model trajectory. Defaults to 15
#' @param LIKELIHOOD_FUNCTION a valid pointer to an R function which returns a single, log likelihood of the data given the current parameters
#' @param MODEL_FUNCTION a valid pointer to an R function which is used to evaluate the model for the current set of parameters
#' @param VERBOSE flag for additional output. Defaults to FALSE
#' @param infection_times an optional vector of infection event labels and times to be plotted on the model graph
#' @param SKIP_ZEROES optional flag indicating whether or not all zero value groups should be skipped. Defaults to FALSE
#' @return returns nothing!
#' @export
#' @seealso \code{\link{MCMC_fit_1.2}}
MCMC_fit_1.1 <- function(top_dir,
                         all_data,
                         output_file=NULL,
                         file_strains,
                         param_table,
                         iterations=1000,
                         opt_freq=50,
                         thin=50,
                         burnin=100,
                         adaptive_period=1,
                         nchain=1,
                         popt=0.44,
                         group=NULL,
                         lower_plot_bound=0,
                         upper_plot_bound=15,
                         LIKELIHOOD_FUNCTION=NULL,
                         MODEL_FUNCTION=NULL,
                         VERBOSE=FALSE,
                         infection_times=NULL,
                         SKIP_ZEROES=FALSE
                         ){
    
    no_cores<- detectCores()

    if(.Platform$OS.type=="unix"){
        registerDoParallel(cores=(no_cores-1))
    } else {
        c1 <- makeCluster(no_cores-1)
        registerDoParallel(c1)
    }
    
    # Pipe output to given file
    if(!is.null(output_file)) sink(output_file,append=TRUE)

    # Check given data. If any all zeroes, remove from further analysis
    if(SKIP_ZEROES){
        tmp <- remove_zero_data(all_data,file_strains)
        all_data <- tmp[[1]]
        file_strains <- tmp[[2]]
    }

    if(VERBOSE){
        print("==================================================")
        print(" COMMENCING PERSONAL MCMC FITS")
        print("==================================================")
    }

    # Get in top level working directory, and then create/move into this group's directory
    setwd(top_dir)
    dir.create(group,showWarnings=FALSE)
    results_MCMC <- NULL
   # We will then fit the model to each strain
    number_strains <- unique(all_data$strain)
    
    results_MCMC <- foreach(i=1:length(number_strains),.packages='mcmcJH') %dopar%{
        setwd(group)
        # Get subset of data
        temp_dat <- all_data[all_data$group == group & all_data$strain == number_strains[i], c("variable","value")]
        
                                        # Just to show what is going on
        if(VERBOSE){
            print("")
            print("-----------------------------------------------")
            print(paste("MCMC fit for ", group, ", ", number_strains[i]," underway:",sep=""))
        }
        
        tmp_filename <- paste(group,"_",file_strains[i],sep="")
        
        # Returns a list of MLE parameters, corresponding data time points (for plots) and a list of files containing the MCMC chains
        final <- MCMC_fit_1.2(temp_dat,
                              param_table,
                              iterations,
                              opt_freq,
                              thin,
                              burnin,
                              adaptive_period,
                              nchain,
                              popt,
                              tmp_filename,
                              LIKELIHOOD_FUNCTION,
                              MODEL_FUNCTION,
                              VERBOSE,
                              topdir=paste(top_dir,"/tmp",sep="")
                              )
    }
    # Fairly elaborate way of generation 95% prediction intervals. Probably room for improvement...
    prediction_intervals<- generate_prediction_intervals_1.1(results_MCMC, (burnin+adaptive_period), param_table, all_data, MODEL_FUNCTION)
    plot_model_fits(prediction_intervals$model_data,
                    prediction_intervals$prediction_data,
                    group,
                    c(lower_plot_bound,upper_plot_bound),
                    prediction_intervals$lower_prediction_bounds,
                    prediction_intervals$upper_prediction_bounds,
                    infection_times=infection_times)

    if(!is.null(output_file)) sink()
    
    if(.Platform$OS.type=="unix"){
        registerDoParallel(cores=1)
    } else {
        stopCluster(c1)
    }
    return(list(prediction_intervals$model_data,
                    prediction_intervals$prediction_data,
                    group,
                    c(lower_plot_bound,upper_plot_bound),
                    prediction_intervals$lower_prediction_bounds,
                    prediction_intervals$upper_prediction_bounds,
                    infection_times))
}

#' Second level wrapper for the MCMC fitting procedure. HIGHEST, GENERALISED LEVEL
#'
#' This is the highest level wrapper function for the MCMC algorithm that can be applied to any generic model, data and likelihood function. Calls another function, \code{\link{MCMC_fit_single}} to carry out the random walk. This function takes care of formatting the results and carrying out some MCMC diagnostics and plots.
#' @param temp_dat the data over which to calculate likelihoods
#' @param param_table table of parameter data as specified in \code{\link{load_param_table}}
#' @param iterations number of iterations for the MCMC algorithm
#' @param opt_freq how frequently the acceptance rate is adapted. Defaults to 50
#' @param thin thinning value for the MCMC chain. Default is 1
#' @param burnin the length of the burn in period. Defaults to 100
#' @param adaptive_period length of the adaptive period. Defaults to 1
#' @param nchain number of chains to run
#' @param popt the desired acceptance rate. Defaults to 0.44
#' @param filename the generic file name/location at which to save plots and diagnostics. Note that the file extension will be appended to this string
#' @param LIKELIHOOD_FUNCTION a valid pointer to an R function which returns a single, log likelihood of the data given the current parameters
#' @param MODEL_FUNCTION a valid pointer to an R function which is used to evaluate the model for the current set of parameters
#' @param VERBOSE boolean flag for additional output. Defaults to FALSE
#' @param PARALLEL OPTIONAL boolean flag indicating if each MCMC chain should be run in parallel using doPar. Defaults to FALSE
#' @param OPTIM_PARAMS OPTIONAL boolean flag indicating if optim should be used to find the "best fit" parameters. Otherwise the maximum likelihood parameters from the MCMC chain will be used. Defaults to FALSE
#' @param topdir the top level directory where the MCMC chains will be saved. By default, a tmp folder will be created in the home folder
#' @return returns a list containing a vector of the maximum likelihood parameters, the time points of the fitted data, and a list of the MCMC chain files
#' @export
#' @seealso \code{\link{MCMC_fit_single}}
MCMC_fit_1.2 <- function(temp_dat,
                         param_table,
                         iterations=1000,
                         opt_freq=50,
                         thin=1,
                         burnin=100,
                         adaptive_period=1,
                         nchain=1,
                         popt=0.44,
                         filename,
                         LIKELIHOOD_FUNCTION,
                         MODEL_FUNCTION,
                         VERBOSE=FALSE,
                         PARALLEL=FALSE,
                         OPTIM_PARAMS=FALSE,
                         topdir=paste(getwd(),"/tmp/",sep="")
){
    # Set up quicker tables
    # Turns the data frame table into a matrix that can allow faster indexing
    param_transform_table <- as.matrix(param_table[,c("use_logistic","use_log","lower_bound","upper_bound","step")])


    if(!is.null(topdir)){
        oldwd <- getwd()
        dir.create(topdir)
        setwd(topdir)
    }
    
    tmp_mcmc_results <- MCMC_fit_single(temp_dat,
                                        param_table,
                                        iterations,
                                        opt_freq,
                                        thin,
                                        burnin,
                                        adaptive_period,
                                        nchain,
                                        popt,
                                        filename,
                                        LIKELIHOOD_FUNCTION,
                                        MODEL_FUNCTION,
                                        VERBOSE,
                                        PARALLEL
                                        )
    
    if(VERBOSE){
        print("MCMC fit successful! Summary:")
        print(tmp_mcmc_results$summary)
    }
    
  # Get effective sample size
    effective_sample_size <- effectiveSize(as.mcmc.list(tmp_mcmc_results$chains))
    if(VERBOSE){
        print("Effective sample size:")
        print(effective_sample_size)
    }
  
  # Do optim to find maximum likelihood estimate
    if(OPTIM_PARAMS){
        optim_start <- tmp_mcmc_results$summary$statistics[1:nrow(param_table),1]
        optim_start <- transform_params_logit(optim_start,param_transform_table)
        
        if(VERBOSE){
            print("Fixed for optim: ")
            print(param_table$names[c(bitwOr(param_table$nuisance, param_table$fixed)==1)])
        }
        
        tmp_optim_result <- optifix(par=optim_start,
                                    fixed=c(bitwOr(param_table$nuisance, param_table$fixed)==1),
                                    fn=LIKELIHOOD_FUNCTION,
                                    data=temp_dat,
                                    param_table=param_transform_table,
                                    optimisation_direction=-1,
                                    MODEL_FUNCTION=MODEL_FUNCTION,
                                    lower=transform_params_logit(param_table$lower_bound,param_transform_table),
                                    upper=transform_params_logit(param_table$upper_bound,param_transform_table),
                                    method="L-BFGS-B",
                                    control=list(maxit=10000)
                                    )
        if(VERBOSE){
            print("Results from OPTIM:")
            print(tmp_optim_result[c("fullpars","value","counts")])
        }
        
        best_pars <- transform_params_logistic(tmp_optim_result$fullpars,param_transform_table)
        mode_pars <- best_pars
    }
    else {
        best_pars_plot <- tmp_mcmc_results$params
        best_pars <- as.numeric(tmp_mcmc_results$params)
        mode_pars <- tmp_mcmc_results$mode_params
    }

  # Save MCMC summary to the appropriate csv file
    tmp_table <- tmp_mcmc_results$summary$statistics
    tmp_table <- cbind(tmp_table, tmp_mcmc_results$summary$quantiles,effective_sample_size,c(mode_pars,NA),c(best_pars,NA))
    colnames(tmp_table)[ncol(tmp_table)] <- "maximum likelihood"
    colnames(tmp_table)[ncol(tmp_table)-1] <- "modal params"
    
    # Plot densities and iter plots
    tmp_chains <- NULL
    for(i in 1:length(tmp_mcmc_results$files)){
        tmp_chains[[i]] <- read.csv(tmp_mcmc_results$files[[i]], header=1)
        tmp_chains[[i]] <- tmp_chains[[i]][,c(1,which(param_table$fixed==0)+1)]
    }
    
    if(!is.null(topdir)){
        setwd(oldwd)
    }
    write.table(file=paste(filename, "_mcmc_summary.csv",sep=""),tmp_table,sep=",",col.names=NA)
    remove(tmp_table)
    
    # Only pass parameters that were included in MCMC fitting. First index should be "sampno" for iteration number
    #'mcmc_all_plots_multi(filename,tmp_chains,param_table,(burnin+adaptive_period),best_pars_plot)
    remove(tmp_chains)
    
    # Carry out Gelman diagnostics
    diagnostics_error <- mcmc_diagnostics(tmp_mcmc_results$chains,filename, param_table)
    
    if(!is.null(diagnostics_error)){
        if(VERBOSE){
            print(paste("Error in diagnostics: ",diagnostics_error,sep=""))
        }
    }
    
    files <- tmp_mcmc_results$files
    remove(tmp_mcmc_results)
    if(VERBOSE){
        print(files)
    }

    return(list(params=best_pars, times=temp_dat[,1],files=files))
}


#' Bottom level MCMC wrapper function
#' 
#' Calls the \code{\link{run_metropolis_MCMC}} function to carry out the MCMC algorithm. This function handles parallelisation of chains and some basic formatting.
#' @param data the data over which to calculate likelihoods
#' @param param_table table of parameter data as specified in \code{\link{load_param_table}}
#' @param iterations number of iterations for the MCMC algorithm
#' @param opt_freq how frequently the acceptance rate is adapted. Defaults to 50
#' @param thin thinning value for the MCMC chain. Default is 50
#' @param burnin the length of the burn in period. Defaults to 100
#' @param adaptive_period length of the adaptive period. Defaults to 1
#' @param nchain number of chains to run
#' @param popt the desired acceptance rate. Defaults to 0.44
#' @param tmp_filename the generic file name/location at which to save the MCMC chains. Note that the file extension will be appended to this string
#' @param LIKELIHOOD_FUNCTION a valid pointer to an R function which returns a single, log likelihood of the data given the current parameters
#' @param MODEL_FUNCTION a valid pointer to an R function which is used to evaluate the model for the current set of parameters
#' @param VERBOSE boolean flag for additional output. Defaults to FALSE
#' @param PARALLEL OPTIONAL boolean flag indicating if each MCMC chain should be run in parallel using doPar. Defaults to FALSE
#' @return returns a list containing the whole MCMC chains, MCMC summary objects (from coda), the file locations of the MCMC chains, and the maximum likelihood parameters
#' @export
#' @seealso \code{\link{run_metropolis_MCMC}}
MCMC_fit_single <- function(data,
                            param_table,
                            iterations=1000,
                            opt_freq=50,
                            thin=50,
                            burnin=100,
                            adaptive_period=1,
                            nchain=1,
                            popt=0.44,
                            tmp_filename,
                            LIKELIHOOD_FUNCTION,
                            MODEL_FUNCTION,
                            VERBOSE=FALSE,
                            PARALLEL=FALSE
                            ){
        
    # Parallel setup
    if(PARALLEL){
        if(nchain > 1){
            no_cores<- detectCores()
        }
        else {
            no_cores <- 2
        }
        if(.Platform$OS.type=="unix"){
            registerDoParallel(cores=(no_cores-1))
        } else {
            c1 <- makeCluster(no_cores-1)
            registerDoParallel(c1)
        }
    }
    
    # For each chain
    if(VERBOSE){
        print(paste("Number of chains: ", nchain,sep=""))
    }
    
    if(PARALLEL){
        tmp_chains <- foreach(i=1:nchain, .packages='mcmcJH') %dopar%{
            x <- paste(getwd(), "/", run_metropolis_MCMC(startvalue=rand.params(param_table),
                                                         iterations=iterations,
                                                         data=data,
                                                         param_table=param_table,
                                                         popt=popt,
                                                         opt_freq=opt_freq,
                                                         thin=thin,
                                                         burnin=burnin,
                                                         adaptive_period=adaptive_period,
                                                         LIKELIHOOD_FUNCTION,
                                                         MODEL_FUNCTION,
                                                         paste(tmp_filename,"_",i,sep="")
                                                         ),sep="")
        }
    }
    else {
        tmp_chains <- NULL
        for(i in 1:nchain){
            tmp_chains[[i]] <-  paste(getwd(), "/", run_metropolis_MCMC(startvalue=rand.params(param_table),
                                                                        iterations=iterations,
                                                                        data=data,
                                                                        param_table=param_table,
                                                                        popt=popt,
                                                                        opt_freq=opt_freq,
                                                                        thin=thin,
                                                                        burnin=burnin,
                                                                        adaptive_period=adaptive_period,
                                                                        LIKELIHOOD_FUNCTION,
                                                                        MODEL_FUNCTION,
                                                                        paste(tmp_filename,"_",i,sep="")
                                                                        ), sep="")
        }
    }
    
    
    # Reload all chains to return
    final_chains <- NULL
    # Find maximum likelihood parameters
    best_pars <- NULL
    best_lnlike <- -Inf

    #' Container to rbind all chains
    tmp_big_chain <- NULL
    
    for(i in 1:length(tmp_chains)){
        tmp <- read.csv(tmp_chains[[i]], header=1)
        # Remove burnin and only save parameter values (ie. remove sampleno and lnlikelihood)
        tmp_big_chain <- rbind(tmp_big_chain,tmp[(burnin+adaptive_period):nrow(tmp),2:ncol(tmp)])        
        final_chains[[i]] <- as.mcmc(tmp[(burnin+adaptive_period):nrow(tmp),2:ncol(tmp)])
        if(max(tmp[,ncol(tmp)]) > best_lnlike){
            best_pars <- tmp[which.max(tmp[,ncol(tmp)]),2:(ncol(tmp)-1)]
            best_lnlike <- max(tmp[,ncol(tmp)])
        }
    }

    #' Use the rbound chain to get density estimates and therefore modes
    colnames(tmp_big_chain) <- colnames(final_chains[[i]])

    modal_pars <- NULL

    for(index in 1:ncol(tmp_big_chain)){
        tmp_den <- tmp_big_chain[,index]
        z <- density(tmp_den)
        
        mode_par <- z$x[which.max(z$y)]
        modal_pars[index] <- mode_par
    }

    # Convert results to returnable format
    combined_mcmc <- as.mcmc.list(final_chains)
    mcmc_summary <- summary(combined_mcmc)
    
    if(PARALLEL){
        if(.Platform$OS.type=="unix"){
            registerDoParallel(cores=old_workers)
        } else {
            stopCluster(c1)
        }
    }
    
    return(list(chains=final_chains,summary=mcmc_summary,files=tmp_chains,params=best_pars, mode_params=modal_pars))
}

