MCMC_main <- function(data_file, param_files, mcmc_param_file, main_dir="",infection_labels=NULL,consolidate_param_files){

    mcmc_params <- read.csv(mcmc_param_file,header=1)

    raw_data <- read.csv(data_file, stringsAsFactors=FALSE,na.strings=c("NA","-","?"),header=0)
    colnames(raw_data) <- raw_data[1,]
    raw_data <- raw_data[2:nrow(raw_data),]

    #' Data
    melted.data <- melt(raw_data)
    melted.data$variable <- as.integer(as.character(melted.data$variable))

    #' Save current WD and then move to specified WD
    oldwd <- getwd()
    dir.create(main_dir)
    setwd(main_dir)

    # Match up correct likelihood function
    mcmc_params$likelihood_function <- match.fun(as.character(mcmc_params$likelihood_function))
    mcmc_params$model_function <- match.fun(as.character(mcmc_params$model_function))

    file_names <- tolower(unlist(lapply(unique(melted.data$strain),function(x) unlist(strsplit(x,"/"))[2])))
    file_names <- gsub("[[:punct:]]|\\s+","_",file_names)
    names(file_names) <- unique(melted.data$strain)
    
    final <- MCMC_fit_1.1(
        main_dir,
        melted.data,
        file_names,
        param_files,
        mcmc_params$iterations,
        mcmc_params$opt_freq,
        mcmc_params$thin,
        mcmc_params$burnin,
        mcmc_params$adaptive_period,
        mcmc_params$nchain,
        mcmc_params$popt,
        mcmc_params$lower_plot_bound,
        mcmc_params$upper_plot_bound,
        mcmc_params$model_function,
        infection_times,
        FALSE,
        FALSE
        )
    
    consolidate_data(consolidate_param_files,getwd())
    setwd(oldwd)
    return(final)
    
}

MCMC_fit_1.1 <- function(
    top_dir,
    data,
    file_names,
    param_files,
    iterations,
    opt_freq,
    thin,
    burnin,
    adaptive_period,
    nchain,
    popt,
    lower_plot_bound,
    upper_plot_bound,
    MODEL_FUNCTION,
    infection_labels,
    SKIP_ZEROES = FALSE,
    SIMPLE_PLOTS=TRUE
    ){
    print(mem_used())
    #' Detect number of available cores for parallelism
    no_cores<- detectCores()

    #' Get operating system and call appropriate function
    if(.Platform$OS.type=="unix"){
        registerDoParallel(cores=(no_cores-1))
    } else {
        c1 <- makeCluster(no_cores-1)
        registerDoParallel(c1)
    }

    # Check given data. If any all zeroes, remove from further analysis
    if(SKIP_ZEROES){
        tmp <- remove_zero_data(data,file_strains)
        data <- tmp[[1]]
        file_strains <- tmp[[2]]
    }
    
    print(mem_used())
    #' For each group in the data
    index <- 1
    for(group in unique(data$group)){
        print("Fitting group: ")
        print(group)
        print(mem_used())
        #' Load parameter table for this group
        param_table <- load_param_table(param_files[[index]])
        infection_times <- infection_labels[infection_labels$group==group,c("time","infection")]
        #' Make a directory for that group, and a tmp directory for chains
        setwd(top_dir)
        dir.create(group,showWarnings=FALSE)
        tmp_dir <- getwd()
        setwd(group)
        chain_dir <- paste(getwd(),"/tmp/",sep="")
        dir.create(chain_dir)
        setwd(tmp_dir)
        
        #' Get number of unique strains for loops
        number_strains <- unique(data[data$group == group, "strain"])
        file_strains <- as.character(file_names[intersect(names(file_names), unique(data[data$group==group,"strain"]))])
        results_MCMC <- NULL
        print(file_strains)        
        #' Parallel fitting for each strain, within each group
        results_MCMC <- foreach(i =1:length(number_strains)) %dopar%{
            print("Strain: ")
            print(number_strains[i])

            setwd(group)

            #' Get subset of data for this group and strain
            temp_dat <- data[data$group==group & data$strain==number_strains[i],c("variable","value")]
            print(temp_dat[,1])
            #' Create a filename to be used for all MCMC files
            tmp_filename <- paste(group,"_",file_strains[i],sep="")

            #' The MCMC fitting procedure
            final <- MCMC_fit_single(
                temp_dat,
                param_table,
                iterations,
                opt_freq,
                thin,
                burnin,
                adaptive_period,
                nchain,
                popt,
                tmp_filename,
                chain_dir
                )

            setwd(tmp_dir)
            final            
        }
        print(mem_used())
        i <- 1
        for(result in results_MCMC){
            filenames <- result$files
            tmp_filename <- paste(group,"_",file_strains[i],sep="")
            best_pars_plot<- result$params
            modal_pars <- result$mode_pars

            tmp_chains <- read_chain_list(filenames, c("sampno",param_table$names, "lnlike"), 0, 100, c(1,(which(param_table$fixed==0)+1),(nrow(param_table)+2)))

            if(SIMPLE_PLOTS) simple_mcmc_plots(tmp_filename,tmp_chains,(burnin+adaptive_period))
            else mcmc_all_plots_multi(tmp_filename,tmp_chains,param_table,(burnin+adaptive_period),best_pars_plot[(which(param_table$fixed==0))])
            
            #' Carry out Gelman diagnostics
            tmp_chains <- read_chain_list(filenames, c("sampno",param_table$names, "lnlike"), 0, 1, seq(2,nrow(param_table)+1,by=1))
          #'  diagnostics_error <- mcmc_diagnostics(tmp_chains,tmp_filename, param_table)
            
            remove(tmp_chains)
            
          #'  if(!is.null(diagnostics_error)){
          #'      print(paste("Error in diagnostics: ",diagnostics_error,sep=""))
          #'  }
            i <- i + 1
        }
        print(mem_used())
        #' Generate prediction intervals from the MCMC chains
        prediction_intervals<- generate_prediction_intervals_1.1(
            results_MCMC,
            (burnin+adaptive_period),
            param_table,
            data[data$group==group,],
            MODEL_FUNCTION,
            group
            )

        #' Plot model fit and prediction intervals against data
        plot_model_fits(
            prediction_intervals$model_data,
            prediction_intervals$prediction_data,
            group,
            c(lower_plot_bound,upper_plot_bound),
            prediction_intervals$lower_prediction_bounds,
            prediction_intervals$upper_prediction_bounds,
            infection_times=infection_times
            )
        remove(results_MCMC)
        print(mem_used())
        
        index <- index + 1
    }
    print(mem_used())
    if(.Platform$OS.type=="unix"){
        registerDoParallel(cores=old_workers)
    } else {
        stopCluster(c1)
    }
    
    return(0)
    
}

MCMC_fit_single <- function(data, param_table, iterations, opt_freq, thin, burnin, adaptive_period, nchain,popt,filename,chain_dir){
    print(mem_used())
    #' FUNCTION TO FIND NUMBER OF INFECTIONS FROM PARAM TABLE?
    no_infections <- floor(length(which(param_table$nuisance==0))/6)
    print(paste("Number of infections: ", no_infections,sep=""))

    #' Get current directory
    tmp_dir <- getwd()

    #' Go into chain directory for chains
    setwd(chain_dir)
    
    tmp_chains <- NULL
    #' For each chain
    for(i in 1:nchain){
        print(paste("Chain number: ",i,sep=""))
        tmp_filename <- paste(filename, "_",i,sep="")
        #' Get number of rows to disregard for covariance matrix
        disgard <- (burnin+adaptive_period)/thin
        start <- rand.params(param_table)

        #' Runs a chain with univariate, uniform proposal to get starting covariance matrix for final fitting
        chain1_file <- paste(
            getwd(),
            "/",
            run_MCMC(
                start,
                as.matrix(data),
                as.matrix(param_table[,c("value","fixed","lower_bound","upper_bound","step","log_proposal")]),
                iterations,
                0.1,
                opt_freq,
                thin,
                burnin,
                adaptive_period,
                tmp_filename,
                opt_freq),
            sep=""
            )
        print("First chain done")

        #' Read in this chain, disregard burnin, and take parameter columns
        first_chain <- read.csv(chain1_file)
        first_chain <- first_chain[disgard:nrow(first_chain),2:(ncol(first_chain)-1)]

        #' Calculate covariance matrix
        cov_matrix <- cov(first_chain)

        #' Run main chain. Note that covariance matrix weight and initial step sizes for infection blocks are hardcoded here

        tmp_chains[[i]] <- paste(
            getwd(),
            "/",
            run_MCMC_test(
                startvalue=rand.params(param_table),
                data=as.matrix(data),
                param_table=as.matrix(param_table[,c("value","fixed","lower_bound","upper_bound","step","log_proposal")]),
                iterations,
                popt,
                opt_freq,
                thin,
                burnin,
                adaptive_period,
                tmp_filename,
                opt_freq,
                cov_matrix,
                rep(2.38,no_infections+1),
                0.8),
            sep=""
            )
        print("Second chain done")
    }
    
    setwd(tmp_dir)
    
    final_chains <- NULL
    tmp_big_chain <- NULL

    #' Using the saved filenames, read in all of the chains
    for(i in 1:length(tmp_chains)){
        tmp <- read.csv(tmp_chains[[i]],header=0)
        if(any(is.na(tmp))){
            tmp <- !is.na(tmp)
        }
        colnames(tmp) <- c("sampno",param_table$names,"lnlike")
        #' Disregard burnin and adaptive period
        disregard <- (burnin+adaptive_period)/thin
        tmp_big_chain <- rbind(tmp_big_chain,tmp[disregard:nrow(tmp),2:ncol(tmp)])

        #' Save chain as part of a list
        final_chains[[i]] <- as.mcmc(tmp[disregard:nrow(tmp),2:ncol(tmp)])
    }
    if(any(is.na(tmp_big_chain))){
        tmp_big_chain <- !is.na(tmp_big_chain)
    }
    
    #' Set column names of the big chain to match individual chain
    colnames(tmp_big_chain) <- colnames(final_chains[[i]])

    #' From all of the chains, get the MLE parameters
    best_pars <- as.numeric(tmp_big_chain[which.max(tmp_big_chain[,ncol(tmp_big_chain)]),])
    best_lnlike <- as.numeric(max(tmp_big_chain[,ncol(tmp_big_chain)]))
    print("Best pars: ")
    print(best_pars)

    #' Get modal parameters (excluding sampno)
    modal_pars <- NULL
    for(index in 1:ncol(tmp_big_chain)){
        tmp_den <- tmp_big_chain[,index]
        z <- density(tmp_den)
        mode_par <- z$x[which.max(z$y)]
        modal_pars[index] <- mode_par
    }

    #' Get summary of the MCMC list
    combined_mcmc <- as.mcmc.list(final_chains)
    mcmc_summary <- summary(combined_mcmc)

    #' Get and print effective sample size for chain
    effective_sample_size <- effectiveSize(combined_mcmc)
    print("Effective sample sizes: ")
    print(effective_sample_size)

    #' Save MCMC summary to the appropriate csv file
    tmp_table <- mcmc_summary$statistics

    print("here")
    
    tmp_table <- cbind(tmp_table, mcmc_summary$quantiles,effective_sample_size,c(as.numeric(modal_pars),NA),c(as.numeric(best_pars),NA))

    colnames(tmp_table)[ncol(tmp_table)] <- "maximum likelihood"
    colnames(tmp_table)[ncol(tmp_table)-1] <- "modal params"
    write.table(file=paste(filename, "_mcmc_summary.csv",sep=""),tmp_table,sep=",",col.names=NA)

    remove(combined_mcmc)
    remove(final_chains)
    remove(tmp_big_chain)
    remove(tmp)
    gc()
    print(mem_used())
    return(list(files=tmp_chains,params=best_pars, mode_params=modal_pars,times=data[,1]))
}
