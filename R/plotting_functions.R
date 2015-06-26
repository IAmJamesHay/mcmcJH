#' Model plotting function
#'
#' Provides a neat plot of the provided data, "best fit" line and prediction intervals.
#' @param model_data a data frame containing the best fit model trajectory. Should be three columns: variable (time), value and strain. Strain is simply a grouping variable and could be anything (eg. age group, location)
#' @param actual_data a data frame containing the actual data. Columns are same as model_data
#' @param group this is simply used for the file name of the plot. _fit.pdf will be appended
#' @param lower OPTIONAL coordinates of the lower prediction bound
#' @param upper OPTIONAL coordinates of the upper prediction bound
#' @param infection_times OPTIONAL data frame to mark vertical lines representing events. A data frame with columns time and infection, with the time and label respectively
#' @param plot_title OPTIONAL plot title
#' @return returns the ggplot object
#' @export
#' @seealso \code{\link{restrain_bounds}}, \code{\link{create_ploygons}}
plot_model_fits <- function(model_data, actual_data, group,plot_bounds,lower=NULL,upper=NULL,infection_times=NULL,plot_title=NULL){
    # Plot bounds
    low <- plot_bounds[1]
    high <- plot_bounds[2]
 
    # Filename - might need to change file type
    filename <- paste(group,"_fit.pdf",sep="")

    title <- ""
    if(!is.null(plot_title)){
        title <- plot_title
    }
    
    # Just check that model data is in correct range.
    model_data$value[model_data$value < low] <- low
    model_data$value[model_data$value > high] <- high

    xscale <- c(0,21,37,49,70, infection_times$time)
    xlabels <- c("0","21","37","49","70",paste("\n\n",infection_times$infection,sep=""))
    xlabel_colours <- c(rep("gray20",5),rep("red",nrow(infection_times)))
    xlabel_sizes <- c(rep(14,5),rep(10,nrow(infection_times)))
    
                                        # Get base plot
    plot <- ggplot() +
        scale_fill_brewer(palette="Dark2") +
            scale_colour_brewer(palette="Dark2") + 
                geom_point(data=actual_data,aes(x=variable,y=value,fill=strain,colour=strain,ymax=13),size=4,position=position_jitter(w=0.5,h=0.1)) + 
                    geom_line(data=model_data,aes(x=variable,y=value,fill=strain,colour=strain),size=0.8) +
                        xlab("Time (days)") +
                            scale_y_continuous(breaks=seq(low,high,by=1),limits=c(low,high),expand=c(0,0))+
                                ylab("Log Titre") +
                                    ggtitle(title)+
                                        theme(
                                            legend.justification=c(0,1),
                                            legend.position=c(0,1),
                                            text=element_text(size=16,colour="gray20"),
                                            plot.title=element_text(size=28),
                                            legend.text=element_text(size=14,colour="gray20"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            axis.line=element_line(colour="gray20"),
                                            axis.line.x = element_line(colour = "gray20"),
                                            axis.line.y=element_line(colour="gray20"),
                                            axis.text.x=element_text(colour=xlabel_colours,size=xlabel_sizes),
                                            plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),
                                            panel.background=element_blank(),
                                            axis.text.y=element_text(colour="gray20",size=14))+
                                                scale_x_continuous(breaks=xscale,labels=xlabels)
                                              
    # If specified, add vertical lines for infection times
    if(!is.null(infection_times)){
        for(i in 1:nrow(infection_times)){
            plot <- plot + geom_vline(xintercept=infection_times[i,"time"],colour="red",linetype="longdash",angle="90")
        }
        #plot <- plot +  annotate("text",x=(infection_times$time),y=15,label=infection_times$infection,hjust=0,color="red",size=4)
    }

    # If specified, constrain passed bounds and add as polygon to the plot
    if(!is.null(lower) & !is.null(upper)){
        lower <- restrain_bounds(lower,low,high)
        upper <- restrain_bounds(upper,low,high)
        bounds <- create_polygons(lower,upper)
        plot <- plot + geom_polygon(data=bounds,aes(x=x,y=y,fill=strain),alpha=0.2)
    }

    # Saves the plot as the specified file. DPI SHOULD BE CHANGED IF FILE IS TOO BIG
    ggsave(plot=plot,filename=filename, width=14, height=10)
    return(plot)
}

#' Polygon for prediction intervals
#'
#' Takes upper and lower bounds as data frames and uses these to produce a polygon for
# ggplot. The idea here is simply to combine the two bounds into a single set of coordinates. Note
#' that this function expects a LIST of bounds (ie. a different polygon for each index in the list),
#' so a single polygon should come from lower[[1]] and upper[[1]]
#' @param lower a LIST of lower bound coordinates
#' @param upper a LIST of upper bound coordinates
#' @return a set of coordinates to draw a polygon representing the given prediction intervals
#' @export
#' @seealso \code{\link{plot_model_fits}}
create_polygons <- function(lower,upper){
    bounds <- NULL
    for(i in 1:length(upper)){
            if(length(upper) > 0 & length(lower) > 0){
                tmp <- as.data.frame(upper[[i]])
                tmp <- tmp[rev(rownames(tmp)),]
                tmp1 <- rbind(as.data.frame(lower[[i]]),tmp)
                bounds <- rbind(bounds, tmp1)
            }
        }
    colnames(bounds) <- c("x","y","strain","group")
    return(bounds)    
}

#' Restrains values in a data frame by the given bounds
#' 
#' Takes a list of data frames (ie. titre data) and makes sure that all y values are between lower and upper
#' @param dat a data frame to be bounded (with a column y)
#' @param lower the lower bound
#' @param upper the upper bound
#' @return the same data frame but with bounded y values
#' @export
restrain_bounds <- function(dat, lower, upper){
    for(j in 1:length(dat)){
        tmp <- dat[[j]]$y
        for(i in 1:length(tmp)){
            if(tmp[i] < lower){
                tmp[i] <- lower
            }
            else if(tmp[i] > upper){
                tmp[i] <- upper
            }
        }
        dat[[j]]$y <- tmp
    }
    return(dat)
}

#' Individual MCMC density plots with multiple chains
#'
#' Plots the MCMC density plots for the given parameter and data
#' @param name the name of the parameter from the data to be plotted
#' @param data the whole set of density data. Has a column to separate data by chain number
#' @param xlims the x axis bounds for the density plot
#' @param prior OPTIONAL set of data to plot the prior distribution
#' @return a single ggplot object of the density plot
#' @export
#' @seealso \code{\link{mcmc_all_plots_mcmc}}
mcmc_density_multi <- function(name, data, xlims, prior=NULL,best_fits=NULL){
    dat <- data[data$variable==name,]
    z <- density(dat[,2])
    mean_line <- mean(dat[,2])
    mode_line <- z$x[which.max(z$y)]

    q <- ggplot(data[data$variable==name,],aes(x=value,fill=chain,group=chain,y=..density..)) + geom_density(size=1,alpha=0.5) + ggtitle(paste(name, " Density Plot", sep="")) + scale_x_continuous(limits=xlims) +
        geom_vline(xintercept=mean_line,colour="red") +
            geom_text(aes_q(x=mean_line,label="\nMean",y=max(z$y/2)),colour="red",angle=90,text=element_text(size=6)) +
                geom_vline(xintercept=mode_line,colour="blue") +
                    geom_text(aes_q(x=mode_line,label="\nMode",y=max(z$y/2)),colour="blue",angle=90,text=element_text(size=6))
    if(!is.null(best_fits)){
        mle_line <- as.numeric(best_fits[which(names(best_fits)==name)])
        q <- q + geom_vline(xintercept=mle_line,colour="purple") +
                    geom_text(aes_q(x=mle_line,label="\nMLE",y=max(z$y/2)),colour="purple",angle=90,text=element_text(size=6))
    }
    if(!is.null(prior)){
        prior <- rbind(c(xlims[1],0.0,"prior"),prior,c(xlims[2],0,"prior"))
        prior$variable <- as.numeric(prior$variable)
        prior$value <- as.numeric(prior$value)
        q <- q + geom_polygon(data=prior, aes(x=variable, y=value,fill="prior"),alpha=0.5,size=1,color="black")
    }
    return(q)
}

#' Individual MCMC iteration plots with multiple chains
#'
#' Plots the MCMC chain for a given parameter and data
#' @param name the name of the parameter to be plotted
#' @param data the whole set of data (the MCMC chain). Has a column to separate data by chain number
#' @param burnin the length of the burnin period for a vertical line
#' @return a single ggplot object of the iter plot
#' @export
#' @seealso \code{\link{mcmc_all_plots_mcmc}}
mcmc_iter_multi <- function(name, data,burnin,best_fit=NULL){
    tmp_dat <- data[,c("iteration",name,"chain")]
    colnames(tmp_dat) <- c("iteration","value","chain")
    
    z <- density(tmp_dat[,"value"])
    mean_line <- mean(tmp_dat[,"value"])
    mode_line <- z$x[which.max(z$y)]
    
    q <- ggplot(tmp_dat,aes(x=iteration,y=value,colour=chain)) + geom_line() + ggtitle(paste(name, " Iter Plot",sep="")) + geom_vline(xintercept=burnin, colour="green", linetype="longdash")+
        geom_hline(yintercept=mean_line,colour="red") +
            geom_text(aes_q(y=mean_line,label="\nMean",x=max(z$x/2)),colour="red",text=element_text(size=6)) +
                geom_hline(yintercept=mode_line,colour="blue") +
                    geom_text(aes_q(y=mode_line,label="\nMode",x=max(z$x/2)),colour="blue",text=element_text(size=6))
    if(!is.null(best_fit)) {
        mle_line <- as.numeric(best_fit[which(names(best_fit)==name)])
        q <- q + geom_hline(yintercept=mle_line,colour="purple") +            
            geom_text(aes_q(y=mle_line,label="\nMLE",x=max(z$x/2)),colour="purple",text=element_text(size=6))
    }
    
}


#' Generates data for prior density plots
#'
#' Given a vector of parameter names and the parameter table as loaded by \code{\link{load_param_table}}, returns a data frame that is used to generate density plots of the prior distributions
#' @param names a vector of parameter names
#' @param the parameter table as loaded by \code{\link{load_param_table}}
#' @return a data frame of combined prior densities, with column names "variable (x coord)", "value (y coord)", "param (parameter name)", "chain" (all set to 'prior')
#' @export
generate_prior_data <- function(names, param_table){
    prior_dat <- NULL
    to_check <- which(param_table$fixed==0)
    for(i in 1:length(to_check)){
        lower_lim <- param_table[to_check[i],"lower_bound"]
        upper_lim <- param_table[to_check[i],"upper_bound"]
        tmp <- NULL
        samples <- seq(lower_lim, upper_lim, by=((upper_lim-lower_lim)/100))
        tmp_pars <- param_table$value
        for(j in 1:length(samples)){
            tmp_pars[to_check[i]] <- samples[j]
            tmp[j] <- prior_wrapper(param_table$prior_func[[to_check[i]]], tmp_pars, to_check[i], param_table$prior_args[[to_check[i]]])
        }
        if(sum(tmp) > 1 | sum(tmp) < 0){
            tmp <- exp(tmp)
        }
        tmp <- data.frame(variable=samples,value=tmp, param=param_table$names[to_check[i]],chain="prior")
        prior_dat <- rbind(prior_dat, tmp)
    }
    return(prior_dat)
}

#' Plot all MCMC density and iteration plots
#'
#' Saves a pdf file as 'filename'_MCMC_plots.pdf, which contains density and iteration plots for all of the provided parameters.
#' @param filename the overall file name/location to be used, appended with _MCMC_plots.pdf
#' @param mcmc_chains the actual MCMC chains to be used to create density plots and iteration plots
#' @param param_table OPTIONAL the parameter table as loaded by \code{\link{load_param_table}}, which is then used to create prior density plots
#' @param burnin OPTIONAL the burn in period, used to draw a vertical line
#' @return nothing
#' @export
#' @seealso \code{\link{mcmc_iter_multi}}, \code{\link{mcmc_density_multi}}, \code{\link{mcmc_density_single}}, \code{\link{mcmc_iter_single}}
mcmc_all_plots_multi <- function(filename, mcmc_chains, param_table=NULL,burnin=NULL,best_fit=NULL){
    tmp_filename <- paste(filename, "_MCMC_plots.pdf", sep="")

    # For iter
    tmp_all <- NULL
    for(i in 1:length(mcmc_chains)){
        tmp <- as.data.frame(mcmc_chains[[i]])
        colnames(tmp)[1] <- "iteration"
        tmp$chain <- as.character(i)
        tmp_all <- rbind(tmp_all, tmp)
    }
    colnames(tmp_all) <- colnames(tmp)

    if(!is.null(param_table)){
       # Generate data for prior plots
        prior_dat <- generate_prior_data(colnames(mcmc_chains[[1]]),param_table)
    }
    # For densities
    melted <- NULL
    for(i in 1:length(mcmc_chains)){
        tmp <- as.data.frame(mcmc_chains[[i]])
        tmp_melt <- melt(tmp)
        tmp_melt$chain <- as.character(i)
        melted <- rbind(melted, tmp_melt)
    }

    pdf(tmp_filename, height=6,width=16)
    for(i in 2:ncol(mcmc_chains[[1]])){
        print(suppressWarnings(grid.arrange(
            mcmc_iter_multi(colnames(mcmc_chains[[1]])[i],tmp_all,burnin,
                            best_fit),
            mcmc_density_multi(colnames(mcmc_chains[[1]])[i],melted,
                               c(param_table[param_table$names==colnames(mcmc_chains[[1]])[i],"lower_bound"],param_table[param_table$names==colnames(mcmc_chains[[1]])[i],"upper_bound"]),
                               prior_dat[prior_dat$param==colnames(mcmc_chains[[1]])[i],c("variable","value","chain")],
                               best_fit
                               )
            ,ncol=2)))
    }
    dev.off()
}

#' Individual MCMC density plots for a single chain
#'
#' Plots the MCMC density plots for the given parameter and data
#' @param name the name of the parameter from the data to be plotted
#' @param data the whole set of density data
#' @return a single ggplot object of the density plot
#' @export
#' @seealso \code{\link{mcmc_all_plots_mcmc}}
mcmc_density_single <- function(name,data){
    q <- ggplot(data[data$variable==name,],aes(x=value)) + geom_histogram(aes(y=..density..),fill="deepskyblue") + geom_density(alpha=1,colour="black",size=1) + ggtitle(paste(name, " Density Plot", sep=""))
    return(q)
}

#' Individual MCMC iter plots for a single chain
#'
#' Plots the MCMC iter plots for the given parameter and data
#' @param name the name of the parameter from the data to be plotted
#' @param data the whole set of iteration data
#' @return a single ggplot object of the density plot
#' @export
#' @seealso \code{\link{mcmc_all_plots_mcmc}}
mcmc_iter_single <- function(name, data){
    tmp_dat<- data.frame(value = data,iteration=seq(1,length(data),by=1))
    q <- ggplot(tmp_dat,aes(x=iteration,y=value),colour="deepskyblue") + geom_line(colour="deepskyblue") + ggtitle(paste(name," Iter Plot", sep=""))
    return(q)
}

#' Plot all MCMC density and iteration plots with a single chain
#'
#' Saves a pdf file as 'filename'_MCMC_plots.pdf, which contains density and iteration plots for all of the provided parameters.
#' @param filename the overall file name/location to be used, appended with _MCMC_plots.pdf
#' @param mcmc_chains the actual MCMC chains to be used to create density plots and iteration plots
#' @return nothing
#' @export
#' @seealso \code{\link{mcmc_iter_multi}}, \code{\link{mcmc_density_multi}}, \code{\link{mcmc_density_single}}, \code{\link{mcmc_iter_single}}
mcmc_all_plots_single <- function(filename, mcmc_chain){
    tmp_filename <- paste(filename, "_MCMC_plots.pdf",sep="")
    tmp_dat <- mcmc_chain
    tmp_melt <- melt(tmp_dat)
    
    pdf(tmp_filename,height=5,width=12)
    
    for(i in 1:ncol(mcmc_chain)){
        print(suppressWarnings(grid.arrange(mcmc_iter_single(colnames(tmp_dat)[i],tmp_dat[,i]), mcmc_density_single(colnames(tmp_dat)[i],tmp_melt),ncol=2)))
    }
    dev.off()
}


