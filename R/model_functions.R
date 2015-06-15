#' ODE for ferret model
#'
#' ODE for Boosting and Waning model, used by the ode solver, deSolve.
#' @param t vector of time points over which dY should be calculated
#' @param x nothing - just required by deSolve I think...
#' @param params the vector of parameters for the model
titre_model <- function(t, x, params){
    mu <- params[1]
    dRF <- params[2]
    tp <- params[3]
    ts <- params[4]
    m <- params[5]
    dY <- ifelse(t < tp, mu/tp, 
                 ifelse(t >= tp & t< ts, -((1-dRF)*mu)/(ts-tp), 
                        ifelse(t >= ts, -m,0)))
    list(c(dY))
    
}



#' Ferret model function with lower bounded titres
#'
#' Uses the Ferret boost/wane model in the form y=mx+c to produce a matrix of model results from a given vector of parameters. This particular version allows the titres to be truncated from below, with the first argument of the parameter vector providing the lower bound
#' @param params the vector of parameters used to solve the model. The first argument should be the lower titre bound
#' @param times a vector of times over which to solve the model. If in doubt, use a sequence from 0 to max time at steps of 1
#' @return a matrix containing the titre values predicted from the given parameters at each time point specified by times
#' @export
#' @seealso \code{\link{predict.titre.fast.bounded}}, \code{\link{predict.titre.fast}}, \code{\link{predict_titre_universal_m_2}}, \code{\link{predict_titre_universal_m}}
predict.titre.fast.bounded <- function(params,times){
    trunc_lower <- params[1]
    params <- params[-1]
    first_infection <- params[1]
    out <- matrix(data=c(times,rep(0,length(times))),nrow=length(times),ncol=2)
    max_t <- times[length(times)]
    no_infections <- as.integer(length(params)/6)
    y0 <- 0
    
    # Get dynamics between each infection
    for(i in 1:no_infections){
      # Check if this is final infection. If so, get dynamics up to end. If not, get dynamics up until next infection
      if(i==no_infections) final_t <- max_t
      else final_t <- params[6*(i)+1]

      # Get current infection time
      t_i <- params[6*(i-1)+1]

      # Generate time frame, where length is between the current and next infection
      t <- c(times[times <= final_t & times > t_i],final_t)

      # Get parameters
      mu <- params[6*(i-1)+2]
      dRF <- params[6*(i-1)+3]
      tp <- params[6*(i-1)+4]
      ts <- params[6*(i-1)+5]
      m <- params[6*(i-1)+6]

      #' Use y=mx + c model to get predictions
      Y <-(
          (t <= t_i)*0
          ) +
              (
                  (t > t_i)*(t <= (tp + t_i))*((mu/tp)*t-(mu/tp)*t_i)
                  ) +
                      (
                          (t > (tp + t_i))*(t <= (ts + tp + t_i))*((-(dRF*mu)/ts)*t+((mu*dRF)/ts)*(t_i+tp) + mu)
                          ) +
                              (
                                  (t > (ts + tp + t_i))*(-m*t+m*(t_i+tp+ts)+(1-dRF)*mu)
                                  ) + y0
      
      #'Y <- (((t <= (tp + t_i))*(mu/tp)*(t-t_i)) +
      #'   ((t > (tp + t_i))*(t <= (ts + t_i))*((-(((1-dRF)*mu)/(ts-tp))*(t-t_i) + (((1-dRF)*mu)/(ts-tp))*tp + mu))) +
      #' ((t > (ts + t_i))*((m*ts-m*(t-t_i)) + (dRF*mu)))) + y0

      
      Y[Y<trunc_lower] <- trunc_lower
      y0 <- Y[length(t)]
      present <- out[,1] %in% t
      
      Y <- Y[1:length(out[present,2])]

      out[present,2] <- Y
  }
    return(out)
}


#' Ferret model function with lower bounded titres and a universal waning parameter, m
#'
#' Uses the Ferret boost/wane model in the form y=mx+c to produce a matrix of model results from a given vector of parameters. This particular version allows the titres to be truncated from below, and a universal waning parameter, m. The first argument should be m, and  the second argument of the parameter vector should provide the lower bound
#' @param params the vector of parameters used to solve the model. The first argument should be m, and the second argument the lower titre bound
#' @param times a vector of times over which to solve the model. If in doubt, use a sequence from 0 to max time at steps of 1
#' @return a matrix containing the titre values predicted from the given parameters at each time point specified by times
#' @export
#' @seealso \code{\link{predict.titre.fast.bounded}}, \code{\link{predict.titre.fast}}, \code{\link{predict_titre_universal_m_2}}, \code{\link{predict_titre_universal_m}}
predict_titre_universal_m_2<- function(params, times){
    m <- params[1]
    trunc_lower <- params[2]
    first_infection <- params[3]

    out <- matrix(data=c(times,rep(0,length(times))),nrow=length(times),ncol=2)
    max_t <- times[length(times)]
    params <- params[c(-1,-2)]
    no_infections <- as.integer((length(params))/5)
    y0 <- 0
    
                                        # Get dynamics between each infection
    for(i in 1:no_infections){
                                        # Check if this is final infection. If so, get dynamics up to end. If not, get dynamics up until next infection
        if(i==no_infections) final_t <- max_t
        else final_t <- params[5*(i)+1]

                                        # Get current infection time
        t_i <- params[5*(i-1)+1]

                                        # Generate time frame, where length is between the current and next infection
        t <- c(times[times <= final_t & times > t_i],final_t)

                                        # Get parameters
        mu <- params[5*(i-1)+2]
        dRF <- params[5*(i-1)+3]
        tp <- params[5*(i-1)+4]
        ts <- params[5*(i-1)+5]

        #' Use y=mx + c model to get predictions
        Y <-(
            (t <= t_i)*0
            ) +
                (
                    (t > t_i)*(t <= (tp + t_i))*((mu/tp)*t-(mu/tp)*t_i)
                    ) +
                        (
                            (t > (tp + t_i))*(t <= (ts + tp + t_i))*((-(dRF*mu)/ts)*t+((mu*dRF)/ts)*(t_i+tp) + mu)
                            ) +
                                (
                                    (t > (ts + tp + t_i))*(-m*t+m*(t_i+tp+ts)+(1-dRF)*mu)
                                    ) + y0
        
        #'Y <- (((t <= (tp + t_i))*(mu/tp)*(t-t_i)) +
        #'   ((t > (tp + t_i))*(t <= (ts + t_i))*((-(((1-dRF)*mu)/(ts-tp))*(t-t_i) + (((1-dRF)*mu)/(ts-tp))*tp + mu))) +
        #' ((t > (ts + t_i))*((m*ts-m*(t-t_i)) + (dRF*mu)))) + y0


        Y[Y<trunc_lower] <- trunc_lower
        y0 <- Y[length(t)]
        present <- out[,1] %in% t
        
        Y <- Y[1:length(out[present,2])]

        out[present,2] <- Y
    }
    return(out)
}

#' Old titre prediction function
#'
#' No longer active - this used ODEs rather than y=mx+c which was much too slow
#' @seealso \code{\link{predict.titre.fast.bounded}}, \code{\link{predict.titre.fast}}, \code{\link{predict_titre_universal_m_2}}, \code{\link{predict_titre_universal_m}}
predict.titre.OLD <- function(params,times){
    first_infection <- params[1]
    t <- c(times[times < first_infection], first_infection)

    out <- data.frame(time=t,y=0)
    max_t <- max(times)

    no_infections <- as.integer(length(params)/6)

    # Get dynamics between each infection
    for(i in 1:no_infections){
      # Check if this is final infection. If so, get dynamics up to end. If not, get dynamics up until next infection
      if(i==no_infections) final_t <- max_t
      else final_t <- params[6*(i)+1]

      # Get current infection time
      infection <- params[6*(i-1)+1]

      # Generate time frame, where length is between the current and next infection
      t <- c(times[times < final_t & times > infection], final_t)
      t <- t-infection
      # Starting titre is titre at time of infection as predicted by previous infection
      y0 <- out[out$time==infection,2]

      # Get parameters
      mu <- params[6*(i-1)+2]
      dRF <- params[6*(i-1)+3]
      tp <- params[6*(i-1)+4]
      ts <- params[6*(i-1)+5]
      m <- params[6*(i-1)+6]

      
      # Use y=mx + c model to get predictions
      Y <- ifelse(t <= tp, (mu/tp)*t, 
                  ifelse(t > tp & t<= ts, -(((1-dRF)*mu)/(ts-tp))*t + (((1-dRF)*mu)/(ts-tp))*tp + mu, 
                         ifelse(t > ts, m*ts-m*t + dRF*mu,0)))

      # Make into data frame. Make sure to add initial titre, y0!
      tmp <- data.frame(time=t+infection,y=Y + y0)

      out <- rbind(out, tmp)
    }
    return(out)
}



#' Ferret model function with a universal waning parameter, m
#'
#' Uses the Ferret boost/wane model in the form y=mx+c to produce a matrix of model results from a given vector of parameters. This particular version allows a universal waning parameter, m. The first argument should be m
#' @param params the vector of parameters used to solve the model. The first argument should be m
#' @param times a vector of times over which to solve the model. If in doubt, use a sequence from 0 to max time at steps of 1
#' @return a matrix containing the titre values predicted from the given parameters at each time point specified by times
#' @export
#' @seealso \code{\link{predict.titre.fast.bounded}}, \code{\link{predict.titre.fast}}, \code{\link{predict_titre_universal_m_2}}, \code{\link{predict_titre_universal_m}}
predict_titre_universal_m <- function(params, times){
    m <- params[1]
    first_infection <- params[2]
    out <- matrix(data=c(times,rep(0,length(times))),nrow=length(times),ncol=2)
    max_t <- times[length(times)]
    params <- params[-1]
    no_infections <- as.integer((length(params))/5)
    y0 <- 0
    
                                        # Get dynamics between each infection
    for(i in 1:no_infections){
                                        # Check if this is final infection. If so, get dynamics up to end. If not, get dynamics up until next infection
        if(i==no_infections) final_t <- max_t
        else final_t <- params[5*(i)+1]

                                        # Get current infection time
        t_i <- params[5*(i-1)+1]

                                        # Generate time frame, where length is between the current and next infection
        t <- c(times[times <= final_t & times > t_i],final_t)

                                        # Get parameters
        mu <- params[5*(i-1)+2]
        dRF <- params[5*(i-1)+3]
        tp <- params[5*(i-1)+4]
        ts <- params[5*(i-1)+5]

        #' Use y=mx + c model to get predictions
        Y <-(
            (t <= t_i)*0
            ) +
                (
                    (t > t_i)*(t <= (tp + t_i))*((mu/tp)*t-(mu/tp)*t_i)
                    ) +
                        (
                            (t > (tp + t_i))*(t <= (ts + tp + t_i))*((-(dRF*mu)/ts)*t+((mu*dRF)/ts)*(t_i+tp) + mu)
                            ) +
                                (
                                    (t > (ts + tp + t_i))*(-m*t+m*(t_i+tp+ts)+(1-dRF)*mu)
                                    ) + y0
        
        #'Y <- (((t <= (tp + t_i))*(mu/tp)*(t-t_i)) +
        #'   ((t > (tp + t_i))*(t <= (ts + t_i))*((-(((1-dRF)*mu)/(ts-tp))*(t-t_i) + (((1-dRF)*mu)/(ts-tp))*tp + mu))) +
        #' ((t > (ts + t_i))*((m*ts-m*(t-t_i)) + (dRF*mu)))) + y0

        Y[Y<0] <- 0
        y0 <- Y[length(t)]
        present <- out[,1] %in% t
        
        Y <- Y[1:length(out[present,2])]

        out[present,2] <- Y
    }
    return(out)
}

#' Basic ferret model function
#'
#' Uses the Ferret boost/wane model in the form y=mx+c to produce a matrix of model results from a given vector of parameters. 
#' @param params the vector of parameters used to solve the model. The first argument should be the lower titre bound first infection time
#' @param times a vector of times over which to solve the model. If in doubt, use a sequence from 0 to max time at steps of 1
#' @return a matrix containing the titre values predicted from the given parameters at each time point specified by times
#' @export
#' @seealso \code{\link{predict.titre.fast.bounded}}, \code{\link{predict.titre.fast}}, \code{\link{predict_titre_universal_m_2}}, \code{\link{predict_titre_universal_m}}
predict.titre.fast <- function(params,times){
    first_infection <- params[1]
    out <- matrix(data=c(times,rep(0,length(times))),nrow=length(times),ncol=2)
    max_t <- times[length(times)]
    no_infections <- as.integer(length(params)/6)
    y0 <- 0
    
    # Get dynamics between each infection
    for(i in 1:no_infections){
      # Check if this is final infection. If so, get dynamics up to end. If not, get dynamics up until next infection
      if(i==no_infections) final_t <- max_t
      else final_t <- params[6*(i)+1]

      # Get current infection time
      t_i <- params[6*(i-1)+1]

      # Generate time frame, where length is between the current and next infection
      t <- c(times[times <= final_t & times > t_i],final_t)

      # Get parameters
      mu <- params[6*(i-1)+2]
      dRF <- params[6*(i-1)+3]
      tp <- params[6*(i-1)+4]
      ts <- params[6*(i-1)+5]
      m <- params[6*(i-1)+6]
      
      #' Use y=mx + c model to get predictions
      Y <-(
          (t <= t_i)*0
          ) +
              (
                  (t > t_i)*(t <= (tp + t_i))*((mu/tp)*t-(mu/tp)*t_i)
                  ) +
                      (
                          (t > (tp + t_i))*(t <= (ts + tp + t_i))*((-(dRF*mu)/ts)*t+((mu*dRF)/ts)*(t_i+tp) + mu)
                          ) +
                              (
                                  (t > (ts + tp + t_i))*(-m*t+m*(t_i+tp+ts)+(1-dRF)*mu)
                                  ) + y0
      
      #'Y <- (((t <= (tp + t_i))*(mu/tp)*(t-t_i)) +
      #'   ((t > (tp + t_i))*(t <= (ts + t_i))*((-(((1-dRF)*mu)/(ts-tp))*(t-t_i) + (((1-dRF)*mu)/(ts-tp))*tp + mu))) +
      #' ((t > (ts + t_i))*((m*ts-m*(t-t_i)) + (dRF*mu)))) + y0
      

      Y[Y<0] <- 0
      y0 <- Y[length(t)]
      present <- out[,1] %in% t
      
      Y <- Y[1:length(out[present,2])]

      out[present,2] <- Y
  }
    return(out)
}
