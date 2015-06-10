#GLOBAL VARIABLES
MAX_TITRE <- 15
EPSILON <- 0.01
S <- 0.8
EA <- 0.08

# Generates an observation error matrix using global parameters.
# MAX_TITRE is the maximum possible log titre that can be observed
# S is the sensitivity of the observation. That is, the proportion of
# measurements that give the true titre value
# EA is the proportion of measurements that fall within one log-titre
# level of the true measurement
generate_observation_error_matrix <- function(){
    obs.matrix <- matrix(ncol=MAX_TITRE+1,nrow=MAX_TITRE+1)
    for(i in 1:(MAX_TITRE+1)){
        for(j in 1:(MAX_TITRE+1)){
            if(i == (MAX_TITRE+1) & j ==(MAX_TITRE+1)) obs.matrix[i,j] <- S + EA/2 - (1/(MAX_TITRE-2))*(1-S-EA)
            else if(i == 1 & j == 1) obs.matrix[i,j] <- S + EA/2 - (1/(MAX_TITRE-2))*(1-S-EA)
            else if(i==j) obs.matrix[i,j] <- S
            else if(i == (j+1) | i == (j-1)) obs.matrix[i,j] <- EA/2
            else obs.matrix[i,j] <- (1/(MAX_TITRE-2))*(1-S-EA)
        }
    }
    return(obs.matrix)
}

OBS.MATRIX <- generate_observation_error_matrix()


# Probability distribution of titres. Poisson distributed titres. Assume that titres must be between 0 and max.
# If passed a true titre value (lambda) outside of this, return 0, as we assume that this is not possible.
# If the observed value is on the border (ie. max), sum probability of everything from max onwards, as we assume
# that k could have been anything from 8 upwards.
f <- function(k, lambda){
    if(k != MAX_TITRE) final <- dpois(k,lambda)
    else{
        final <- sum(vapply(seq(MAX_TITRE,50,by=1),function(x) dpois(x,lambda),1))
    }
    return(final)
}

# Likelihood function for an observed titre value, obs, given a true mean titre, lambda
likelihood_complex<- function(lambda, obs){
    sum <- NULL
    sum <- sum(vapply(seq(0,MAX_TITRE,by=1),function(x) (f(x,lambda)*OBS.MATRIX[x+1,obs+1]),1))
    return(sum)
}

# Simplified version of the above function requiring no function calls
likelihood_simple<- function(lambda, obs){
    return((1-((MAX_TITRE+1)*EPSILON)/MAX_TITRE)*f(obs,lambda) + EPSILON/MAX_TITRE)
}


# Reads observation error from observation matrix
o <- function(k, obs){
    return(OBS.MATRIX[k+1,obs+1])
}

# Uniform observation model accounting for the potential to observing a titre different to the true one
# Based on model provided by Kucharski et al.
observation_model_kucharski<- function(k, obs){
    epsilon <- EPSILON
    if(k == obs) return(1-EPSILON)
    else return(EPSILON/MAX_TITRE)
}
