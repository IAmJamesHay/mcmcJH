
fnProposeParamUpdatesSingle <- function(ptab,index) {

# Takes biological parameters, a mask which is a list of the parameters

# being fitted, a table of the fitted parameters and their ranges and

# the number of parameters to be fitted.

# Returns a proposed vector of parameters.

# Comment in and out bouncing and cyclical boundary conditions

# for the random walk

bps <- ptab[,"val"]

rtn <- bps

if (index < 0 || index > dim(ptab)[1]) stop("index must be less than
or equal to number of parameters")

# Set up and transform to unit scale

rv <- runif(1)

rv <- (rv-0.5)* ptab[index,"step"]

x <- bps[index]

x <- SR_to_unit(x,min=ptab[index,"min"],max=ptab[index,"max"],logflag=ptab[index,"log"])

x <- x + rv

# Cyclical boundary conditions

if (x < 0) x <- 1 + x

if (x > 1) x <- x - 1

# Test for errors and return to originl scales

if (x < 0 || x > 1) stop("problem here")

rtn[index] <- SR_from_unit(x,min=ptab[index,"min"],max=ptab[index,"max"],logflag=ptab[index,"log"])

rtn

}

Here's an example:

# Parameter ordering is: R0, pc, I0

psi.A2.sim.fit.params <- function(

simdat,

nsamps=5) {

# Define bounds for parameters

ptab <- data.frame(

max = c( 5.0, 1.0, 10000 ),

min = c( 0.5, 0.01, 1 ),

log = c( FALSE, TRUE, TRUE ),

step = c( 0.1, 0.1, 0.1 ),

val = c( 1.4, 0.01, 10 )

)

# Simulate one single

maxtime <- length(simdat)

nps <- dim(ptab)[1]

# Set initial conditions

curpars <- ptab$val

initsol <- psi.sim.modelA(R0=ptab$val[1],pC=ptab$val[2],seed=ptab$val[3],stochastic=FALSE,tps=0:maxtime,plot=FALSE)

lnlike <- psi.mbm.lnlike(initsol,simdat)

# Set up measuring table

nMcmc <- nps * nsamps

tabMcmc <- data.frame(

index=1:nMcmc,

lnlike=rep(NA,nMcmc),

R0=rep(NA,nMcmc),

pc=rep(NA,nMcmc),

I0=rep(NA,nMcmc))

# Start the sampling loop

for (i in 1:nsamps) {

for (j in 1:nps) {

sampno <- (i-1)*nps + j

# Propose solution

newparams <- fnProposeParamUpdatesSingle(ptab,j)

# Calculate new loglike

newsol <- psi.sim.modelA(R0=newparams[1],pC=newparams[2],seed=newparams[3],stochastic=FALSE,tps=0:maxtime,plot=FALSE)

# Calculate loglike

newlnlike <- psi.mbm.lnlike(newsol,simdat)

diff_like <- newlnlike - lnlike

if (diff_like > 0) {

accept <- TRUE

} else if (exp(diff_like) > runif(1)) {

accept <- TRUE

} else accept <- FALSE

if (accept) {

ptab[,"val"] <- newparams

lnlike <- newlnlike

}

# Update the record of the chain

tabMcmc[sampno,"lnlike"] <- lnlike

tabMcmc[sampno,"R0"] <- ptab[1,"val"]

tabMcmc[sampno,"pc"] <- ptab[2,"val"]

tabMcmc[sampno,"I0"] <- ptab[3,"val"]

}

}

tabMcmc

}
