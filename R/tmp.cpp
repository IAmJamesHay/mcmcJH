//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}
//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}

//[[Rcpp::export]]
double proposalfunction(double current, double lower, double upper, double step){
  double x,rv,rtn;
  x = toUnitScale(current, lower, upper);
  
  rv = (R::runif(0,1)-0.5)*step;
  x = x + rv;
  if (x < 0) x = -x;
  if(x > 1) x = 2-x;

  rtn = fromUnitScale(x,lower,upper);
  return(rtn);
}

// This is just the central part of the MCMC algorithm
// fixed is the vector of model parameters - it's the same vector as "current_params"
for(int i = 0; i < (iterations + adaptive_period + burnin); ++i){
  jj = randNumber(fixed.size()); // Choose a range index from the model parameters
  j = fixed[jj];
  proposal = clone(current_params); // Create a copy of the current parameters
  proposal(j) = proposalfunction(proposal(j),lower_bounds(j),upper_bounds(j),step_sizes(j)); // Generate a proposed move
  newprobab = posterior(proposal, data); // Calculate the posterior of this new move
  
  difflike=newprobab-probab;

  // NOTE: I am maximising the positive log likelihood, whereas you are minimising the negative log likelihood. I don't think this matters, as the signs and inequalities in our respective if statements make it equivalent
  if(R::runif(0,1) < exp(difflike) || difflike > 0 && proposal(0) < proposal(1)){
    current_params = clone(proposal);
    probab = newprobab;
    tempaccepted(j) += 1;
  }
  tempiter(j) += 1; 
 }

// Rest of MCMC code follows....
