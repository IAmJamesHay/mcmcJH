//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <iostream>
#include <fstream>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace std;
using namespace Rcpp;

int MAX_TITRE = 15;
//[[Rcpp::export]]
NumericMatrix predict_titres(NumericVector params, NumericVector times){
  double y0 = 0;
  double final_t;
  double t_i;
  double mu, dRF, tp, ts, m;
  double tmp;
  NumericVector::iterator t = times.begin();
  NumericVector Y(times.size());

  int j = 0;
  double trunc_lower = params[0];
  double first_infection = params[1];
  double max_t = times[times.size()-1];
  int no_infections=params.size()/6;
  NumericMatrix out(times.size(),2);

  for(int i=1;i <= no_infections; ++i){
    if(i==no_infections){
      final_t = max_t;
    }
    else {
      final_t = params[6*(i)+1];
    }

    t_i = params[6*(i-1)+1];
    mu = params[6*(i-1)+2];
    dRF = params[6*(i-1)+3];
    tp = params[6*(i-1)+4];
    ts = params[6*(i-1)+5];
    m = params[6*(i-1)+6];

    while(t != times.end() && *t <= final_t){
      if(*t <= t_i) tmp = 0;
      else if(*t > t_i && *t <= (t_i+tp)) tmp = (mu/tp)*(*t-t_i);
      else if(*t > (tp+t_i) && *t <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(*t) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(*t)+m*(t_i+tp+ts)+(1-dRF)*mu);
      tmp += y0;
      Y[j] = tmp;
      ++t;
      ++j;
    }
    
    y0 = Y[j-1];
    
    if(y0 < trunc_lower){
      y0 = trunc_lower;
    }
  }
  
  out(_,0) = times;
  out(_,1)=Y;
  return out;
}



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


//[[Rcpp::export]]
double scaletuning(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  if(step > 1) step = 0.99;
  if(step < 0) step = 0.01;
  return(step);
}

//[[Rcpp::export]]
int randNumber(const int n){return floor(unif_rand()*n);}

//[[Rcpp::export]]
double obs_error(int actual, int obs){
  int MAX_TITRE = 20;
  double S = 0.8;
  double EA = 0.1;
  if(actual == obs) return(S);
  if(actual == obs+1 || actual == obs-1) return(EA/2);
  double tmp = (1.0/(MAX_TITRE-2.0))*(1.0-S-EA);
  return(tmp);
}

//[[Rcpp::export]]
double posterior(NumericVector params, NumericMatrix data){
  NumericMatrix y(data.nrow(),data.ncol());
  double ln = 0;
  y = predict_titres(params, data(_,0));
  for(int i = 0; i < y.nrow();++i){
    if(y(i,1) < 0) y(i,1) = 0;
    if(y(i,1) > 20) y(i,1) = 20;
    ln += log(obs_error(floor(y(i,1)), floor(data(i,1))));
  }
  return(ln);
}



//[[Rcpp::export]]
string run_MCMC(NumericVector startvalue,
		NumericMatrix data,
		NumericMatrix param_table,
		int iterations, 
		double popt, 
		int opt_freq, 
		int thin, 
		int burnin, 
		int adaptive_period, 
		string filename,
		int save_block
		){
  ofstream csv_write; // Output stream for csv results
  double TUNING_ERROR = 0.1; // Allowable tuning error for adaptive step sizes
  double newprobab, probab, difflike; // Doubles for log-likelihoods and weights
  int no_recorded = 0; // Record number of recorded iterations
  int sampno = 1; // Record current sample number (note this may be different to no_recorded if thin != 1)
  int j = 0; // Record current parameter index
  string mcmc_chain_file = filename + "_chain.csv";

  // Vectors for parameters
  NumericVector step_sizes(startvalue.length());
  NumericVector current_params(startvalue.length());
  NumericVector sampnos(save_block);
  NumericVector lnlikes(save_block);
  NumericVector proposal(startvalue.length());
  vector<int> fixed;
  NumericVector log_proposal(startvalue.length());
  
  // Vectors for adaptive steps
  NumericVector tempaccepted(startvalue.length()); // Store total number of accepted proposals for each parameter
  NumericVector tempiter(startvalue.length()); // Store total number of proposals for each parameter
  NumericVector reset(startvalue.length()); // Vector of zeroes for reset of above vectors
  NumericVector pcur(startvalue.length()); // Vector of acceptance rates for each parameter

  NumericVector lower_bounds(param_table.nrow()); // Lower bounds for each parameter
  NumericVector upper_bounds(param_table.nrow()); // Upper bounds for each parameter

  // Matrices to store the parameter set at each time point
  NumericMatrix empty_chain(save_block, startvalue.length());
  NumericMatrix chain(save_block, startvalue.length());
  
  // Get control parameters from parameter table
  for(int i =0; i < startvalue.length();++i){
    if(param_table(i,1) == 0){
      fixed.push_back(i);
    }
    step_sizes[i] = param_table(i,4);
    lower_bounds[i] = param_table(i,2);
    upper_bounds[i] = param_table(i,3);
    cout << "Fixed: " << param_table(i,1) << endl;
    cout << "Step size: " << step_sizes[i] << endl;
    cout << "Lower bound: " << lower_bounds[i] << endl;
    cout << "Upper bound: " << upper_bounds[i] << endl;
  }

  // Open stream to csv file
  csv_write.open(mcmc_chain_file.c_str());
  
  // Get initial parameter values and likelihood
  current_params = clone(startvalue);
  probab = posterior(current_params, data);
  
  // Write first line to csv file
  csv_write << sampno << ",";
  for(int i = 0; i < (startvalue.length()); ++i){
    csv_write << startvalue[i] << ",";
  }
  csv_write << probab << endl;
  sampno += 1;
  int jj = 0;
  for(int i = 0; i < (iterations + adaptive_period + burnin); ++i){
    jj = randNumber(fixed.size());
    j = fixed[jj];
    proposal = clone(current_params);
    proposal(j) = proposalfunction(proposal(j),lower_bounds(j),upper_bounds(j),step_sizes(j));

    newprobab = posterior(proposal, data);
    
    difflike=newprobab-probab;
    if(R::runif(0,1) < exp(difflike) || difflike > 0){
      current_params = clone(proposal);
      probab = newprobab;
      tempaccepted(j) += 1;
    }
    tempiter(j) += 1; 
      
    // If an iteration to be saved
    if(sampno%thin ==0){
      // Write sampno, current params and probab to chain
      sampnos(no_recorded) = sampno;
      for(int x = 0; x < (current_params.size()); ++x){
	chain(no_recorded,x) = current_params[x];
      }
      lnlikes(no_recorded) = probab;
      no_recorded++;
      
      if(no_recorded >= save_block){
	// Write chain to csv file
	for(int x = 0; x < chain.nrow(); ++x){
	  csv_write << sampnos(x) << ",";
	  for(int q =0; q < chain.ncol();++q){
	    csv_write << chain(x,q) << ",";
	  }
	  csv_write << lnlikes(x);
	  csv_write << endl;
	}
	chain = empty_chain;
	no_recorded = 0;
      }
    }
    sampno++;

    // If in the adaptive period and an update step
    if(opt_freq != 0 && i < adaptive_period && (i+1)%opt_freq==0){

      // Update step sizes for normal parameters
      cout << "Current step sizes: " << endl;
      for(int x = 0; x < fixed.size();++x){
	int jjj = fixed[x];
	// get current acceptance rate
	if(tempiter(jjj) > 0){

	  pcur[jjj] = tempaccepted(jjj)/tempiter(jjj);
	  cout << "Pcur: " << pcur[jjj] << endl;
	  // If outside acceptable range:
	  if(pcur(jjj) < popt - (TUNING_ERROR*popt) || pcur(jjj) > popt + (TUNING_ERROR*popt)){
	    step_sizes(jjj) = scaletuning(step_sizes(jjj), popt, pcur(jjj));
	  }
	  cout << "Step size: " << step_sizes(jjj) << endl;
	}
      }
      cout << endl;
      tempaccepted = clone(reset);
      tempiter = clone(reset);
    }
  }
  csv_write.close();

  return(mcmc_chain_file);
}
