// rcpp_functions.cpp
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace std;
using namespace Rcpp;

int MAX_TITRE = 15;


// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   bool success = false;
   arma::mat Y = arma::randn(n, ncols);
   arma::mat R;
   while(!success){
     success = arma::chol(R, sigma);
     if(!success){
       cout << "Error in mvnorm" << endl;
       sigma += arma::eye(sigma.n_rows, sigma.n_rows)*1e-3;
     }
   }
   return arma::repmat(mu, 1, n).t() + Y * R;
}


//[[Rcpp::export]]
arma::mat subset_mat(arma::mat a, arma::uvec b){
  return(a.submat(b,b));
}


//[[Rcpp::export]]
NumericMatrix predict_titres(NumericVector params, NumericVector times){
  double y0 = 0;
  double final_t;
  double t_i;
  double mu, dRF, tp, ts, m;
  double tmp;
  double tmp2;
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
      tmp = 0;
      if(*t <= t_i) tmp = 0;
      else if(*t > t_i && *t <= (t_i+tp)) tmp = (mu/tp)*(*t-t_i);
      else if(*t > (tp+t_i) && *t <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(*t) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(*t)+m*(t_i+tp+ts)+(1-dRF)*mu);
      tmp += y0;
      
      Y[j] = tmp;
      ++t;
      ++j;
    }
    
    if((6*i + 1) < params.size()){
      tmp2 = params[6*i + 1];
      if(tmp2 <= t_i) tmp = 0;
      else if(tmp2 > t_i && tmp2 <= (t_i+tp)) tmp = (mu/tp)*(tmp2-t_i);
      else if(tmp2 > (tp+t_i) && tmp2 <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(tmp2) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(tmp2)+m*(t_i+tp+ts)+(1-dRF)*mu);
      y0 = y0 + tmp;
    }
    
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
double proposal_function(double current, double lower, double upper, double step){
  double update;
  double move;
  double new1;
  new1 = toUnitScale(current,lower,upper);
  
  do {
    new1 = toUnitScale(current,lower,upper);
    update = R::rnorm(0, 1);
    new1 = new1 + update*step;
  } while(new1 > 1 || new1 < 0);
  
  new1 = fromUnitScale(new1,lower,upper);
  
  return(new1);
}


//[[Rcpp::export]]
NumericVector stats_proposal(NumericVector current, NumericVector upper, NumericVector lower, NumericVector steps, NumericVector fixed){
  double x, rv;
  NumericVector proposal(current.size());
  proposal = clone(current);

  for(int i =0; i < fixed.size();++i){
    x = toUnitScale(current(fixed[i]),lower(fixed[i]),upper(fixed[i]));

    rv = R::runif(0,1)-0.5;
    x = x + rv*steps(fixed[i]);
    
    if(x < 0) x = -x;
    if(x > 1) x = 2-x;

    proposal(fixed(i)) = fromUnitScale(x,lower(fixed[i]),upper(fixed[i]));
  }
  return(proposal);    
}


//[[Rcpp::export]]
NumericVector s_and_ea_proposal(NumericVector current, NumericVector upper, NumericVector lower, NumericVector steps){
  double S = current(0);
  double EA = current(1);
  NumericVector proposal(current.size());
  proposal = clone(current);
  S = S + (R::runif(0,1)-0.5)*steps(0);
  EA = EA + (R::runif(0,1)-0.5)*steps(1);

  if(S < 0) S = -S;
  if(S > 1) S = 2-S;

  if(EA < 0) EA = -EA;
  if(EA > 1) EA = 2-EA;
  
  proposal(0) = S;
  proposal(1) = EA;

  return(proposal);
}

//[[Rcpp::export]]
NumericVector mvnorm_proposal(NumericVector current, arma::mat covar_mat, NumericVector upper, NumericVector lower, NumericVector fixed){
  NumericVector proposed(current.size());
  proposed = clone(current);
  NumericVector x(fixed.size());
  NumericVector zeroes(fixed.size());
  proposed = clone(current);
  int j = 0;
  int ii = 0;
  while(j==0){
    j = 1;
    x = mvrnormArma(1,zeroes,covar_mat);
    for(int i = 0; i < fixed.size();++i){
      ii = fixed[i];
      proposed(ii) = current(ii) + x(i);
    }  
  }
  return(proposed);
}


//[[Rcpp::export]]
double scaletuning2(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  return(step);
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
double obs_error(int actual, int obs, double S, double EA){
  int MAX_TITRE = 13;
  //  double S = 0.95;
  //  double EA = 0.02;
  if(actual == (MAX_TITRE) && obs == (MAX_TITRE)) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==0 && obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==obs) return(S);
  else if(actual == (obs + 1) || actual==(obs-1)) return(EA/2.0);
  return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
}

//[[Rcpp::export]]
double discrete_norm(int a, int b, double sd){
  double c;
  c = R::pnorm((a+0.5),b,sd,1,0) - R::pnorm((a-0.5),b,sd,1,0);
  return(c);
}

//[[Rcpp::export]]
NumericVector tail_test(NumericVector a, int b){
  return(tail(a,a.size()-b));
}


//[[Rcpp::export]]
double posterior(NumericVector params, NumericMatrix data){
  NumericMatrix y(data.nrow(),data.ncol());
  double ln = 0;
  y = predict_titres(tail(params,params.size()-2), data(_,0));
  //cout << y.nrow() << endl;
  for(int i = 0; i < y.nrow();++i){
    //cout << y(i,1) << " " << data(i,1) << endl;
    if(y(i,1) < 0) y(i,1) = 0;
    if(y(i,1) >= 13) y(i,1) = 13;
    //ln += log(discrete_norm(floor(y(i,1)),floor(data(i,1)),params(0)));
    ln += log(obs_error(floor(y(i,1)), floor(data(i,1)),params(0),params(1)));
    //cout << ln << endl;
    //ln += log(R::dnorm(y(i,1),data(i,1),1,0));
    //ln += R::dpois(floor(data(i,1)),y(i,1),1);
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
    //proposal(j) = proposal_function(proposal(j),lower_bounds(j), upper_bounds(j),step_sizes(j));
    newprobab = posterior(proposal, data);
    
    difflike=newprobab-probab;
    if(R::runif(0,1) < exp(difflike) || difflike > 0 && proposal(0) < proposal(1)){
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
      for(int x = 0; x < fixed.size();++x){
	int jjj = fixed[x];
	// get current acceptance rate
	if(tempiter(jjj) > 0){

	  pcur[jjj] = tempaccepted(jjj)/tempiter(jjj);
	  // If outside acceptable range:
	  if(pcur(jjj) < popt - (TUNING_ERROR*popt) || pcur(jjj) > popt + (TUNING_ERROR*popt)){
	    step_sizes(jjj) = scaletuning(step_sizes(jjj), popt, pcur(jjj));
	  }
	}
      }
      tempaccepted = clone(reset);
      tempiter = clone(reset);
    }
  }
  csv_write.close();

  return(mcmc_chain_file);
}

//[[Rcpp::export]]
NumericMatrix create_blocks(NumericVector a){
  int no_infections = (a.size()-2)/6;
  NumericMatrix blocks(no_infections+1,2);
  blocks(0,0) = 0;
  blocks(0,1) = 2;
  for(int i = 1; i < blocks.nrow();++i){
    blocks(i,0) = 3 + 6*(i-1);
    blocks(i,1) = 8 + 6*(i-1);
  }
  return(blocks);
}
  

//[[Rcpp::export]]
string run_MCMC_test(NumericVector startvalue,
		     NumericMatrix data,
		     NumericMatrix param_table,
		     int iterations, 
		     double popt, 
		     int opt_freq, 
		     int thin, 
		     int burnin, 
		     int adaptive_period, 
		     string filename,
		     int save_block,
		     NumericMatrix start_covmatrix,
		     NumericVector scales,
		     double w
		     ){
  ofstream csv_write; // Output stream for csv results
  double TUNING_ERROR = 0.2; // Allowable tuning error for adaptive step sizes
  double newprobab, probab, difflike; // Doubles for log-likelihoods and weights
  int no_recorded = 0; // Record number of recorded iterations
  int sampno = 1; // Record current sample number (note this may be different to no_recorded if thin != 1)
  int j = 0; // Record current parameter index
  string mcmc_chain_file = filename + "_chain.csv";
  bool negative = false;
  int no_infections = (startvalue.size()-2)/6;
  int rand_i = 0;
  bool no_stats = false;
  
  // Create matrix that stores the indices of the blocks (one per infection and one for stats parameters)

  NumericMatrix blocks(no_infections+1,2);
  blocks = create_blocks(startvalue);

  // Vectors for parameters
  NumericVector current_params(startvalue.length());
  NumericVector sampnos(save_block);
  NumericVector lnlikes(save_block);
  NumericVector proposal(startvalue.length());
  NumericVector unfixed;
  NumericVector log_proposal(startvalue.length());
  
  // Vectors for adaptive steps
  NumericVector tempaccepted(blocks.nrow()); // Store total number of accepted proposals for each parameter
  NumericVector tempiter(blocks.nrow()); // Store total number of proposals for each parameter
  NumericVector reset(blocks.nrow()); // Vector of zeroes for reset of above vectors
  NumericVector pcur(blocks.nrow()); // Vector of acceptance rates for each parameter

  NumericVector lower_bounds(param_table.nrow()); // Lower bounds for each parameter
  NumericVector upper_bounds(param_table.nrow()); // Upper bounds for each parameter
  NumericVector steps(param_table.nrow());

  // Matrices to store the parameter set at each time point
  NumericMatrix empty_chain(save_block, startvalue.length());
  NumericMatrix chain(save_block, startvalue.length());
  
  // Get control parameters from parameter table
  for(int i =0; i < startvalue.length();++i){
    if(param_table(i,1) == 0){
      unfixed.push_back(i);
    }
    lower_bounds[i] = param_table(i,2);
    upper_bounds[i] = param_table(i,3);
    steps[i] = param_table(i,4);
  }
 
  // one 3x3 covar matrix for S, EA and lower bound. Then NxN matrix for each infection, where N is number unfixed parameters per infection
  // Cov matrix for stats parameters
  NumericVector tmp = unfixed[unfixed <= 2];
  int cov_size = tmp.size();
  
  // Stats matrix might all be fixed, so this size would be 0. If so, need to artificially increase this value to prevent errors.
  if(cov_size <= 0){
    cov_size = 3;
    no_stats = true;
    cout << "No stats parameters to be explored" << endl;
  }
  

  // Vectors for adaptive steps of stats parameters
  //NumericVector tempaccepted_stats(cov_size);
  //NumericVector tempiter_stats(cov_size);
  // NumericVector pcur_stats(cov_size);
  // NumericVector reset_stats(cov_size);

  arma::mat stats_cov_matrix(cov_size,cov_size);
  arma::mat stats_cov_matrix_scaled(cov_size,cov_size);
  stats_cov_matrix = (as<arma::mat>(start_covmatrix)).submat(0,0,cov_size-1,cov_size-1);
  stats_cov_matrix_scaled = stats_cov_matrix*scales(0);

  // Cov matrix for infection parameters
  tmp = unfixed [unfixed > 2 & unfixed <= 8];
  cov_size = tmp.size();  
  tmp = unfixed[unfixed <= 2];

  arma::cube cov_matrix(cov_size,cov_size,no_infections);
  arma::cube cov_matrix_scaled(cov_size,cov_size,no_infections);

  for(int i = 0; i < no_infections;++i){
    cov_matrix.slice(i) = (as<arma::mat>(start_covmatrix)).submat(as<arma::uvec>(unfixed[unfixed <= blocks(i+1,1) & unfixed >= blocks(i+1,0)]),as<arma::uvec>(unfixed[unfixed <= blocks(i+1,1) & unfixed >= blocks(i+1,0)]));
    cov_matrix_scaled.slice(i) = cov_matrix.slice(i)*scales(i+1);
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
    // Get the block to be updated
    if(no_stats){
      jj = randNumber(blocks.nrow()-1)+1;
    }
    else {
      jj = randNumber(blocks.nrow());
    }
    proposal = clone(current_params);
    if(jj == 0){
      proposal = mvnorm_proposal(proposal,stats_cov_matrix_scaled,upper_bounds,lower_bounds,unfixed[unfixed>=blocks(jj,0)&unfixed<=blocks(jj,1)]);
    }
    else{
      proposal = mvnorm_proposal(proposal,cov_matrix_scaled.slice(jj-1),upper_bounds,lower_bounds,unfixed[unfixed >= blocks(jj,0) & unfixed <= blocks(jj,1)]);
    }
    newprobab = posterior(proposal, data);
    difflike=newprobab-probab;

    // Checks if proposed parameters are outside allowable range, and automatically rejects if so
    negative = false;
    for(int qq = 0; qq < unfixed.size();++qq){
      int iiii = unfixed(qq);
      if(proposal(iiii) < lower_bounds(iiii) || proposal(iiii) > upper_bounds(iiii)){
	negative = true;
      }
    }
    if((proposal(0) + proposal(1)) >= 1) negative = true;
    if(proposal(1) > proposal(0)) negative = true;
    
    if((R::runif(0,1) < exp(difflike) || difflike > 0) && !negative){
      current_params = clone(proposal);
      probab = newprobab;
      tempaccepted(jj) += 1.0;
    }
    tempiter(jj) += 1.0; 
    
    // If in the adaptive period and an update step
    if(opt_freq != 0 && i < adaptive_period && (i+1)%opt_freq==0){
      if(tempiter(0) > 0 && !no_stats){
	pcur(0) = tempaccepted(0)/tempiter(0);
	if(pcur(0) < popt - (TUNING_ERROR*popt) || pcur(0) > popt + (TUNING_ERROR*popt)){
	  scales(0) = scaletuning2(scales(0),popt,pcur(0));
	  arma::uvec greb;
	  greb = as<arma::uvec>(unfixed[unfixed <= blocks(0,1) & unfixed >= blocks(0,0)]);
	  arma::mat tmp1 = cov(as<arma::mat>(chain));
	  stats_cov_matrix = (1-w)*stats_cov_matrix + w*subset_mat(tmp1,greb);
	  stats_cov_matrix_scaled = stats_cov_matrix*scales(0);
	}
      }
     
      // Update proposal parameters
      for(int iii = 1; iii < no_infections+1;++iii){
	if(tempiter(iii) > 0){
	  pcur(iii) = tempaccepted(iii)/tempiter(iii);
	  if(pcur(iii) < popt - (TUNING_ERROR*popt) || pcur(iii) > popt + (TUNING_ERROR*popt)){
	    scales(iii) = scaletuning2(scales(iii), popt,pcur(iii));
	    arma::uvec greb;
	    greb = as<arma::uvec>(unfixed[unfixed <= (blocks(iii,1)) & unfixed >= (blocks(iii,0))]);
	    arma::mat tmp1 = cov(as<arma::mat>(chain));
	    cov_matrix.slice(iii-1) = (1-w)*cov_matrix.slice(iii-1) + w*subset_mat(tmp1,greb);
	    cov_matrix_scaled.slice(iii-1) = cov_matrix.slice(iii-1)*scales(iii);
	  }
	}
      }
      tempaccepted = clone(reset);
      tempiter = clone(reset);
      pcur = clone(reset);
    }
    
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
  }
  csv_write.close();
  
  return(mcmc_chain_file);
}
