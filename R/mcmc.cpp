#include "mcmc.hpp"
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
double single_strain_model(double mu, double tp, double m, double y0, double ti, double t){
  double y = y0;
  if(ti == -1){
    return(y);
  }
  else {
    
    if(t <= ti){
      y = 0;
    }
    else if(t > ti && t <= (ti+tp)){
      y = (mu/tp)*t - (mu/tp)*ti;
    }
    else {
      y = -m*t + m*(ti+tp) + mu;
    }
    y = y + y0;
  }
  return(y);
}


//[[Rcpp::export]]
NumericMatrix single_infection_model(NumericVector mu_pars, NumericVector tp_pars, NumericVector m_pars, NumericVector y0_pars, int ti, NumericVector t){
  NumericMatrix y(t.size(), mu_pars.size());

  for(int j = 0;j < t.size();++j){
    for(int i =0; i < mu_pars.size();++i){
      y(j,i) = single_strain_model(mu_pars(i),tp_pars(i),m_pars(i),y0_pars(i),ti,t(j));
    }
  }

  return(y);
}

//[[Rcpp::export]]
NumericMatrix matrix_add(NumericMatrix a, NumericMatrix b){
  NumericMatrix c(a.nrow(),a.ncol());
  for(int i=0;i < c.nrow();++i){
    for(int j = 0; j < c.ncol();++j){
      c(i,j) = a(i,j) + b(i,j);
    }
  }
  return(c);
}


//[[Rcpp::export]]
NumericMatrix array_to_matrix(NumericVector a){
  NumericMatrix b(a.size()/3,a.size()/3);

  b(0,0) = a(0);
  b(0,1) = a(1);
  b(0,2) = a(2);
  b(1,0) = a(3);
  b(1,1) = a(4);
  b(1,2) = a(5);
  b(2,0) = a(6);
  b(2,1) = a(7);
  b(2,2) = a(8);

  return(b);
}


//[[Rcpp::export]]
NumericMatrix multiple_infection_model(NumericMatrix mu_matrix, NumericMatrix tp_matrix, NumericMatrix m_matrix, NumericVector y0_pars, NumericVector ti_pars, NumericVector t){
  NumericMatrix y(t.size(),ti_pars.size());
  for(int i = 0;i<ti_pars.size();++i){
    y = matrix_add(y,single_infection_model(mu_matrix(i,_),tp_matrix(i,_),m_matrix(i,_),y0_pars, ti_pars(i),t));
  }

  return(y);
}

double proposal_function(double current, double step){
  double update;
  double move;
  update = R::rnorm(1,current, step*step);
  return(update);
  
}

double posterior(NumericVector params, NumericMatrix t1_titres, NumericMatrix infection_times, NumericMatrix t0_titres, double final_time){
  double sd = 1.0;
  int N = t0_titres.nrow();
  double ln = 0;
  int index = 0;

  NumericVector mu_pars(9);
  NumericVector tp_pars(9);
  NumericVector m_pars(9);
  NumericVector times(1);
  NumericVector single_tmp1(1);
  NumericVector lns;

  NumericMatrix mu_matrix(3,3);
  NumericMatrix tp_matrix(3,3);
  NumericMatrix m_matrix(3,3);
  NumericMatrix tmp(1,3);
  NumericMatrix titres(N,3);

  times[0] = final_time;

  for(int i =0;i <9;++i){
    mu_pars[index++] = params[i];
  }
  mu_matrix = array_to_matrix(mu_pars);

  index = 0;
  for(int i =9;i <18;++i){
    m_pars[index++] = params[i];
  }
  m_matrix = array_to_matrix(m_pars);  

  tp_pars = NumericVector::create(21,21,21,21,21,21,21,21,21);
  tp_matrix = array_to_matrix(tp_pars);

  for(int i=0;i<N;++i){
    tmp = multiple_infection_model(mu_matrix, tp_matrix, m_matrix, t0_titres(i,_),infection_times(i,_),times);
    titres(i,_) = tmp(0,_);
  }
    
  for(int i =0;i < titres.nrow();++i){
    for(int j = 0;j < titres.ncol();++j){
      single_tmp1[0] = titres(i,j);
      lns = dnorm(single_tmp1,t1_titres(i,j),sd,true);
      ln = ln +  std::accumulate(lns.begin(),lns.end(),0.0);
    }
  }

  return(ln);  
}

string run_MCMC(
				  NumericVector startvalue, 
				  int iterations, 
				  NumericMatrix data, 
				  NumericVector step_sizes,
				  double popt, 
				  int opt_freq, 
				  int thin, 
				  int burnin, 
				  int adaptive_period, 
				  string filename,
				  int save_block,
				  bool VERBOSE,
				  NumericMatrix infection_times,
				  NumericMatrix t0_titres,
				  double final_time
				  ){

  ofstream csv_write;
  double TUNING_ERROR = 0.1;
  double newprobab, probab, difflike;
  int no_recorded = 0;
  int sampno = 0;
  string mcmc_chain_file = filename + "_chain.csv";
  NumericMatrix empty_chain(save_block, startvalue.size());
  NumericMatrix chain(save_block, startvalue.size());
  NumericVector current_params(startvalue.size());
  NumericVector sampnos(save_block);
  NumericVector lnlikes(save_block);
  NumericVector tempaccepted(startvalue.size());
  NumericVector tempiter(startvalue.size());
  NumericVector reset(startvalue.size());
  NumericVector pcur(startvalue.size());
  NumericVector temp_transform(startvalue.size());
  
  // Open stream to csv file
  csv_write.open(mcmc_chain_file.c_str());
  
  probab = posterior(startvalue, data, infection_times, t0_titres, final_time);

  // Write first line to csv file
  csv_write << sampno << ",";
  for(NumericVector::iterator it = startvalue.begin();it != startvalue.end();++it){
    csv_write << *it << ",";
  }
  csv_write << probab << endl;
    
  // Print outs if asked
  if(VERBOSE){
    if(opt_freq==0) cout << "Not running adaptive MCMC - opt_freq set to 0" << endl;
    else cout << "Adaptive MCMC -will adapt step size during burnin period" << endl;
  }
  
  // Calculate initial likelihood
  probab = posterior(startvalue, data, infection_times, t0_titres, final_time);
 
  // For each iteration, go through each parameter
  for(int i = 0;i < (iterations + adaptive_period + burnin); ++i){
    for(int j =0; j < current_params.size();++j){

      // Propose new value for current parameter
      proposal = current_params;
      proposal(j) = proposal_function(current_params(j),step_sizes(j));
      
      // Calculate likelihood with this new set of parameters
      newprobab = posterior(proposal, data, infection_times, t0_titres, final_time);

      // Get difference between this and old likelihood
      difflike = newprobab - probab;

      // Metropolis acceptance
      if(R::runif(0,1) < exp(difflike) || difflike > 0){
	current_params = proposal;
	probab = newprobab;
	tempaccepted(j) += 1;
      }
      tempiter(j) += 1; 

      // If an iteration to be saved
      if(sampno%thin ==0){
	// Write sampno, current params and probab to chain
	sampnos(no_recorded) = sampno;
	for(int x = 0; x < current_params.size();++x){
	  chain(no_recorded,x) = current_params(x);
	}
	probabs(no_recorded) = probab;
	no_recorded++;
	
	if(no_recorded >= save_block){
	  // Write chain to csv file
	  for(int x = 0; x < chain.nrow(); ++x){
	    for(int q =0; q < chain.ncol();++q){
	      csv_write << chain(x,q) << ",";
	    }
	    csv_write << endl;
	  }
	  chain = empty_chain;
	  no_recorded = 0
	}
      }
      sampno++;
    }
    if(opt_freq != 0 && i < adaptive_period && i%%opt_freq==0){
      for(int x = 0; x < current_params.size();++x){
	  pcur[x] = tempaccepted(x)/tempiter(x);
	  tempaccepted = reset;
	  tempiter = reset;
	  tmp_transform = step_sizes;
	  if(pcur(x) < popt - (TUNING-ERROR*popt) || pcur(x) > popt + (TUNING_ERROR*popt)){
	    tmp_transform[x] = scaletuning(tmp_transform(x), popt, pcur(x));
	  }
      }
      step_sizes = tmp_transform;
    }
  }
  csv_write.close();
  return(mcmc_chain_file);
}
