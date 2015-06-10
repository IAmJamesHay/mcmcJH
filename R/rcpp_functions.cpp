#include <Rcpp.h>
using namespace Rcpp;

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
    ts = params[6*(i-1)+5] + tp;
    m = params[6*(i-1)+6];

    while(t != times.end() && *t <= final_t){
      tmp = (((*t <= (tp + t_i))*(mu/tp)*(*t-t_i)) + ((*t > (tp + t_i))*(*t <= (ts + t_i))*((-(((1-dRF)*mu)/(ts-tp))*(*t-t_i) + (((1-dRF)*mu)/(ts-tp))*tp + mu))) + ((*t > (ts + t_i))*((m*ts-m*(*t-t_i)) + (dRF*mu)))) + y0;
      if(tmp < trunc_lower){
	tmp = trunc_lower;
      }
      Y[j] = tmp;
      ++t;
      ++j;
    }

    y0 = (((final_t <= (tp + t_i))*(mu/tp)*(final_t-t_i)) + ((final_t > (tp + t_i))*(final_t <= (ts + t_i))*((-(((1-dRF)*mu)/(ts-tp))*(final_t-t_i) + (((1-dRF)*mu)/(ts-tp))*tp + mu))) + ((final_t > (ts + t_i))*((m*ts-m*(final_t-t_i)) + (dRF*mu)))) + y0;

    if(y0 < trunc_lower){
      y0 = trunc_lower;
    }
  }
  
  out(_,0) = times;
  out(_,1)=Y;
  return out;
}

///[[Rcpp::export]]
double mu_prior_norm(NumericVector values, int index, NumericVector args){
  NumericVector yy = dnorm(values,args[1],args[2],args[3]);
  return yy(index-1);
}

