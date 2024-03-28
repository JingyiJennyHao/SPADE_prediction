#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double log_factorial(int x) {
  return lgamma(x + 1);
}

// [[Rcpp::export]]
NumericVector expit(NumericVector k) {
  int n = k.size();
  NumericVector temp(n);
  for (int i = 0; i < n; ++i) {
    temp[i] = 150 * exp(k[i]) / (1 + exp(k[i])) + 0.5;
  }
  return temp;
}

// [[Rcpp::export]]
double li(int ni1,int ni2,int ni3,NumericVector p,NumericVector alpha){
  //Environment ld = Environment::namespace_env("LaplacesDemon");
  //Function ddirichlet = ld["ddirichlet"];
  Function ddirichlet("ddirichlet");
  double temp = log_factorial(ni1 + ni2 + ni3) - log_factorial(ni1) - log_factorial(ni2) - log_factorial(ni3) + ni1 * log(p[0]) + ni2 * log(p[1]) + ni3 * log(p[2]) + as<double>(ddirichlet(Named("x") = p, Named("alpha")=alpha, Named("log", true)));
  return temp;
}

// [[Rcpp::export]]
double llhc(NumericVector t_alpha, NumericVector t_params, DataFrame simulated_dataset) {
  // Create an R function object for the li function
  Function rdirichlet("rdirichlet");
  
  NumericVector alpha = expit(t_alpha);
  NumericVector params = expit(t_params);
  double sum_temp_mean = 0.0;
  
  NumericVector sum1_list = simulated_dataset["sum1_list"];
  NumericVector sum2_list = simulated_dataset["sum2_list"];
  NumericVector sum3_list = simulated_dataset["sum3_list"];
  
  for (int i = 0; i < simulated_dataset.nrows(); ++i) {
    int ni1 = sum1_list[i];
    int ni2 = sum2_list[i];
    int ni3 = sum3_list[i];
    
    NumericMatrix pis = rdirichlet(1000, params);
    double temp = 0.0;
    
    for (int j=0; j<pis.nrow();++j){
      NumericVector row = pis(j,_);
      double log_likelihood = li(ni1, ni2, ni3, row, alpha);
      temp += log_likelihood;
      }
    
    double temp_mean = temp/1000;
    sum_temp_mean += temp_mean;
  }
  
  return -sum_temp_mean;
}
