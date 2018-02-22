#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector normal_prop(
    const NumericVector & x,
    const NumericVector & lb,
    const NumericVector & ub,
    const NumericVector & scale,
    const LogicalVector & fixed
) {
  
  int K = x.size();
  
  // Proposal
  NumericVector ans(K); // = x + rnorm(K, 0, 1)*scale;
  
  for (int k=0; k<K; k++) {
    
    // Is it fixed?
    if (fixed.at(k)) {
      ans.at(k) = x.at(k);
      continue;
    }
      
    // Proposal
    ans.at(k) = x.at(k) + norm_rand()*scale.at(k);
    
    
    // Reflection adjustment
    while( (ans[k] > ub[k]) | (ans[k] < lb[k]) ) {
      
      if (ans[k] > ub[k]) {
        ans[k] = 2.0*ub[k] - ans[k];
      } else {
        ans[k] = 2.0*lb[k] - ans[k];
      }  
      
    }
    
  }
  
  
  return ans;
}

// [[Rcpp::export]]
NumericMatrix MCMCcpp(
    Function & fun,
    NumericVector theta0,
    int nbatch,
    const NumericVector & lb,
    const NumericVector & ub,
    const NumericVector & scale,
    const LogicalVector & fixed
) {
  
  int K = lb.size();
  
  NumericMatrix ans(nbatch, K);
  NumericVector theta1(K);
  NumericVector f0 = fun(theta0), f1;
  
  // Checking values
  if (is_na(f0)[0] || is_nan(f0)[0])
    stop("fun(par) is undefined. Check either -fun- or the -lb- and -ub- parameters.");
  
  // Using sugar to generate the random values for the hastings ratio.
  NumericVector R = runif(nbatch);

  for (int i = 0; i < nbatch; i++) {
  
    // Generating proposal
    theta1 = normal_prop(theta0, lb, ub, scale, fixed);
    f1     = fun(theta1);
  
    // Checking values
    if (is_na(f1)[0] || is_nan(f1)[0])
      stop("fun(par) is undefined. Check either -fun- or the -lb- and -ub- parameters.");
    
    // Metropolis-Hastings ratio
    if (R.at(i) < fmin(1, exp( f1.at(0) - f0.at(0) ))) {
      
      theta0 = theta1;
      f0     = f1;
      
    }
    
    // Storing the current state
    ans(i,_) = theta0;
    
  }
  
  return ans;
}