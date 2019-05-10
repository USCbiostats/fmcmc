#include <Rcpp.h>
using namespace Rcpp;

void normal_prop_void(
    NumericVector & ans,
    const NumericVector & x,
    const NumericVector & lb,
    const NumericVector & ub,
    const NumericVector & scale,
    const IntegerVector & fixed
) {
  
  int K = x.size();
  
  // Proposal
  GetRNGstate();
  ans = x + rnorm(x.size())*scale;
  // print(ans);
  PutRNGstate();
  
  for (int k=0; k<K; ++k) {
    
    // Is it fixed?
    if (fixed.at(k)==1) {
      ans.at(k) = x.at(k);
      continue;
    }
    
    // Reflection adjustment
    double d = ub[k] - lb[k];
    while (true) {
      
      if (ans[k] > ub[k]) {
        ans[k] = 2.0*ub[k] - ans[k];
      } else if (ans[k] < lb[k]) {
        ans[k] = 2.0*lb[k] - ans[k];
      } else
        break;
      
    }
    
  }
  
  return;
}

// [[Rcpp::export(rng=false)]]
NumericVector normal_prop(
    const NumericVector & x,
    const NumericVector & lb,
    const NumericVector & ub,
    const NumericVector & scale,
    const IntegerVector & fixed
) {
  
  // Proposal
  NumericVector ans(x.length());
  normal_prop_void(ans, x, lb, ub, scale, fixed);
  
  return ans;
}



