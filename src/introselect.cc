#include <vector>
#include <algorithm>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double introselect(NumericVector data_, int k) {
  // copy the vector since nth_element modifies the vector in place
  // and we'd rather not have side-effects
  vector<double> data = as< vector<double> >(data_);

  // Map R's 1-based indexing to C++'s 0-based indexing
  k = k-1;
  
  nth_element(data.begin(), data.begin() + k, data.end());
  return data[k];
}
