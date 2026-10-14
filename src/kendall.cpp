#include <Rcpp.h>
#include "kendall.h"

// Kendall's tau-b of two numeric vectors (no missing values), O(m log m).
// [[Rcpp::export]]
double kendall_tau_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  if (x.size() != y.size()) Rcpp::stop("'x' and 'y' must have the same length");
  return kendall_tau_b(&x[0], &y[0], x.size());
}
