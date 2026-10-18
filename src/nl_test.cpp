// Test helpers (internal): the ast quantile table against the exact t quantile.
#include "nl.h"

// z, dz/dV and dz/dnu from the table at the base-scale values v
// [[Rcpp::export]]
NumericMatrix nl_ast_table_cpp(NumericVector v, double nu) {
  AstZTab tab;
  tab.build(nu, true);
  NumericMatrix out(v.size(), 3);
  for (int i = 0; i < v.size(); i++) tab.eval(v[i], out(i, 0), out(i, 1), out(i, 2));
  return out;
}
