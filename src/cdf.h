#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

// Auxiliar function
#include <algorithm>
#include <iostream>
#include <vector>
// must export since it is used on the quantile routine
Rcpp::NumericVector psep(
  Rcpp::NumericVector
  ,double
  ,double
  ,double
  ,double
);
