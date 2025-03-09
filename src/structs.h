#ifndef STRUCTS_H
#define STRUCTS_H

struct multimin_params {
  double step_size;
  double tol;
  unsigned maxiter;
  double epsabs;
  double maxsize;
  unsigned method;
};

struct g_params {

  // sample to estimate the parameters/variables (rows)
  Rcpp::NumericVector data;

  // dimension of the problem (number of parameters/variables) (columns)
  size_t n;

  // which transformation to do for each parameter/variable, if any
  Rcpp::IntegerVector type;

  // vector with minimum value for each parameter
  Rcpp::NumericVector xmin;

  // vector with maximum value for each parameter
  Rcpp::NumericVector xmax;

  // depending on the algorithm, one of the following
  // functions is used for the minimization process

  // *f - function to be minimized
  // Arguments, on this order, are
  // vector of data
  // size of vector
  // initial guess pointer - is substituted by the final solution's position
  // parameters to the minimization problem
  // solution - minimum value obtained
  void (*f) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *);

  // *df - derivative of the function to be minimized
  // Arguments, on this order, are
  // vector of data
  // size of vector
  // initial guess pointer - is substituted by the final solution's position
  // parameters to the minimization problem
  // solution - minimum value obtained
  void (* df) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, Rcpp::NumericVector);

  // combination of f and df
  void (* fdf) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *, Rcpp::NumericVector);

  // parameters for the functions
  void *fparams;
};


struct multimin_algorithm {

  const gsl_multimin_fdfminimizer_type *Tfdf;
  const gsl_multimin_fminimizer_type *Tf;
  const char *Tname;
};

// used in the sep quantile function
struct newton_args{
  double x;
  double x_guess;
  Rcpp::NumericVector func_args;
  double step_size;
  double tol;
  int max_iter;
  int verb;
  double (*f)(newton_args) = NULL;
  double (*dfdx)(newton_args) = NULL;
};

#endif /* STRUCTS_H */
