


struct multimin_params {
  double step_size;
  double tol;
  unsigned maxiter;
  double epsabs;
  double maxsize;
  unsigned method;
  unsigned verbosity;
};

struct g_params {

  // sample to estimate the parameters/variables (rows)
  std::vector<double> data;

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
  void (*f) (std::vector<double>&, const size_t, std::vector<double>&, void *, double *);

  // *df - derivative of the function to be minimized
  // Arguments, on this order, are
  // vector of data
  // size of vector
  // initial guess pointer - is substituted by the final solution's position
  // parameters to the minimization problem
  // solution - minimum value obtained
  void (* df) (std::vector<double>&, const size_t, std::vector<double>&, void *, std::vector<double>&);

  // combination of f and df
  void (* fdf) (std::vector<double>&, const size_t, std::vector<double>&, void *, double *, std::vector<double>&);

  // parameters for the functions
  void *fparams;
};


struct multimin_algorithm {

  const gsl_multimin_fdfminimizer_type *Tfdf;
  const gsl_multimin_fminimizer_type *Tf;
  const char *Tname;
};
