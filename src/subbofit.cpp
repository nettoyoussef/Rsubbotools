/*
  subbofit (ver. 1.2.1) -- Fit a power exponential density via maximum likelihood

  Copyright (C) 2007-2014 Giulio Bottazzi
  Copyright (C) 2020-2021 Elias Youssef Haddad Netto

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  (version 2) as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/



//Main headers

#include "common.h"

/*------------------*/


/* Output Functions */
/*----------------- */

double geta(Rcpp::NumericVector& data, const double b, const double mu){

  double sum = 0.0;
  unsigned utmp1;
  int size = data.size();

  for(utmp1 = 0;utmp1<size;utmp1++){
    sum+=pow(fabs(mu-data[utmp1]),b);
  }

  return(pow(sum/size, 1./b));
}

/*----------------- */
/*
  par   array of paramaters
  N     number of observations
  dim   dimension of the matrix: 2 m known; 3 m unknown
  I     the variance-covariance matrix
*/
RcppGSL::Matrix subbo_varcovar(const Rcpp::NumericVector par, const size_t N, const size_t dim){

  RcppGSL::Matrix I(dim, dim);
  const double b = par[0];
  const double a = par[1];
  //const double m = par[2]; // declared but not used in the original
  const double dtmp1= log(b)+gsl_sf_psi(1+1./b);

  size_t i, j;

  /* needed for inversion */
  RcppGSL::Matrix J(dim, dim);
  //gsl_matrix *J = gsl_matrix_calloc (dim, dim);
  gsl_permutation *P = gsl_permutation_alloc (dim);
  int signum;


  /* define the matrix */

  /* b - b */
  I(0,0) = (dtmp1*dtmp1+gsl_sf_psi_1(1+1/b)*(1+1/b)-1)/pow(b, 3);

  /* b - a */
  I(0,1) = -dtmp1/(a*b);
  I(1,0) = -dtmp1/(a*b);

  /* a - a */
  I(1,1) = b/(a*a);


  if(dim == 3){

    /* b - m */
    I(0,2) = 0;
    I(2,0) = 0;

    /* a - m */
    I(1,2) = 0;
    I(2,1) = 0;

    /* m - m */
    I(2,2) = (pow(b,-2/b+1)*gsl_sf_gamma(2-1/b))/
                    (gsl_sf_gamma(1+1/b)*a*a);
  }

  // print matrix
  //for(i = 0;i<dim;i++){
  //  for(j = 0;j<dim;j++){
  //    Rprintf("%.4e ", gsl_matrix_get(I,i,j));
  //  }
  //  Rprintf("\n");
  //}
  //Rprintf("-------------\n");

  /* invert I; in J store temporary LU decomp. */
  gsl_matrix_memcpy (J, I);
  gsl_linalg_LU_decomp (J, P,&signum);
  gsl_linalg_LU_invert (J, P,I);

  /* free allocated memory */
  gsl_permutation_free(P);
  /* freeing the below matrix gets an error - double freeing*/
  /*gsl_matrix_free(J);*/

  /* set the var-covar matrix */
  for(i = 0;i<dim;i++){
    for(j = 0;j<i;j++){
      I(i,j) = I(i,j)/sqrt(I(i,i)* I(j,j));
    }
  }

  for(i = 0;i<dim;i++){
    for(j = i;j<dim;j++){
      I(i,j) = I(i,j)/N;
    }
  }

  return I;
}


/* Method of Moments */
/*------------------ */


/*

The b parameter is determined by the equation:

sqrt( gamma[3/b]*gamma[1/b] )/gamma[2/b] = stdev/adev

taking the logarithm l = log(stdev/adev) and with x= -log(b)
one obtains the equation:

.5*lgamma(3*exp(x))+.5*lgamma(exp(x))-lgamma(2*exp(x))-l ==0

*/

/* this store the stdev on the average absolute deviation */
struct mm_params {
  double l;  /* l = log(stdev/adev) */
};


// function for calculating b from the method
// of moments
/* x = -log(b) */
double mm_f(double x, void *params){

  struct mm_params *p= (struct mm_params *) params;
  const double dtmp1 = exp(x);
  const double y = .5*gsl_sf_lngamma(3.*dtmp1)+.5*gsl_sf_lngamma(dtmp1)-
    gsl_sf_lngamma(2.*dtmp1)-(p->l);

  return(y);
}

// first derivative of the function for the method of
// moments
double mm_df(double x, void *params){

  const double dtmp1 = exp(x);
  const double dy = dtmp1*(1.5*gsl_sf_psi(3.*dtmp1)+.5*gsl_sf_psi(dtmp1)-
                           2.*gsl_sf_psi(2.*dtmp1));

  return(dy);
}


void mm_fdf(double x, void *params, double *f, double *df){


  struct mm_params *p= (struct mm_params *) params;
  const double dtmp1 = exp(x);

  *f = .5*gsl_sf_lngamma(3.*dtmp1)+.5*gsl_sf_lngamma(dtmp1)-
    gsl_sf_lngamma(2.*dtmp1)-(p->l);

  *df = dtmp1*(1.5*gsl_sf_psi(3.*dtmp1)+.5*gsl_sf_psi(dtmp1)-
             2.*gsl_sf_psi(2.*dtmp1));
}



Rcpp::List moment(Rcpp::NumericVector& data){

  double average, aad, sdev, var;
  unsigned j;
  double s = 0.0, p = 0.0;
  int size = data.size();


  for (j = 0;j<size;j++){
    s += data[j];
  }

  average = s/size;
  aad=(var) = 0.0;

  for (j = 0;j<size;j++) {

    // calculates mean (or average) absolute deviation
    aad += fabs(s = data[j]-(average));

    // calculates variance
    var += (p = s*s);
  }


  aad /= size;
  var = (var)/(size-1);

  // calculates standard deviation
  sdev = sqrt(var);

  Rcpp::List ans =
    Rcpp::List::create(
                        Rcpp::Named("average") = average
                       ,Rcpp::Named("aad")     = aad
                       ,Rcpp::Named("sdev")    = sdev
                       ,Rcpp::Named("var")     = var
                       );

  return ans;
}


/* determination of b parameter from the ratio stdev/ave */
/* returns the parameter b */
// [[Rcpp::export]]
double mm_subbotin(const double std_over_aad, int verb = 0){

  /* default initial value for -log(b)*/
  double mlogb = 0.0;

  // this solver find roots using methods
  // that require derivatives
  // documentation available at
  // https://www.gnu.org/software/gsl/doc/html/roots.html
  gsl_root_fdfsolver *d1solver;
  gsl_function_fdf FdF;

  struct mm_params params;

  // number of iterations
  unsigned iter;
  double dtmp1;
  int status;
  int max_iter = 500;

  /* initialize the parameter */
  params.l = log(std_over_aad);

  /* allocate the solver */
  // From the documentation
  // The Steffenson Method provides the fastest convergence of all the routines. It
  // combines the basic Newton algorithm with an Aitken “delta-squared” acceleration.
  // If the Newton iterates are x_i then the acceleration procedure generates a new
  // sequence which converges faster than the original sequence under reasonable
  // conditions.
  // Notice: as all acceleration procedures, it can become unstable if the function is
  // not well-behaved
  d1solver = gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_steffenson);
  //d1solver = gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_newton);

  /* set the function to solve */
  FdF.f = &mm_f;       // function to solve the root finding problem
  FdF.df = &mm_df;     // first derivative

  // this is a function that receives a pointer from the two previous functions
  // as, apparently, it is faster to calculate both values numerically at
  // the same time
  FdF.fdf = &mm_fdf;

  // this function receives the parameters from the function to be
  // minimized
  FdF.params = &params;

  /* initial point */
  // initial value for mlogb = 0
  gsl_root_fdfsolver_set (d1solver,&FdF, mlogb);

  // number of iterations
  iter = 0;

  do{
    iter++;
    status = gsl_root_fdfsolver_iterate (d1solver);

    if(verb > 1){
      Rprintf("# iteration: %d; b=%f; lb=%f:\n",
              iter, exp(-mlogb), -mlogb);
    }

    if(status){
      Rprintf("# WARNING in 1d solver: %s\n", gsl_strerror (status));
      Rprintf("# the error 'problem with user-supplied function; means that:\n");
      Rprintf("# 1. the function value is not finite (bad convergence):\n");
      Rprintf("# 2. the function derivative is not finite (bad convergence):\n");

      return(exp(-mlogb));
    }

    // salves the result from the previous iteration
    dtmp1 = mlogb;

    // updates mlogb with the results from the current iteration
    mlogb = gsl_root_fdfsolver_root (d1solver);

    // This function tests for the convergence of the sequence
    // dtmp1, mlogb  with absolute error 0.0 and relative error 1e-4
    status = gsl_root_test_delta (dtmp1, mlogb, 0.0, 1e-4);

    if(status == GSL_SUCCESS)
      if(verb > 1){
        Rprintf("# Converged after %d iterations lb=%f:\n",
                iter,-mlogb);
      }
    if(iter>=max_iter){
      Rprintf("# WARNING in 1d solver  : exceeded max. num. of iterations %d\n"
              ,max_iter);
      break;
    }
  }
  while (status == GSL_CONTINUE);
  gsl_root_fdfsolver_free (d1solver);

  return(exp(-mlogb));
}


/*------------------ */



/* Objective Functions */
/* ================ */


/* Functions of a and m */
/* -------------------- */

/* reduced log likelyhood x[0] = b x[1] = mu */
void subbo_objf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f){

  size_t size = data.size();
  const double b = x[0]; // default to 2
  const double mu = x[1]; // default to 0

  double sum = 0.0;
  int i = 0;
  int status;
  gsl_sf_result result;


  for(i = 0; i<size; ++i){

    sum+=pow(fabs(mu-data[i]),b);
  }
  // error message
  if( (status = gsl_sf_lngamma_e(1./b+1.,&result)) ){
    Rprintf("subbo_objf: lngamma(%e)=%e status  = %s \n", 1./b+1, result.val, gsl_strerror(status));
    Rprintf("b=%e mu=%e\n",b, mu);
  }

  *f = log(2.)+log(sum/size)/b + log(b)/b+result.val+1./b;
}


/* derivative of the reduced log likelyhood x[0] = b x[1] = mu */
void subbo_objdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, Rcpp::NumericVector df){

  size_t size = data.size();
  const double b = x[0];
  const double mu = x[1];

   double sum1 = 0.0;
   double sum2 = 0.0;
   double sum3 = 0.0;
   unsigned utmp1;

   int status;
   gsl_sf_result result;

   for(utmp1 = 0;utmp1<size;utmp1++){
     const double dtmp1 = fabs(mu-data[utmp1]);
     sum1+=pow(dtmp1, b);
     if(dtmp1!=0.0){
       sum2+=pow(dtmp1, b-1.)*(mu>data[utmp1]?1.:-1.);
       sum3+=pow(dtmp1, b)*log(dtmp1);
     }
   }

   status = gsl_sf_psi_e(1./b+1.,&result);

   if(status){
     Rprintf("subbo_objdf [psi] status  = %s\n", gsl_strerror(status));
     Rprintf("b=%e mu=%e\n",b, mu);
     Rcpp::stop("Error.");
   }

   df[0] = sum3/(b*sum1)-(log(sum1/size)+log(b)+result.val)/(b*b);
  df[1] = sum2/sum1;

}

/* reduced likelyhood and derivatives x[0] = b x[1] = mu */
void subbo_objfdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f, Rcpp::NumericVector df){

  int size = data.size();
  const double b = x[0];
  const double mu = x[1];

  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  unsigned utmp1;

  int status;
  gsl_sf_result result;

  for(utmp1 = 0;utmp1<size;utmp1++){
    const double dtmp1 = fabs(mu-data[utmp1]);
    sum1+=pow(dtmp1, b);
    if(dtmp1!=0.0){
      sum2+=pow(dtmp1, b-1.)*(mu>data[utmp1]?1.:-1.);
      sum3+=pow(dtmp1, b)*log(dtmp1);
    }
  }


  if( (status = gsl_sf_lngamma_e(1./b+1.,&result)) ){
    Rprintf("subbo_objfdf [lngamma] status  = %s\n", gsl_strerror(status));
    Rprintf("b=%e mu=%e\n",b, mu);
  }


  *f = log(2.)+log(sum1/size)/b + log(b)/b+result.val+1./b;

  if( (status = gsl_sf_psi_e(1./b+1.,&result)) ){
    Rprintf("subbo_objfdf [psi] status  = %s\n", gsl_strerror(status));
    Rprintf("b=%e mu=%e\n",b, mu);
    Rcpp::stop("Error.");
  }

  df[0] = sum3/(b*sum1)-(log(sum1/size)+log(b)+result.val)/(b*b);
  df[1] = sum2/sum1;

}


// [[Rcpp::export]]
Rcpp::List optim_method_moments(
                                Rcpp::NumericVector data
                                ,double fmin
                                ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
                                ,int verb = 0
                                ){

  // variables
  Rcpp::NumericVector par(3);
  Rcpp::NumericVector x(2);
  double ave, adev, sdev, var;

  /* compute statistics */
  Rcpp::List moments =  moment(data);
  ave  = moments[0];
  adev = moments[1];
  sdev = moments[2];
  var  = moments[3];


  /* method of moments estimation for mu */
  if(provided_m_.isNotNull()){

    // casting the null value
    // necessary according to
    // https://gallery.rcpp.org/articles/optional-null-function-arguments/
    Rcpp::NumericVector provided_m(provided_m_);
    par[2] = provided_m[0];
  }
  else{
    par[2] = ave;
  }

  if(verb > 1){
    Rprintf("METHOD OF MOMENTS\n");
    Rprintf("# Moments of the data\n");
    Rprintf("#  par    mean  var   sd    aad\n");
    Rprintf("#  value  %.2f  %.2f  %.2f  %.2f\n", ave, var, sdev, adev);
    Rprintf("\n");
  }

  /* method of moments estimation for b */
  // this finds the root of
  par[0] = mm_subbotin(sdev/adev, verb);

  /*  method of moment estimation for a */
  par[1] = sdev*pow(par[0],-1./par[0])*
    exp(.5*(gsl_sf_lngamma(1./par[0])-gsl_sf_lngamma(3./par[0])));

  /* reduced log likelyhood x[0] = b x[1] = mu */
  x[0]= par[0];
  x[1]= par[2];
  subbo_objf(data, 2, x, NULL, &fmin);

  if(verb > 1){
    Rprintf("Method of moments results\n");
    Rprintf("#  par    b     a     m     ll\n");
    Rprintf("#  value  %.2f  %.2f  %.2f  %.2f\n", par[0], par[1], par[2], fmin);
    Rprintf("\n");
    Rprintf("END OF METHOD OF MOMENTS\n");
  }

  Rcpp::List ans =
    Rcpp::List::create(
                        Rcpp::Named("par") = par
                       ,Rcpp::Named("fmin") = fmin
                       );
  return ans;
}



//' Fit a power exponential density via maximum likelihood
//'
//' \code{subbofit} returns the parameters, standard errors. negative
//' log-likelihood and covariance matrix of the Subbotin Distribution for a
//' sample. The process can execute three steps, dependending on the level of
//' accuracy required. See details below.
//'
//'
//' The Subbotin distribution is a exponential power distribution controlled
//' by three parameters, with formula:
//' \deqn{f(x;a,b,m) = \frac{1}{A} e^{-\frac{1}{b} |\frac{x-m}{a}|^b}}
//' with:
//' \deqn{A = 2ab^{1/b}\Gamma(1+1/b)}
//' where \eqn{a} is a scale parameter, \eqn{b} controls the tails (lower values
//' represent fatter tails), and \eqn{m} is a location parameter. Due to its
//' simmetry, the equations are simple enough to be estimated by the method of
//' moments, which produce rough estimations that should be used only for first
//' explorations. The maximum likelihood global estimation improves on this
//' initial guess by using a optimization routine, defaulting to the
//' Broyden-Fletcher-Goldfarb-Shanno method. However, due to the lack of
//' smoothness of this function on the \eqn{m} parameter (derivatives are zero
//' whenever \eqn{m} equals a sample observation), an exhaustive search must be
//' done by redoing the previous step in all intervals between two observations.
//' For a sample of \eqn{n} observations, this would lead to \eqn{n-1}
//' optimization problems. Given the computational cost of such procedure,
//' an interval search is used, where the optimization is repeated in the
//' intervals at most the value of the \emph{interv_step} from the last
//' minimum found. Details on this method are available on the package vignette.
//'
//' @param data (NumericVector) - the sample used to fit the distribution.
//' @param verb (int) - the level of verbosity. Select one of:
//' * 0  just the final result
//' * 1  headings and summary table
//' * 2  intermediate steps results
//' * 3  intermediate steps internals
//' * 4+  details of optim. routine
//' @param method int - the steps that should be used to estimate the
//' parameters.
//' * 0 no optimization perform - just return the log-likelihood from initial guess.
//' * 1 initial estimation based on method of moments
//' * 2 global optimization not considering lack of smoothness in m
//' * 3 interval optimization taking non-smoothness in m into consideration
//' @param interv_step  int - the number of intervals to be explored after
//' the last minimum was found in the interval optimization. Default is 10.
//' @param provided_m_ NumericVector - if NULL, the m parameter is estimated
//' by the routine. If numeric, the estimation fixes m to the given value.
//' @param par NumericVector - vector containing the initial guess for
//' parameters b, a and m, respectively. Default values of are c(2, 1, 0).
//' @param g_opt_par NumericVector - vector containing the global optimization
//' parameters.
//' The optimization parameters are:
//' * step  - (num) initial step size of the searching algorithm.
//' * tol   - (num) line search tolerance.
//' * iter  - (int) maximum number of iterations.
//' * eps   - (num) gradient tolerance. The stopping criteria is \eqn{||\text{gradient}||<\text{eps}}.
//' * msize - (num) simplex max size. stopping criteria given by \eqn{||\text{max edge}||<\text{msize}}
//' * algo  - (int) algorithm. the optimization method used:
//'   * 0 Fletcher-Reeves
//'   * 1 Polak-Ribiere
//'   * 2 Broyden-Fletcher-Goldfarb-Shanno
//'   * 3 Steepest descent
//'   * 4 Nelder-Mead simplex
//'   * 5 Broyden-Fletcher-Goldfarb-Shanno ver.2
//'
//' Details for each algorithm are available on the [GSL Manual](https://www.gnu.org/software/gsl/doc/html/multimin.html).
//' Default values are c(.1, 1e-2, 100, 1e-3, 1e-5, 3,0).
//' @param itv_opt_par NumericVector - interval optimization parameters. Fields
//' are the same as the ones for the global optimization. Default values
//' are c(.01, 1e-3, 200, 1e-3, 1e-5, 5, 0).
//' @return a list containing the following items:
//' * "dt" - dataset containing parameters estimations and standard deviations.
//' * "log-likelihood" - negative log-likelihood value.
//' * "matrix" - the covariance matrix for the parameters.
//'
//' @examples
//' sample_subbo <- rpower(1000, 1, 2)
//' subbofit(sample_subbo)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List subbofit(
              Rcpp::NumericVector data
             ,int verb = 0
             ,int method = 3
             ,int interv_step = 10
             ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
             ,Rcpp::NumericVector par = Rcpp::NumericVector::create(2.,1.,0.)
             ,Rcpp::NumericVector g_opt_par = Rcpp::NumericVector::create(.1, 1e-2, 100, 1e-3, 1e-5, 3)
             ,Rcpp::NumericVector itv_opt_par = Rcpp::NumericVector::create(.01, 1e-3, 200, 1e-3, 1e-5, 5)
                    ){

  // check arguments
  int check_par_size = par.size();
  if(check_par_size != 3){
    Rcpp::stop("Number of default values for parameters must be equal to three.");
  }

  /* initial values */
  /* -------------- */

  /* store possibly provided values for parameters */
  unsigned is_m_provided = 0;

  // define optimization parameters
  struct multimin_params global_oparams =
    {(double)       g_opt_par[0]
     ,(double)       g_opt_par[1]
     ,(unsigned int) g_opt_par[2]
     ,(double)       g_opt_par[3]
     ,(double)       g_opt_par[4]
     ,(unsigned int) g_opt_par[5]
    };

  struct multimin_params interv_oparams =
    { (double)       itv_opt_par[0]
     ,(double)       itv_opt_par[1]
     ,(unsigned int) itv_opt_par[2]
     ,(double)       itv_opt_par[3]
     ,(double)       itv_opt_par[4]
     ,(unsigned int) itv_opt_par[5]
  };

  // casting the null value
  // necessary according to
  // https://gallery.rcpp.org/articles/optional-null-function-arguments/

  // initial value for m parameter
  if(provided_m_.isNotNull()){

    is_m_provided = 1;

  }

  /* store data */
  unsigned int Size = data.size();/*the number of data*/

  // Fisher Matrix or negative log-likelihood?
  // log_likelihood
  double fmin = 0;

  /* customized admin. of errors */
  /* --------------------------- */
  gsl_set_error_handler_off ();


  /* sort data */
  /* --------- */
  sortRcpp(data);
  //qsort(data, Size, sizeof(double),sort_by_value);

  /* reduced log likelyhood x[0] = b x[1] = mu */
  Rcpp::NumericVector x = Rcpp::NumericVector::create(par[0], par[2]);
  subbo_objf(data, 2, x, NULL, &fmin); // this updates the value of the reduced log likelyhood in fmin

  if(verb > 0){

    /* output of initial values */
    /* ------------------------ */
    Rprintf("INITIAL VALUES\n");
    Rprintf("#  par    b     a     m     ll\n");
    Rprintf("#  value  %.2f  %.2f  %.2f  %.2f\n", par[0], par[1], par[2], fmin);
    Rprintf("\n");

  }

  /* no minimization */
  /* --------------- */
  if(method == 0){

    // Name of parameters
    Rcpp::CharacterVector param_names =
      Rcpp::CharacterVector::create("b", "a", "m");

    // main dataframe with coefficients
      Rcpp::DataFrame dt =
        Rcpp::DataFrame::create( Rcpp::Named("param") = param_names
                                 ,Rcpp::Named("coef") = par
                                 );

      Rcpp::List ans = Rcpp::List::create(
                                          Rcpp::Named("dt")              = dt
                                          ,Rcpp::Named("log-likelihood") = fmin
                                          );

      return ans;

  }

  /* method of moments estimation */
  /* ---------------------------- */
  if(method >= 1){

    // calculate parameters with the method of moments
    Rcpp::List moment_results  = optim_method_moments(data
                                                      ,fmin
                                                      ,provided_m_
                                                      ,verb
                                                      );
    par = moment_results["par"];
    fmin = moment_results["fmin"];

  }

  /* global optimization  */
  /* -------------------- */
  if(method >= 2){

    x = Rcpp::NumericVector::create(par[0], par[2]);

    // notice: x is updated by reference
    Rcpp::List g_opt_results =
      global_optim(
                   data
                   ,fmin
                   ,global_oparams
                   ,x
                   ,(unsigned) 2
                   ,(unsigned) 1
                   ,subbo_objf
                   ,subbo_objdf
                   ,subbo_objfdf
                   ,provided_m_
                   ,verb
                   );

    /* store final values */
    /* ------------------ */
    par[0] = x[0]; //sets b to the new minimum
    par[2] = x[1]; //sets m to the new minimum
    par[1] = geta(data, par[0], par[2]); // updates 'a' value

    /* interval optimization  */
    /* ---------------------- */

    if(method >= 3){

      g_opt_results =
        interval_optim(
                       data
                       ,g_opt_results["type"]
                       ,g_opt_results["xmin"]
                       ,g_opt_results["xmax"]
                       ,x
                       ,g_opt_results["fmin"]
                       ,interv_oparams /* interval optimization parameters */
                       ,interv_step /* interval optimization step */
                       ,(unsigned) 2
                       ,(unsigned) 1
                       ,subbo_objf
                       ,subbo_objdf
                       ,subbo_objfdf
                       ,verb
                       );

      /* store final values */
      /* ------------------ */
      par[0] = x[0]; //sets b to the new minimum
      par[2] = x[1]; //sets m to the new minimum
      par[1] = geta(data, par[0], par[2]); // updates 'a' value

      // updates log-likelihood
      fmin = Rcpp::as<double>(g_opt_results["fmin"]);
    }
  }

  // generate outputs
  // variables
  const size_t dim=(is_m_provided?2:3); /* set the size of the var-covar matrix */

  /* allocate var-covar matrix */
  // The V matrix has on its main diagonal the variance of parameters
  // on the lower diagonal the correlation coefficients
  // on the upper diagonal the covariances
  RcppGSL::Matrix V =  subbo_varcovar(par, Size, dim);
  // this matrix in its upper diagonal presents the covariances
  // and on its lower diagonal presents the correlation coefficients between the parameters

  // vector of standard errors
  Rcpp::NumericVector std_error =
    Rcpp::NumericVector::create(
                                 sqrt(V(0,0)) // b
                                ,sqrt(V(1,1)) // a
                                ,sqrt(V(2,2)) // m
                                );

   // Name of parameters
   Rcpp::CharacterVector param_names = Rcpp::CharacterVector::create("b", "a", "m");

   // main dataframe with coefficients
   Rcpp::List dt =
     Rcpp::DataFrame::create( Rcpp::Named("param")     = param_names
                             ,Rcpp::Named("coef")      = par
                             ,Rcpp::Named("std_error") = std_error
                       );

  // convert matrix
  Rcpp::NumericMatrix matrix = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(V));

  // add names for the matrices
  colnames(matrix) = param_names;

  Rcpp::List ans = Rcpp::List::create(
                                       Rcpp::Named("dt")             = dt
                                      ,Rcpp::Named("log-likelihood") = fmin
                                      ,Rcpp::Named("matrix")         = matrix
                                     );


  return ans;

}
