
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


void subbo_printcumul(Rcpp::NumericVector data, double param[]){

  int size = data.size();
  unsigned i;
  double dtmp1; 
  const double b = param[0];
  const double a = param[1];
  const double m = param[2];

  for(i = 0;i<size;i++){
    dtmp1 = pow(fabs(data[i]-m)/a, b)/b;
    dtmp1= .5+.5*(data[i]>m?1.:-1.)*gsl_sf_gamma_inc_P(1./b, dtmp1);
    Rprintf("%e %e\n",data[i], dtmp1);
  }

}

void subbo_printdensity(Rcpp::NumericVector data, double param[]){

  int size = data.size();
  unsigned i;
  const double b = param[0];
  const double a = param[1];
  const double m = param[2];

  const double norm = 2*a*gsl_sf_gamma(1./b+1.)*pow(b, 1./b);

  for(i = 0;i<size;i++){
    const double dtmp1 = data[i];
    Rprintf("%e ",dtmp1);
    Rprintf("%e\n",exp(-pow(fabs(dtmp1-m)/a, b)/b)/norm);
  }

}
/*----------------- */

/* 
   par   array of paramaters
   N     number of observations
   dim   dimension of the matrix: 2 m known; 3 m unknown
   I     the variance-covariance matrix
*/
// [[Rcpp::export]]
RcppGSL::Matrix varcovar(const Rcpp::NumericVector par, const size_t N, const size_t dim){

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

  /* invert I; in J store temporary LU decomp. */
  gsl_matrix_memcpy (J, I);
  gsl_linalg_LU_decomp (J, P,&signum);
  gsl_linalg_LU_invert (J, P,I);

  /* free allocated memory */
  gsl_permutation_free(P);


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
double mm(const double std_over_aad){

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

    Rprintf("# iteration: %d; b=%f; lb=%f:\n",
                iter, exp(-mlogb), -mlogb);
    
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
        Rprintf("# Converged after %d iterations lb=%f:\n",
                iter,-mlogb);

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

  int size = data.size();
  const double b = x[0]; // default to 2
  const double mu = x[1]; // default to 0

  double sum = 0.0;
  int i = 0;
  int status;
  gsl_sf_result result;


  for(i = 0; i<size; ++i){

    sum+=pow(fabs(mu-data[i]),b);
  }
  if( (status = gsl_sf_lngamma_e(1./b+1.,&result)) ){
    Rprintf("subbo_objf: lngamma(%e)=%e status  = %s \n", 1./b+1, result.val, gsl_strerror(status));
    Rprintf("b=%e mu=%e\n",b, mu);
  }

  *f = log(2.)+log(sum/size)/b + log(b)/b+result.val+1./b;
}


/* derivative of the reduced log likelyhood x[0] = b x[1] = mu */
void subbo_objdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, Rcpp::NumericVector df){

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


Rcpp::NumericVector optim_method_moments(
                                          Rcpp::NumericVector data
					 ,double *fmin  
	                                 ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
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
  
  Rprintf("METHOD OF MOMENTS\n");

  Rprintf("# Moments of the data\n");
  Rprintf("#  par    mean  var   sd    aad\n");
  Rprintf("#  value  %.2f  %.2f  %.2f  %.2f\n", ave, var, sdev, adev);
  Rprintf("\n");
  

  /* method of moments estimation for b */
  // this finds the root of 
  par[0] = mm(sdev/adev);
  
  /*  method of moment estimation for a */
  par[1] = sdev*pow(par[0],-1./par[0])*
    exp(.5*(gsl_sf_lngamma(1./par[0])-gsl_sf_lngamma(3./par[0])));

  /* reduced log likelyhood x[0] = b x[1] = mu */
  x[0]= par[0];
  x[1]= par[2];
  subbo_objf(data, 2, x, NULL, &fmin);

  Rprintf("Method of moments results\n");
  Rprintf("#  par    b     a     m     ll\n");
  Rprintf("#  value  %.2f  %.2f  %.2f  %.2f\n", par[0], par[1], par[2], fmin);
  Rprintf("\n");
    

  Rprintf("END OF METHOD OF MOMENTS\n");

  return par;
}

  Rcpp::List ans =
    Rcpp::List::create(
		        Rcpp::Named("par") = par
		       ,Rcpp::Named("fmin") = fmin
		       );
  return ans;
}


//' subbofit
//' Fit a power exponential density via maximum likelihood
//' 
//'   Verbosity levels:
//'   0 just the final ouput
//'   1 the results of intermediate steps
//'   2 internal information on intermediate steps
//'   3 gory details on intermediate steps
//' 
//' Fit symmetric power exponential density. Read from files of from standard input 
//' Usage: %s [options] [files]                                                     
//' Options:                                                                       
//'  -O  output type (default 0)                
//'       0  parameter b a m and log-likelihood   
//'       1  the estimated distribution function computed on the provided points     
//'       2  the estimated density function computed on the provided points 
//'       3  parameters b a m and their standard errors 
//'  -x  set initial conditions b, a,m  (default 2, 1,0)
//'  -m  the mode is not estimated but is set to the value provided
//'  -s  number of intervals to explore at each iteration (default 10)
//'  -V  verbosity level (default 0)           
//'       0  just the final result        
//'       1  headings and summary table   
//'       2  intermediate steps results   
//'       3  intermediate steps internals 
//'       4+  details of optim. routine   
//'  -M  active estimation steps. The value is the sum of (default 7)
//'       1  initial estimation based on method of moments
//'       2  global optimization not considering lack of smoothness in m
//'       4  local optimization taking non-smoothness in m into consideration 
//'  -G  set global optimization options. Fields are step, tol, iter, eps, msize, algo.
//'      Empty field implies default (default .1, 1e-2, 100, 1e-3, 1e-5, 3)
//'  -I  set local optimization options. Fields are step, tol, iter, eps, msize, algo.
//'      Empty field implies default (default .01, 1e-4, 200, 1e-4, 1e-5, 5)
//' The optimization parameters are
//'  step  initial step size of the searching algorithm                  
//'  tol  line search tolerance iter: maximum number of iterations      
//'  eps  gradient tolerance : stopping criteria ||gradient||<eps       
//'  msize  simplex max size : stopping criteria ||max edge||<msize     
//'  algo  optimization methods:
//'          0 Fletcher-Reeves
//'          1 Polak-Ribiere     
//'          2 Broyden-Fletcher-Goldfarb-Shanno
//'          3 Steepest descent           
//'          4 Nelder-Mead simplex
//'          5 Broyden-Fletcher-Goldfarb-Shanno ver.2
//' 
//' Examples:
//'  'subbofit -m 1 -M 6 <file'  estimate a and b with m = 1 and skipping initial  
//'                              method of moments estimation
//'
//' methods of estimation
//' 0 - no minimization
//' 1 - method of moments
//' 2 - global optimization
//' 4 - interval optimization
//' interv_step number of intervals to expose 
//' interv_step number of intervals to expose 
//' par - par[0] = b par[1] = a par[2] = m */
//' itv_opt_par - interval optimization parameters */
// [[Rcpp::export]]
Rcpp::List subbofit(
              Rcpp::NumericVector data
             ,int verb = 0
             ,int method = 7
             ,int interv_step = 10 
             ,int output = 0
	     ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
             ,Rcpp::NumericVector par = Rcpp::NumericVector::create(2.,1.,0.) 
             ,Rcpp::NumericVector g_opt_par = Rcpp::NumericVector::create(.1, 1e-2, 100, 1e-3, 1e-5, 3,0)
             ,Rcpp::NumericVector itv_opt_par = Rcpp::NumericVector::create(.01, 1e-3, 200, 1e-3, 1e-5, 5,0) 
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
    { (double)       g_opt_par[0]
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
  
  /* output of initial values */
  /* ------------------------ */
  Rprintf("INITIAL VALUES\n");
  Rprintf("#  par    b     a     m     ll\n");
  Rprintf("#  value  %.2f  %.2f  %.2f  %.2f\n", par[0], par[1], par[2], fmin);
  Rprintf("\n");
  

  /* no minimization */
  /* --------------- */
  if(method == 0){

    // Name of parameters
    Rcpp::CharacterVector param_names = Rcpp::CharacterVector::create("b", "a", "m");

    // main dataframe with coefficients 
      Rcpp::DataFrame dt =
	Rcpp::DataFrame::create( Rcpp::Named("param")     = param_names
				 ,Rcpp::Named("coef")      = par
				 ); 

      Rcpp::List ans = Rcpp::List::create(
					  Rcpp::Named("dt") = dt
					  ,Rcpp::Named("log-likelihood") = fmin
					  );
    
      return ans;

  }

  /* method of moments estimation */
  /* ---------------------------- */
  if(method >= 1){

    // calculate parameters with the method of moments
    Rcpp::List moment_results  = optim_method_moments(data, fmin, provided_m_);
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
		   );

    /* store final values */
    /* ------------------ */
    par[0] = x[0]; //sets b to the new minimum
    par[2] = x[1]; //sets m to the new minimum
    par[1] = geta(data, par[0], par[2]); // updates 'a' value
    
    /* interval optimization  */
    /* ---------------------- */

    if(method >= 4){

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
  RcppGSL::Matrix V =  varcovar(par, Size, dim);
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

