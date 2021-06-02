


#include "common.h"

/*------------------*/


/* Output Functions */
/*----------------- */

void subbola_printcumul(Rcpp::NumericVector data, double param[]){

  int size = data.size();
  int i;
  double dtmp1;
  const double bl=param[0];
  const double br=param[1];
  const double a=param[2];
  const double m=param[3];

  const double gl=gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl);
  const double gr=gsl_sf_gamma(1./br+1.)*pow(br,1./br);
  const double gsum=gl+gr;

  for(i=0;i<size;i++){
    if(data[i]>m){
      dtmp1=pow((data[i]-m)/a,br)/br;
      dtmp1=(gl+gr*gsl_sf_gamma_inc_P(1./br,dtmp1))/gsum;
     }
    else{
      dtmp1=pow((m-data[i])/a,bl)/bl;
      dtmp1=(gl-gl*gsl_sf_gamma_inc_P(1./bl,dtmp1))/gsum;
    }
    Rprintf("%e %e\n", data[i], dtmp1);
  }

}

void subbola_printdensity(Rcpp::NumericVector data, double param[]){

  int size = data.size();
  int i;
  const double bl=param[0];
  const double br=param[1];
  const double a=param[2];
  const double m=param[3];

  const double norm=a*(pow(bl,1/bl)*gsl_sf_gamma(1./bl+1.)
                       +pow(br,1/br)*gsl_sf_gamma(1./br+1));

  for(i=0;i<size;i++){
    double dtmp1=data[i];
    Rprintf("%e ",dtmp1);
    dtmp1-=m;
    if(dtmp1>=0){
      Rprintf("%e\n",exp(-pow(dtmp1/a,br)/br)/norm);
    }
    else{
      Rprintf("%e\n",exp(-pow(-dtmp1/a,bl)/bl)/norm);
    }
  }

}
/*----------------- */


/* Object Function */
/*---------------- */

void subbola_objf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f){

  // data size
  int size = data.size();

  double dtmp1;
  unsigned utmp1;

  double sumL=0.0;
  double sumR=0.0;

  const double bl=x[0];
  const double br=x[1];
  const double a=x[2];
  const double mu=x[3];

  gsl_sf_result result; // Do i need to free this?
  int status;

  if( (status = gsl_sf_gamma_e(1./bl+1.,&result)) ){
    Rprintf("Error - bl=%e\n",bl);
    Rcpp::stop("gamma status  = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_gamma_e(1./br+1.,&result)) ){
    Rprintf("Error - br=%e\n",br);
    Rcpp::stop("gamma status  = %s\n", gsl_strerror(status));
  }


  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]>mu) break;
    sumL+=pow(mu-data[utmp1],bl);
  }
  for(;utmp1<size;utmp1++){
    sumR+=pow(data[utmp1]-mu,br);
  }

  sumL /= size;
  sumR /= size;

  dtmp1 = pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)
    +pow(br,1./br)*gsl_sf_gamma(1./br+1.);

  *f = log(a*dtmp1)
    +sumL/(pow(a,bl)*bl)+sumR/(pow(a,br)*br);

}


void subbola_objdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, Rcpp::NumericVector df){

  // data size
  int size = data.size();

  // temp variables
  double dtmp1;

  // deferencing variables to
  // make code more clear
  const double bl=x[0];
  const double br=x[1];
  const double a=x[2];
  const double mu=x[3];

  unsigned utmp1;
  double sumL=0.0;
  double sumR=0.0;
  double sumL1=0.0;
  double sumR1=0.0;
  double sumLl=0.0;
  double sumRl=0.0;

  gsl_sf_result result; // Do i need to free this?
  int status;


  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]<mu){
      sumL1+= (dtmp1 = pow(mu-data[utmp1],bl-1.) );
      dtmp1*= mu-data[utmp1];
      sumL+= dtmp1;
      sumLl+=log(mu-data[utmp1])*dtmp1;
    }
    else if (data[utmp1]>mu){
      sumR1+=(dtmp1 = pow(data[utmp1]-mu,br-1.));
      dtmp1*= data[utmp1]-mu;
      sumR+=dtmp1;
      sumRl+=log(data[utmp1]-mu)*dtmp1;
    }
  }

  sumL /= size;
  sumR /= size;
  sumL1 /= size;
  sumR1 /= size;
  sumLl /= size;
  sumRl /= size;


  if( (status = gsl_sf_gamma_e(1./bl+1.,&result)) ){
    Rprintf("Error - bl=%e\n",bl);
    Rcpp::stop("gamma status  = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_gamma_e(1./br+1.,&result)) ){
    Rprintf("Error - br=%e\n",br);
    Rcpp::stop("gamma status  = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_psi_e(1./bl+1.,&result)) ){
    Rprintf("Error - bl=%e\n",bl);
    Rcpp::stop("psi status = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_psi_e(1./br+1.,&result)) ){
    Rprintf("Error - br=%e\n",br);
    Rcpp::stop("psi status = %s\n", gsl_strerror(status));
  }


  dtmp1 = pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)+pow(br,1./br)*gsl_sf_gamma(1./br+1.);

  df[0] = (1.-log(bl)-gsl_sf_psi(1./bl+1.))*gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl-2.)/dtmp1
    -(1./(bl*bl) + log(a)/bl)*sumL/pow(a,bl) + sumLl/(pow(a,bl)*bl);

  df[1] = (1.-log(br)-gsl_sf_psi(1./br+1.))*gsl_sf_gamma(1./br+1.)*pow(br,1./br-2.)/dtmp1
    -(1./(br*br) + log(a)/br)*sumR/pow(a,br) + sumRl/(pow(a,br)*br);

  df[2] = 1./a - sumL/pow(a,bl+1.) - sumR/pow(a,br+1.);

  df[3] = sumL1/pow(a,bl) - sumR1/pow(a,br);

}


void subbola_objfdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f, Rcpp::NumericVector df){


  // data size
  int size = data.size();

  // temp variables
  double dtmp1;

  // deferencing variables to
  // make code more clear
  const double bl=x[0];
  const double br=x[1];
  const double a=x[2];
  const double mu=x[3];

  unsigned utmp1;
  double sumL=0.0;
  double sumR=0.0;
  double sumL1=0.0;
  double sumR1=0.0;
  double sumLl=0.0;
  double sumRl=0.0;


  gsl_sf_result result;
  int status;



  if( (status = gsl_sf_gamma_e(1./bl+1.,&result)) ){
    Rprintf("Error - bl=%e\n",bl);
    Rcpp::stop("gamma status  = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_gamma_e(1./br+1.,&result)) ){
    Rprintf("Error - br=%e\n",br);
    Rcpp::stop("gamma status  = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_psi_e(1./bl+1.,&result)) ){
    Rprintf("Error - bl=%e\n",bl);
    Rcpp::stop("psi status = %s\n", gsl_strerror(status));
  }

  if( (status = gsl_sf_psi_e(1./br+1.,&result)) ){
    Rprintf("Error - br=%e\n",br);
    Rcpp::stop("psi status = %s\n", gsl_strerror(status));
  }


  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]<mu){
      sumL1+= (dtmp1 = pow(mu-data[utmp1],bl-1.) );
      dtmp1*= mu-data[utmp1];
      sumL+= dtmp1;
      sumLl+=log(mu-data[utmp1])*dtmp1;
    }
    else if (data[utmp1]>mu){
      sumR1+=(dtmp1 = pow(data[utmp1]-mu,br-1.));
      dtmp1*= data[utmp1]-mu;
      sumR+=dtmp1;
      sumRl+=log(data[utmp1]-mu)*dtmp1;
    }
  }

  sumL /= size;
  sumR /= size;
  sumL1 /= size;
  sumR1 /= size;
  sumLl /= size;
  sumRl /= size;

  dtmp1 = pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)+pow(br,1./br)*gsl_sf_gamma(1./br+1.);


  *f = log(a*dtmp1)+sumL/(pow(a,bl)*bl)+sumR/(pow(a,br)*br);

  df[0] = (1.-log(bl)-gsl_sf_psi(1./bl+1.))*gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl-2.)/dtmp1
    -(1./(bl*bl) + log(a)/bl)*sumL/pow(a,bl) + sumLl/(pow(a,bl)*bl);

  df[1] = (1.-log(br)-gsl_sf_psi(1./br+1.))*gsl_sf_gamma(1./br+1.)*pow(br,1./br-2.)/dtmp1
    -(1./(br*br) + log(a)/br)*sumR/pow(a,br) + sumRl/(pow(a,br)*br);

  df[2] = 1./a - sumL/pow(a,bl+1.) - sumR/pow(a,br+1.);

  df[3] = sumL1/pow(a,bl) - sumR1/pow(a,br);

}
/*---------------- */


// subbolafit (ver. 1.2.1) -- Fit a (less) asymmetric power exponential distribution
// Copyright (C) 2003 Giulio Bottazzi
//
//
//
// Verbosity levels:
// 0 just the final ouput
// 1 the results of intermediate steps
// 2 internal information on intermediate steps
// 3 gory details on intermediate steps
//

//' short help*/
//' Fit skewed power exponential density. Read from files or from standard input    \n\n");
//' Usage: %s [options] [files]\n\n",argv[0]);
//'  Options:                                                            \n");
//'  -O  output type (default 0)                \n");
//'       0  parameter bl br a m and log-likelihood   \n");
//'       1  the estimated distribution function computed on the provided points     \n");
//'       2  the estimated density function computed on the provided points \n");
//'  -x  set initial conditions bl,br,a,m  (default 2,2,1,1,0)\n");
//'  -m  the mode is not estimated but is set to the value provided\n");
//'  -s  number of intervals to explore at each iteration (default 10)\n");
//'  -V  verbosity level (default 0)           \n");
//'       0  just the final result        \n");
//'       1  intermediate steps results   \n");
//'       2  intermediate steps internals \n");
//'       3+  details of optim. routine   \n");
//'  -M  active estimation steps. The value is the sum of (default 1)\n");
//'       0  no optimization\n");
//'       1  global and local optimization \n");
//'  -G  set global optimization options. Fields are step,tol,iter,eps,msize,algo.\n");
//'      Empty field implies default (default .1,1e-2,100,1e-3,1e-5,2)\n");
//'  -I  set local optimization options. Fields are step,tol,iter,eps,msize,algo.\n");
//'      Empty field implies default (default .01,1e-3,200,1e-3,1e-5,2)\n");
//' The optimization parameters are");
//'  step  initial step size of the searching algorithm                  \n");
//'  tol  line search tolerance iter: maximum number of iterations      \n");
//'  eps  gradient tolerance : stopping criteria ||gradient||<eps       \n");
//'  msize  simplex max size : stopping criteria ||max edge||<msize     \n");
//'  algo  optimization methods: 0 Fletcher-Reeves, 1 Polak-Ribiere,     \n");
//'        2 Broyden-Fletcher-Goldfarb-Shanno, 3 Steepest descent,           \n");
//'        4 Nelder-Mead simplex, 5 Broyden-Fletcher-Goldfarb-Shanno ver.2   \n");
//' Examples:\n");
//'  'subbolafit -m 1 <file'  estimate bl,br,a with m=1\n");
// [[Rcpp::export]]
Rcpp::List subbolafit(
                      Rcpp::NumericVector data
                     ,int verb = 0
                     ,int method = 7
                     ,int interv_step = 10
                     ,int output = 0
                     ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
                     ,Rcpp::NumericVector par = Rcpp::NumericVector::create(2., 2., 1., 0.)
                     ,Rcpp::NumericVector g_opt_par = Rcpp::NumericVector::create(.1, 1e-2, 100, 1e-3, 1e-5, 3,0)
                     ,Rcpp::NumericVector itv_opt_par = Rcpp::NumericVector::create(.01, 1e-3, 200, 1e-3, 1e-5, 5,0)
                     ){


  // check arguments
  int check_par_size = par.size();
  if(check_par_size != 4){
    Rcpp::stop("Number of default values for parameters must be equal to four.");
  }

  /* initial values */
  /* -------------- */

   // Name of parameters
   Rcpp::CharacterVector param_names = Rcpp::CharacterVector::create("bl", "br", "a", "m");

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

    // casting the null value
    // necessary according to
    // https://gallery.rcpp.org/articles/optional-null-function-arguments/
    Rcpp::NumericVector provided_m(provided_m_);
    par[3] = provided_m[0];


  }

  /* store data */
  //unsigned int Size = data.size();/*the number of data*/

  /* store guess */
  double   fmin=0;

  /* customized admin. of errors */
  /* --------------------------- */
  gsl_set_error_handler_off ();


  /* sort data */
  /* --------- */
  sortRcpp(data);

  /* initial values */
  /* -------------- */
  subbola_objf(data, 4, par, NULL, &fmin);

  /* output of initial values */
  /* ------------------------ */
  if(verb > 0){
    Rprintf("INITIAL VALUES\n");
    Rprintf("INITIAL VALUES\n");
    Rprintf("#  par    bl    br    a     m     ll\n");
    Rprintf("#  value  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n", par[0], par[1], par[2], par[3], fmin);
    Rprintf("\n");
  }

  /* no minimization */
  /* --------------- */
  if(method == 0){

    // main dataframe with coefficients
    Rcpp::DataFrame dt =
      Rcpp::DataFrame::create(Rcpp::Named("param") = param_names
                              ,Rcpp::Named("coef") = par
                             );

    Rcpp::List ans =
      Rcpp::List::create(
                         Rcpp::Named("dt")              = dt
                         ,Rcpp::Named("log-likelihood") = fmin
                         );
    return ans;

  }


  /* global optimization  */
  /* -------------------- */
  if(method >= 2){

    // notice: par is updated by reference
    Rcpp::List g_opt_results =
      global_optim(
                   data
                   ,fmin
                   ,global_oparams
                   ,par
                   ,(unsigned) 4
                   ,(unsigned) 3
                   ,subbola_objf
                   ,subbola_objdf
                   ,subbola_objfdf
                   ,provided_m_
                   ,verb
                   );

    /* interval optimization  */
    /* ---------------------- */

    if(method >= 4){

      g_opt_results =
        interval_optim(
                       data
                       ,g_opt_results["type"]
                       ,g_opt_results["xmin"]
                       ,g_opt_results["xmax"]
                       ,par
                       ,g_opt_results["fmin"]
                       ,interv_oparams /* interval optimization parameters */
                       ,interv_step /* interval optimization step */
                       ,(unsigned) 4
                       ,(unsigned) 3
                       ,subbola_objf
                       ,subbola_objdf
                       ,subbola_objfdf
                       ,verb
                       );

      // updates log-likelihood
      fmin = Rcpp::as<double>(g_opt_results["fmin"]);

    }
  }


  // generate outputs

  // variables
  //const size_t dim=(is_m_provided?4:5); /* set the size of the var-covar matrix */

  /* allocate var-covar matrix */
  // The V matrix has on its main diagonal the variance of parameters
  // on the lower diagonal the correlation coefficients
  // on the upper diagonal the covariances
  //RcppGSL::Matrix V =  subbola_varcovar(par, Size, dim);
  // this matrix in its upper diagonal presents the covariances
  // and on its lower diagonal presents the correlation coefficients between the parameters

  // vector of standard errors
  //Rcpp::NumericVector std_error =
  //  Rcpp::NumericVector::create(
  //                               sqrt(V(0,0)) // bl
  //                              ,sqrt(V(1,1)) // br
  //                              ,sqrt(V(2,2)) // al
  //                              ,sqrt(V(3,3)) // ar
  //                              ,sqrt(V(4,4)) // m
  //                              );

  // main dataframe with coefficients
  Rcpp::List dt =
    Rcpp::DataFrame::create( Rcpp::Named("param")     = param_names
                            ,Rcpp::Named("coef")      = par
                      );
  //,Rcpp::Named("std_error") = std_error

  // convert matrix
  //Rcpp::NumericMatrix matrix = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(V));

  // add names for the matrices
  //colnames(matrix) = param_names;

  Rcpp::List ans = Rcpp::List::create(
                                       Rcpp::Named("dt")             = dt
                                      ,Rcpp::Named("log-likelihood") = fmin
                                     );
                                      //,Rcpp::Named("matrix")         = matrix


  return ans;
}
