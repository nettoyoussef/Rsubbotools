/*
  sepfit (ver. 1.2.1) -- Fit a skewed exponential power density via
  likelihood maximization

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


/*
  Verbosity levels:
  0 just the final ouput
  1 the results of intermediate steps
  2 internal information on intermediate steps
  3 gory details on intermediate steps
*/

#include "common.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

/* Global Variables */
/* ---------------- */


/* Output Functions */
/*----------------- */


/*
   x   array of parameters
   N     number of observations
   dim   dimension of the matrix: 2 m known; 3 m unknown
   I     the variance-covariance matrix
*/

RcppGSL::Matrix sep_varcovar(
  const Rcpp::NumericVector par
  , const size_t N
  , const size_t dim
){

  size_t i,j;

  // const double m=par[0];
  // const double s=par[1];
  // const double l=par[2];
  // const double a=par[3];

  // main matrix
  RcppGSL::matrix<double> I(dim,dim);

  /* needed for inversion */
  RcppGSL::matrix<double> J(dim,dim);
  gsl_permutation *P = gsl_permutation_alloc (dim);
  int signum;

  /* define the matrix, the order is mu sigma lambda alpha */

  /* mu - mu */
  I(0,0) = 1;

  /* mu - si */
  I(0,1) = 0;
  I(1,0) = gsl_matrix_get(I,0,1);

  /* mu - la */
  I(0,2) = 0;
  I(2,0) = gsl_matrix_get(I,0,2);

  /* mu - al */
  I(0,3) = 0;
  I(3,0) = gsl_matrix_get(I,0,3);

  /* si - si */
  I(1,1) = 1;

  /* si - la */
  I(1,2) = 0;
  I(2,1) = gsl_matrix_get(I,1,2);

  /* si - al */
  I(1,3) = 0;
  I(3,1) = gsl_matrix_get(I,1,3);

  /* la - la */
  I(2,2) = 1;

  /* la - al */
  I(2,3) = 0;
  I(3,2) = gsl_matrix_get(I,2,3);

  /* al - al */
  I(3,3) = 1;


  /* invert I; in J store temporary LU decomposition. */
  gsl_matrix_memcpy (J,I);
  gsl_linalg_LU_decomp (J,P,&signum);
  gsl_linalg_LU_invert (J,P,I);

  /* free allocated memory */
  gsl_permutation_free(P);

  /* set the var-covar matrix */
  for(i=0;i<dim;i++)
      for(j=0;j<i;j++)
        I(i,j) = I(i,j)/sqrt(I(i,i)* I(j,j));

    for(i=0;i<dim;i++)
      for(j=i;j<dim;j++)
        I(i,j) = I(i,j)/N;


    return I;
}


/*-------------------- */
/* Objective Function */
/*-------------------- */
void sep_objf(
  Rcpp::NumericVector data
  , const size_t n
  , Rcpp::NumericVector x
  , void *params
  , double *f
){

  unsigned size = data.size();
  double dtmp1=0;
  size_t i;

  const double m=x[0];
  const double a=x[1];
  const double b=x[2];
  const double lambda=x[3];

  for(i=0;i<size;i++){
    const double z = (data[i]-m)/a;
    const double w = sgn(z)*pow(fabs(z),b/2)*lambda*sqrt(2/b);

    dtmp1 += pow(fabs(z),b)/b - log(gsl_cdf_ugaussian_P(w));
  }

  *f = dtmp1/size +(1/b-1)*log(b)+gsl_sf_lngamma (1/b) + log(a);

}


void sep_objdf(
  Rcpp::NumericVector data
  ,const size_t n
  ,Rcpp::NumericVector x
  ,void *params
  ,Rcpp::NumericVector df
){

  unsigned size = data.size();
  size_t i;

  const double m=x[0];
  const double a=x[1];
  const double b=x[2];
  const double lambda=x[3];

  df[0]=df[1]=df[2]=df[3]=0.0;

  for(i=0;i<size;i++){
    const double z  = (data[i]-m)/a;
    const double sz = sgn(z);
    const double az = fabs(z);
    const double w  = sz*pow(az,b/2)*lambda*sqrt(2/b);

    const double dwdz = pow(az,b/2-1)*lambda*sqrt(b/2);

    const double dwdl = sz*pow(az,b/2)*sqrt(2/b);
    const double dwda = w*(log(az)-1/b)/2;

    const double fionfi = gsl_ran_ugaussian_pdf(w)/gsl_cdf_ugaussian_P(w);

    df[0] += fionfi*dwdz - sz*pow(az,b-1);
    df[1] += fionfi*dwdz*z -pow(az,b);
    df[2] += -fionfi*dwda+pow(az,b)*(b*log(az)-1)/(b*b);
    df[3] += -fionfi*dwdl;
  }

  df[0] = df[0]/(size*a);
  df[1] = (df[1]/size+1)/a;
  df[2] = df[2]/size - (log(b)+b-1+gsl_sf_psi(1/b))/(b*b);
  df[3] = df[3]/size;

}


void sep_objfdf(
  Rcpp::NumericVector data
  , const size_t n
  , Rcpp::NumericVector x
  , void *params
  , double *f
  , Rcpp::NumericVector df
){

  unsigned size = data.size();
  size_t i;

  const double m=x[0];
  const double a=x[1];
  const double b=x[2];
  const double lambda=x[3];

  *f=df[0]=df[1]=df[2]=df[3]=0.0;

  for(i=0;i<size;i++){
    const double z  = (data[i]-m)/a;
    const double sz = sgn(z);
    const double az = fabs(z);
    const double w  = sz*pow(az,b/2)*lambda*sqrt(2/b);

    const double dwdz = pow(az,b/2-1)*lambda*sqrt(b/2);
    const double dwdl = sz*pow(az,b/2)*sqrt(2/b);
    const double dwda = w*(log(az)-1/b)/2;

    const double fionfi = gsl_ran_ugaussian_pdf(w)/gsl_cdf_ugaussian_P(w);

    df[0] += fionfi*dwdz - sz*pow(az,b-1);
    df[1] += fionfi*dwdz*z -pow(az,b);
    df[2] += -fionfi*dwda+pow(az,b)*(b*log(az)-1)/(b*b);
    df[3] += -fionfi*dwdl;

    *f += pow(az,b)/b-log(gsl_cdf_ugaussian_P(w));
  }

  df[0] = df[0]/(size*a);
  df[1] = (df[1]/size+1)/a;
  df[2] = df[2]/size - (log(b)+b-1+gsl_sf_psi(1/b))/(b*b);
  df[3] = df[3]/size;

  *f = (*f)/size +(1/b-1)*log(b)+gsl_sf_lngamma (1/b) + log(a);

}
/*---------------- */


//' Fit a Skewed Exponential Power density via maximum likelihood
//'
//' \code{sepfit} returns the parameters, standard errors. negative
//' log-likelihood and covariance matrix of the skewed power exponential
//' for a sample. The process performs a global minimization over the negative
//' log-likelihood function. See details below.
//'
//' The  SEP is a exponential power distribution controlled
//' by four parameters, with formula:
//' \deqn{ f(x; m, b, a, \lambda) = 2 \Phi(w) e^{-|z|^b/b}/(c)}
//' where:
//' \deqn{z = (x-m)/a}
//' \deqn{w = sign(z) |z|^{(b/2)} \lambda \sqrt{2/b}}
//' \deqn{c = 2 ab^{(1/b)-1} \Gamma(1/b)}
//' with \eqn{\Phi} the cumulative normal distribution with mean zero and
//' variance one.
//' Details on this method are available on the package vignette.
//'
//' @param data (NumericVector) - the sample used to fit the distribution.
//' @param verb (int) - the level of verbosity. Select one of:
//' * 0  just the final result
//' * 1  headings and summary table
//' * 2  intermediate steps results
//' * 3  intermediate steps internals
//' * 4+  details of optim. routine
//' @param par NumericVector - vector containing the initial guess for parameters
//' m (location), a (scale), b (shape), lambda (skewness), respectively.
//' Default values of are c(0, 1, 2, 0), i.e. a normal distribution.
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
//' Default values are c(.1, 1e-2, 100, 1e-3, 1e-5, 2).
//' @return a list containing the following items:
//' * "dt" - dataset containing parameters estimations and standard deviations.
//' * "log-likelihood" - negative log-likelihood value.
//' * "matrix" - the covariance matrix for the parameters.
//'
//' @examples
//' sample_subbo <- rpower(1000, 1, 2)
//' sepfit(sample_subbo)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List sepfit(
  Rcpp::NumericVector data
  ,int verb = 0
  ,Rcpp::NumericVector par = Rcpp::NumericVector::create(0., 1., 2., 0.)
  ,Rcpp::NumericVector g_opt_par = Rcpp::NumericVector::create(.1, 1e-2, 100, 1e-3, 1e-5, 2)
){


  // check arguments
  int check_par_size = par.size();
  if(check_par_size != 4){
    Rcpp::stop("Number of default values for parameters must be equal to four.");
  }

  if(par[1]<=0 || par[2]<=0){
    Rcpp::stop("initial values for a and b should be positive\n");

}

  /* initial values */
  /* -------------- */

  /* store guess */
  Rcpp::IntegerVector type(4);
  Rcpp::NumericVector xmin(4);
  Rcpp::NumericVector xmax(4);

  // define optimization parameters
   struct multimin_params global_oparams =
     {(double)        g_opt_par[0]
      ,(double)       g_opt_par[1]
      ,(unsigned int) g_opt_par[2]
      ,(double)       g_opt_par[3]
      ,(double)       g_opt_par[4]
      ,(unsigned int) g_opt_par[5]
     };

  /* store data */
   unsigned int size = data.size();/*the number of data*/

   // Fisher Matrix or negative log-likelihood?
   // log_likelihood
   double fmin = 0.0;

   /* customized admin. of errors */
   /* --------------------------- */
   gsl_set_error_handler_off ();

   /* sort data */
   /* --------- */
   sortRcpp(data);
   //qsort(data, Size, sizeof(double),sort_by_value);

   /* initial values */
   /* -------------- */
   sep_objf(data, 4, par, NULL, &fmin);

   /* output of initial values */
   /* ------------------------ */
   if(verb >0){
     Rprintf("INITIAL VALUES\n");
     Rprintf("#  par    m     a     b    lambda  ll\n");
     Rprintf("#  value  %.2f  %.2f  %.2f  %.2f  %.2f\n"
             , par[0], par[1], par[2], par[3], fmin);
     Rprintf("\n");
   }

   /* ML global maximization */
   /* ---------------------------------- */
   if(verb >0){
     Rprintf("#--- UNCONSTRAINED OPTIMIZATION\n");
   }

   /* set initial minimization boundaries */
   /* ----------------------------------- */
   xmin.fill(0);
   xmax.fill(0);
   type[0] = 0; // location (m)
   type[1] = 4; // scale (a)
   type[2] = 4; // shape (b)
   type[3] = 0; // skewness (lambda)

   /* perform global minimization */
   /* --------------------------- */
   multimin(
     data              // sample
     ,4                // number of parameters
     ,par              // starting guess / returns parameters value
     ,&fmin            // pointer to update the minimum likelihood
     ,type             // type of transformation of the data
     ,xmin             // minimum values for the parameters
     ,xmax             // maximum values for the parameters
     ,sep_objf         // objective function to minimize
     ,sep_objdf        // df/dx of the objective function
     ,sep_objfdf       // objf + objdf
     ,NULL             // fparams
     ,global_oparams   // parameters for the optmization
     ,verb             // set verbosity level
   );

   if(verb >0){
     Rprintf("#>>> mu=%.3e sigma=%.3e lambda=%.3e alpha=%.3e ll=%e\n",
           par[0],par[1],par[2],par[3],fmin);
   }

   // generate outputs
  // variables

  /* allocate var-covar matrix */
  // The V matrix has on its main diagonal the variance of parameters
  // on the lower diagonal the correlation coefficients
  // on the upper diagonal the covariances
    RcppGSL::Matrix V =  sep_varcovar(par, size, 4);
    // this matrix in its upper diagonal presents the covariances
    // and on its lower diagonal presents the correlation coefficients between the parameters

    // vector of standard errors
    Rcpp::NumericVector std_error = Rcpp::NumericVector::create(
      sqrt(V(0,0))  // m (location)
      ,sqrt(V(1,1)) // a (scale)
      ,sqrt(V(2,2)) // b (shape)
      ,sqrt(V(3,3)) // lambda (skewness)
    );

  // Name of parameters
  Rcpp::CharacterVector param_names = Rcpp::CharacterVector::create("m", "a", "b", "lambda");

  // main dataframe with coefficients
  Rcpp::List dt = Rcpp::DataFrame::create(
    Rcpp::Named("param")      = param_names
    ,Rcpp::Named("coef")      = par
    ,Rcpp::Named("std_error") = std_error
  );

  // convert matrix
  Rcpp::NumericMatrix matrix = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(V));

  // add names for the matrices
  colnames(matrix) = param_names;

  Rcpp::List ans = Rcpp::List::create(
    Rcpp::Named("dt")              = dt
    ,Rcpp::Named("log-likelihood") = fmin
    ,Rcpp::Named("matrix")         = matrix
  );

  return ans;
}
