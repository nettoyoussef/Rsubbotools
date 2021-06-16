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


/* Density parameterization

f(x) = 2 \Phi(w) e^{-|z|^alpha/alpha}/ ( sigma C)

where

z = (x-mu)/sigma
w = sing(z) |z|^(alpha/2) lambda \sqrt{2/alpha}
C = 2 alpha^(1/alpha-1) \Gamma(1/alpha)

with Phi the normal distribution with mean zero and variance one.


x[0] = mu             center
x[1] = sigma          width
x[2] = lambda         skewness
x[3] = alpha          kurtosis


 */

#include "common.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

/* Global Variables */
/* ---------------- */


/* Output Functions */
/*----------------- */

double fSEP(double x, void *params){

  const double *X = (double *) params;

  const double m=X[0];
  const double s=X[1];
  const double l=X[2];
  const double a=X[3];

  const double z = (x-m)/s;
  const double w = (z>0 ? 1 : -1)*pow(fabs(z),a/2)*l*sqrt(2/a);
  const double C = 2*pow(a,1/a-1)*gsl_sf_gamma(1./a);

  return 2*gsl_cdf_ugaussian_P(w)*exp(-pow(fabs(z),a)/a)/(s*C);

}

void sepfit_printcumul(Rcpp::NumericVector data, double x[]){

  int size = data.size();
  size_t i;
  double dtmp1;
  double cumul=0;
  double result,error;

  /* allocate at most 1000 subintervals */
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  /* define the function to be integrated */
  gsl_function F;
  F.function = &fSEP;
  F.params = x;


  gsl_integration_qagil(&F,data[0],1e-7,1e-7,1000,w,&result,&error);

  cumul=result;

  Rprintf("%e %e\n",data[0],cumul);

  for(i=1;i<size;i++){
    gsl_integration_qag(&F,data[i-1],data[i],
                        1e-7,1e-7,1000,4,w,&result,&error);
    cumul+=result;
    Rprintf("%e %e\n",data[i],cumul);
  }

}

void sepfit_printdensity(Rcpp::NumericVector data, double x[]){

  int size = data.size();
  size_t i;

  for(i=0;i<size;i++){
    double dtmp1=data[i];
    Rprintf("%e ",dtmp1);
    Rprintf("%e\n",fSEP(dtmp1,x));
  }

}
/*----------------- */



/*
   x   array of parameters
   N     number of observations
   dim   dimension of the matrix: 2 m known; 3 m unknown
   I     the variance-covariance matrix
*/

RcppGSL::Matrix sep_varcovar(const Rcpp::NumericVector par, const size_t N, const size_t dim){

  size_t i,j;

  const double m=par[0];
  const double s=par[1];
  const double l=par[2];
  const double a=par[3];

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
  I(0,1) = I(1,0) = 0;

  /* mu - la */
  I(0,2) = I(2,0) = 0;

  /* mu - al */
  I(0,3) = I(3,0) = 0;

  /* si - si */
  I(1,1) = 1;

  /* si - la */
  I(1,2) = I(2,1) = 0;

  /* si - al */
  I(1,3) = I(3,1) = 0;

  /* la - la */
  I(2,2) = 1;

  /* la - al */
  I(2,3) = I(3,2) = 0;

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



/* Object Function */
/*---------------- */

void sep_objf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f){

  unsigned size = data.size();
  double dtmp1=0;
  size_t i;

  const double m=x[0];
  const double s=x[1];
  const double l=x[2];
  const double a=x[3];

  for(i=0;i<size;i++){
    const double z = (data[i]-m)/s;
    const double w = (z>0 ? 1 : -1)*pow(fabs(z),a/2)*l*sqrt(2/a);

    dtmp1 += pow(fabs(z),a)/a - log(gsl_cdf_ugaussian_P(w));
  }

  *f = dtmp1/size +(1/a-1)*log(a)+gsl_sf_lngamma (1/a) + log(s);

    /*   fprintf(Fmessages,"# ---- objf mu=%.3e sigma=%.3e lambda=%.3e alpha=%.3e ll=%.3e\n", */
    /*    x[0],x[1],x[2],x[3],*f); */
}


void sep_objdf(
               Rcpp::NumericVector data
               ,const size_t n
               ,Rcpp::NumericVector x
               ,void *params
               ,Rcpp::NumericVector df
               ){

  unsigned size = data.size();
  double dtmp1;

  size_t i;

  const double m=x[0];
  const double s=x[1];
  const double l=x[2];
  const double a=x[3];

  df[0]=df[1]=df[2]=df[3]=0.0;

  for(i=0;i<size;i++){
    const double z = (data[i]-m)/s;
    const double sz = (z>0 ? 1 : -1);
    const double az = fabs(z);
    const double w = sz*pow(az,a/2)*l*sqrt(2/a);

    const double dwdz = pow(az,a/2-1)*l*sqrt(a/2);

    const double dwdl = sz*pow(az,a/2)*sqrt(2/a);
    const double dwda = w*(log(az)-1/a)/2;

    const double fionfi = gsl_ran_ugaussian_pdf(w)/gsl_cdf_ugaussian_P(w);

    df[0] += fionfi*dwdz - sz*pow(az,a-1);
    df[1] += fionfi*dwdz*z -pow(az,a);
    df[2] += -fionfi*dwdl;
    df[3] += -fionfi*dwda+pow(az,a)*(a*log(az)-1)/(a*a);
  }

  df[0] = df[0]/(size*s);
  df[1] = (df[1]/size+1)/s;
  df[2] = df[2]/size;
  df[3] = df[3]/size - (log(a)+a-1+gsl_sf_psi(1/a))/(a*a);

  /*   fprintf(Fmessages,"# ---- objdf mu=%.3e sigma=%.3e lambda=%.3e alpha=%.3e\n", */
  /*      x[0],x[1],x[2],x[3]); */

}



void sep_objfdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f, Rcpp::NumericVector df){

  unsigned size = data.size();
  size_t i;

  const double m=x[0];
  const double s=x[1];
  const double l=x[2];
  const double a=x[3];

  *f=df[0]=df[1]=df[2]=df[3]=0.0;

  for(i=0;i<size;i++){
    const double z = (data[i]-m)/s;
    const double sz = (z>0 ? 1 : -1);
    const double az = fabs(z);
    const double w = sz*pow(az,a/2)*l*sqrt(2/a);

    const double dwdz = pow(az,a/2-1)*l*sqrt(a/2);
    const double dwdl = sz*pow(az,a/2)*sqrt(2/a);
    const double dwda = w*(log(az)-1/a)/2;

    const double fionfi = gsl_ran_ugaussian_pdf(w)/gsl_cdf_ugaussian_P(w);

    df[0] += fionfi*dwdz - sz*pow(az,a-1);
    df[1] += fionfi*dwdz*z -pow(az,a);
    df[2] += -fionfi*dwdl;
    df[3] += -fionfi*dwda+pow(az,a)*(a*log(az)-1)/(a*a);

    *f += pow(az,a)/a-log(gsl_cdf_ugaussian_P(w));
  }

  df[0] = df[0]/(size*s);
  df[1] = (df[1]/size+1)/s;
  df[2] = df[2]/size;
  df[3] = df[3]/size - (log(a)+a-1+gsl_sf_psi(1/a))/(a*a);

  *f = (*f)/size +(1/a-1)*log(a)+gsl_sf_lngamma (1/a) + log(s);


  /*   fprintf(Fmessages,"# ---- objfdf mu=%.3e sigma=%.3e lambda=%.3e alpha=%.3e ll=%.3e grad= %.3e %.3e %.3e %.3e \n", */
  /*      x[0],x[1],x[2],x[3],*f,df[0],df[1],df[2],df[3]); */


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
//' \deqn{ f(x; \mu, \alpha, \lambda, \sigma) =
//' 2 \Phi(w) e^{-|z|^\alpha/\alpha}/ ( \sigma C)}
//' where:
//' \deqn{z = (x-\mu)/\sigma}
//' \deqn{w = sign(z) |z|^{(\alpha/2)} \lambda \sqrt{2/\alpha}}
//' \deqn{C = 2 \alpha^{(1/\alpha-1)} \Gamma(1/\alpha)}
//' with \eqn{\Phi} the cumulative normal distribution with mean zero and variance
//' one.
//' Details on this method are available on the package vignette.
//'
//' @param data (NumericVector) - the sample used to fit the distribution.
//' @param verb (int) - the level of verbosity. Select one of:
//' * 0  just the final result
//' * 1  headings and summary table
//' * 2  intermediate steps results
//' * 3  intermediate steps internals
//' * 4+  details of optim. routine
//' @param par NumericVector - vector containing the initial guess for
//' parameters mu, sigma, lambda and alpha, respectively. Default values of are
//' c(0, 1, 0, 2).
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
                  ,Rcpp::NumericVector par = Rcpp::NumericVector::create(0., 1., 0., 2.)
                  ,Rcpp::NumericVector g_opt_par = Rcpp::NumericVector::create(.1, 1e-2, 100, 1e-3, 1e-5, 2)
                  ){


  // check arguments
  int check_par_size = par.size();
  if(check_par_size != 4){
    Rcpp::stop("Number of default values for parameters must be equal to four.");
  }

  if(par[1]<=0 || par[3]<=0){
    Rcpp::stop("initial values for sigma and alpha should be positive\n");

}

  /* initial values */
  /* -------------- */

  /* store guess */
  Rcpp::IntegerVector type(4);
  Rcpp::NumericVector xmin(4);
  Rcpp::NumericVector xmax(4);

  // define optimization parameters
   struct multimin_params global_oparams =
     {(double)       g_opt_par[0]
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
   double fmin = 0;

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
     Rprintf("#  par    mu     sigma     lambda     alpha     ll\n");
     Rprintf("#  value  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n"
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
   type[0] = 0;
   type[1] = 4;
   type[2] = 0;
   type[3] = 4;

   /* perform global minimization */
   /* --------------------------- */
   multimin(data              // sample
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
  Rcpp::NumericVector std_error =
     Rcpp::NumericVector::create(
                                 sqrt(V(0,0))  // mu
                                 ,sqrt(V(1,1)) // si
                                 ,sqrt(V(2,2)) // la
                                 ,sqrt(V(3,3)) // al
                                 );

  // Name of parameters
  Rcpp::CharacterVector param_names = Rcpp::CharacterVector::create("mu", "si", "la", "al");

  // main dataframe with coefficients
  Rcpp::List dt =
  Rcpp::DataFrame::create(Rcpp::Named("param")      = param_names
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
