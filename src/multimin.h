/*insert GNU extensions*/
#define _GNU_SOURCE
/*in particular, use of NAN extension*/

/* C++ ----------- */
#include<RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]
/*- -------------- */

/* GSL ----------- */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/* --------------- */

#include "structs.h"

struct multimin_algorithm {

  const gsl_multimin_fdfminimizer_type *Tfdf;
  const gsl_multimin_fminimizer_type *Tf;
  const char *Tname;
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


void multimin(
	      Rcpp::NumericVector data
              ,size_t n
	      ,Rcpp::NumericVector x
              ,double *fun
              ,Rcpp::IntegerVector type 
	      ,Rcpp::NumericVector xmin
	      ,Rcpp::NumericVector xmax
	      ,void (*f)    (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *)
	      ,void (* df)  (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, Rcpp::NumericVector)
	      ,void (* fdf) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *, Rcpp::NumericVector)
	      ,void *fparams
	      ,const struct multimin_params oparams
              );

//void multimin(size_t, std::vector<double>, double *,
//	 const unsigned *,const double *,const double *,
//	 void (*) (const Rcpp::NumericVector, const size_t,const double *,void *,double *),
//	 void (*) (const Rcpp::NumericVector, const size_t,const double *, void *,double *),
//	 void (*) (const Rcpp::NumericVector, const size_t,const double *, void *,double *,double *),
//	 void *,
//	 const struct multimin_params);
