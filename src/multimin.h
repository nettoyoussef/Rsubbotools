/*insert GNU extensions*/
#define _GNU_SOURCE
/*in particular, use of NAN extension*/

/* C++ ----------- */
#include<Rcpp.h>
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

void multimin(
              std::vector<double> data
              ,size_t n
              ,std::vector<double> x
              ,double *fun
              ,Rcpp::IntegerVector type
              ,Rcpp::NumericVector xmin
              ,Rcpp::NumericVector xmax
              ,void (*f)    (std::vector<double>, const size_t, std::vector<double>, void *, double *)
              ,void (* df)  (std::vector<double>, const size_t, std::vector<double>, void *, std::vector<double> )
              ,void (* fdf) (std::vector<double>, const size_t, std::vector<double>, void *, double *, std::vector<double> )
              ,void *fparams
              ,const struct multimin_params oparams
              ,int verb
              );

//void multimin(size_t, std::vector<double>, double *,
//       const unsigned *,const double *,const double *,
//       void (*) (const Rcpp::NumericVector, const size_t,const double *,void *,double *),
//       void (*) (const Rcpp::NumericVector, const size_t,const double *, void *,double *),
//       void (*) (const Rcpp::NumericVector, const size_t,const double *, void *,double *,double *),
//       void *,
//       const struct multimin_params);
