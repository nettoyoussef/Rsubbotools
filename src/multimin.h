#ifndef MULTIMIN_H
#define MULTIMIN_H

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
#include "common.h"
#include "structs.h"

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
              ,int verb
              );

//void multimin(size_t, std::vector<double>, double *,
//       const unsigned *,const double *,const double *,
//       void (*) (const Rcpp::NumericVector, const size_t,const double *,void *,double *),
//       void (*) (const Rcpp::NumericVector, const size_t,const double *, void *,double *),
//       void (*) (const Rcpp::NumericVector, const size_t,const double *, void *,double *,double *),
//       void *,
//       const struct multimin_params);

#endif /* MULTIMIN_H */
